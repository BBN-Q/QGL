import h5py
import unittest
import numpy as np

from QGL import *


class CompileUtils(unittest.TestCase):
    def setUp(self):
        self.q1gate = Channels.LogicalMarkerChannel(label='q1-gate')
        self.q1 = Qubit(label='q1', gate_chan=self.q1gate)
        self.q1.pulse_params['length'] = 30e-9

        self.q2gate = Channels.LogicalMarkerChannel(label='q2-gate')
        self.q2 = Qubit(label='q2', gate_chan=self.q2gate)
        self.q2.pulse_params['length'] = 30e-9

        self.trigger = Channels.LogicalMarkerChannel(label='trigger')
        self.measq1 = Channels.Measurement(label='M-q1', meas_type='autodyne')
        self.measq1.trig_chan = self.trigger

        ChannelLibrary(library_file=None) # Create a blank ChannelLibrary
        ChannelLibraries.channelLib.channelDict = {'q1': self.q1,
                                                 'q2': self.q2,
                                                 'M-q1': self.measq1}
        ChannelLibraries.channelLib.build_connectivity_graph()

    def test_add_digitizer_trigger(self):
        q1 = self.q1
        seq = [X90(q1), MEAS(q1), Y(q1), MEAS(q1)]

        PatternUtils.add_digitizer_trigger([seq])
        assert (self.trigger in seq[1].pulses.keys())
        assert (self.trigger in seq[3].pulses.keys())

    def test_add_gate_pulses(self):
        q1 = self.q1
        seq = [X90(q1), Y90(q1)]
        PatternUtils.add_gate_pulses(seq)
        assert ([self.q1gate in entry.pulses.keys() for entry in seq] ==
                [True, True])

        q2 = self.q2
        seq = [X90(q1), X90(q2), X(q1) * Y(q2)]
        PatternUtils.add_gate_pulses(seq)
        assert ([self.q1gate in entry.pulses.keys() for entry in seq] ==
                [True, False, True])
        assert ([self.q2gate in entry.pulses.keys() for entry in seq] ==
                [False, True, True])

    def test_add_slave_trigger(self):
        q1 = self.q1
        trigger = self.trigger
        label = BlockLabel.newlabel()
        seq1 = [qwait(), label, X90(q1)]
        seq2 = [qwait(), X90(q1)]

        PatternUtils.add_slave_trigger([seq1], trigger)
        t = TAPulse("TRIG", trigger, trigger.pulse_params['length'], 1.0, 0.0,
                    0.0)
        assert (seq1 == [qwait(), t, label, X90(q1)])

        PatternUtils.add_slave_trigger([seq2], trigger)
        assert (seq2 == [qwait(), X90(q1) * t])

    def test_concatenate_entries(self):
        q1 = self.q1
        seq = [X90(q1, length=20e-9), Y90(q1, length=40e-9)]
        ll = Compiler.compile_sequence(seq)
        entry = Compiler.concatenate_entries(ll[q1][0], ll[q1][1])
        assert entry.length == seq[0].length + seq[1].length
        wf = np.hstack(
            (seq[0].amp * seq[0].shape, 1j * seq[1].amp * seq[1].shape))
        assert all(abs(entry.shape - wf) < 1e-16)

    def test_pull_uniform_entries(self):
        q1 = self.q1
        q1.pulse_params['length'] = 20e-9
        q2 = self.q2
        q2.pulse_params['length'] = 60e-9
        seq = [(X90(q1) + Y90(q1) + X90(q1)) * X(q2)]
        ll = Compiler.compile_sequence(seq)
        entryIterators = [iter(ll[q1]), iter(ll[q2])]
        entries = [next(e) for e in entryIterators]
        entries, max_length = Compiler.pull_uniform_entries(entries, entryIterators)
        self.assertAlmostEqual(max_length, 60e-9)
        assert all(e.length == max_length for e in entries)
        self.assertRaises(StopIteration, next, entryIterators[0])

        q2.pulse_params['length'] = 40e-9
        seq = [(X90(q1) + Z90(q1) + X90(q1)) * Y(q2)]
        ll = Compiler.compile_sequence(seq)
        entryIterators = [iter(ll[q1]), iter(ll[q2])]
        entries = [next(e) for e in entryIterators]
        entries, max_length = Compiler.pull_uniform_entries(entries, entryIterators)
        self.assertAlmostEqual(max_length, 40e-9)
        assert all(e.length == max_length for e in entries)

    def test_pull_uniform_entries2(self):
        q1 = self.q1
        q1.pulse_params['length'] = 30e-9
        q2 = self.q2
        q2.pulse_params['length'] = 40e-9
        seq = [(X90(q1) + Y90(q1) + X(q1) + Y(q1)) * (Y(q2) + X(q2) + Y(q2))]
        ll = Compiler.compile_sequence(seq)
        entryIterators = [iter(ll[q1]), iter(ll[q2])]
        entries = [next(e) for e in entryIterators]
        entries, max_length = Compiler.pull_uniform_entries(entries, entryIterators)
        self.assertAlmostEqual(max_length, 120e-9)
        self.assertTrue(all(e.length == max_length for e in entries))

    def test_merge_channels(self):
        q1 = self.q1
        q1.pulse_params['length'] = 20e-9
        q2 = self.q2
        q2.pulse_params['length'] = 60e-9
        seqs = [[(X90(q1) + Y90(q1) + X90(q1)) * X(q2)]]
        ll = Compiler.compile_sequences(seqs)

        chLL = Compiler.merge_channels(ll, [q1, q2])
        assert len(chLL[0]) == len(ll[q1][0]) - 2
        assert len(chLL[0]) == len(ll[q2][0]) - 1

    def test_num_measurements(self):
        q1 = self.q1
        q2 = self.q2
        mq1 = MEAS(q1).channel
        mq2 = MEAS(q2).channel
        seqs = [[X(q1)*X(q2)]]
        wireSeqs = Compiler.compile_sequences(seqs)
        assert Compiler.count_measurements(wireSeqs) == 0

        seqs = [[MEAS(q1)]]
        wireSeqs = Compiler.compile_sequences(seqs)
        assert Compiler.count_measurements(wireSeqs) == 1

        seqs = [[MEAS(q1)*MEAS(q2)]]
        wireSeqs = Compiler.compile_sequences(seqs)
        assert Compiler.count_measurements(wireSeqs) == 1

        seqs = [[MEAS(q1), MEAS(q2)]]
        wireSeqs = Compiler.compile_sequences(seqs)
        assert Compiler.count_measurements(wireSeqs) == 1

        seqs = [[MEAS(q1)*MEAS(q2), MEAS(q2)]]
        wireSeqs = Compiler.compile_sequences(seqs)
        assert Compiler.count_measurements(wireSeqs) == 2
        wire_meas = Compiler.count_measurements_per_wire(wireSeqs)
        assert wire_meas[mq1] == 1
        assert wire_meas[mq2] == 2

        seqs = [[MEAS(q1)*MEAS(q2), MEAS(q2)],
                [MEAS(q1)],
                [MEAS(q2)],
                [MEAS(q1), MEAS(q2)],
                [X(q1)]]
        wireSeqs = Compiler.compile_sequences(seqs)
        assert Compiler.count_measurements(wireSeqs) == 5
        wire_meas = Compiler.count_measurements_per_wire(wireSeqs)
        assert wire_meas[mq1] == 3
        assert wire_meas[mq2] == 4

    def test_frame_update(self):
        # test that the compiler replaces Z's with frame updates
        q1 = self.q1
        # sequence of five X's separated by Z90s
        seq = [X(q1), Z90(q1), X(q1), Z90(q1), X(q1), Z90(q1), X(q1), Z90(q1), X(q1)]
        # each X should pick up the frame change of the following Z90, except for
        # the last one
        out_seq = Compiler.compile_sequence(seq)[q1]
        expected_frame_change = [-0.5*np.pi]*4 + [0.0]

        assert len(out_seq) == 5
        for p, frame_change in zip(out_seq, expected_frame_change):
            assert p.frameChange == frame_change

if __name__ == "__main__":
    unittest.main()
