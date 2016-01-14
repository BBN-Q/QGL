import h5py
import unittest
import numpy as np

from QGL import *

class CompileUtils(unittest.TestCase):
    def setUp(self):
        self.q1gate = Channels.LogicalMarkerChannel(label='q1-gate')
        self.q1 = Qubit(label='q1', gateChan=self.q1gate)
        self.q1.pulseParams['length'] = 30e-9

        self.q2gate = Channels.LogicalMarkerChannel(label='q2-gate')
        self.q2 = Qubit(label='q2', gateChan=self.q2gate)
        self.q2.pulseParams['length'] = 30e-9

        self.trigger = Channels.LogicalMarkerChannel(label='trigger')
        self.measq1 = Channels.Measurement(label='M-q1')
        self.measq1.trigChan = self.trigger

        Compiler.channelLib = {'q1': self.q1, 'q2': self.q2, 'M-q1': self.measq1}

    def test_add_digitizer_trigger(self):
        q1 = self.q1
        seq = [X90(q1), MEAS(q1), Y(q1), MEAS(q1)]

        PatternUtils.add_digitizer_trigger([seq])
        assert(self.trigger in seq[1].pulses.keys())
        assert(self.trigger in seq[3].pulses.keys())

    def test_add_gate_pulses(self):
        q1 = self.q1
        seq = [X90(q1), Y90(q1)]
        PatternUtils.add_gate_pulses([seq])
        assert([self.q1gate in entry.pulses.keys() for entry in seq] == [True, True])

        q2 = self.q2
        seq = [X90(q1), X90(q2), X(q1)*Y(q2)]
        PatternUtils.add_gate_pulses([seq])
        assert([self.q1gate in entry.pulses.keys() for entry in seq] == [True, False, True])
        assert([self.q2gate in entry.pulses.keys() for entry in seq] == [False, True, True])

    def test_concatenate_entries(self):
        q1 = Qubit(label='q1')
        seq = [X90(q1, length=20e-9), Y90(q1, length=40e-9)]
        ll = Compiler.compile_sequence(seq)
        entry = Compiler.concatenate_entries(ll[q1][0], ll[q1][1])
        assert entry.length == seq[0].length + seq[1].length
        wf = np.hstack((seq[0].amp*seq[0].shape, 1j*seq[1].amp*seq[1].shape))
        assert all(abs(entry.shape - wf) < 1e-16)

    def test_pull_uniform_entries(self):
        q1 = Qubit(label='q1')
        q1.pulseParams['length'] = 20e-9
        q2 = Qubit(label='q2')
        q2.pulseParams['length'] = 60e-9
        seq = [(X90(q1)+Y90(q1)+X90(q1)) * X(q2)]
        ll = Compiler.compile_sequence(seq)
        entryIterators = [iter(ll[q1]), iter(ll[q2])]
        entries = [e.next() for e in entryIterators]
        blocklen = Compiler.pull_uniform_entries(entries, entryIterators)
        self.assertAlmostEqual(blocklen, 60e-9)
        assert all(e.length == blocklen for e in entries)
        self.assertRaises(StopIteration, entryIterators[0].next)

        q2.pulseParams['length'] = 40e-9
        seq = [(X90(q1) + Z90(q1) + X90(q1)) * Y(q2)]
        ll = Compiler.compile_sequence(seq)
        entryIterators = [iter(ll[q1]), iter(ll[q2])]
        entries = [e.next() for e in entryIterators]
        blocklen = Compiler.pull_uniform_entries(entries, entryIterators)
        self.assertAlmostEqual(blocklen, 40e-9)
        assert all(e.length == blocklen for e in entries)

    def test_pull_uniform_entries2(self):
        q1 = Qubit(label='q1')
        q1.pulseParams['length'] = 30e-9
        q2 = Qubit(label='q2')
        q2.pulseParams['length'] = 40e-9
        seq = [(X90(q1) + Y90(q1) + X(q1) + Y(q1)) * (Y(q2) + X(q2) + Y(q2))]
        ll = Compiler.compile_sequence(seq)
        entryIterators = [iter(ll[q1]), iter(ll[q2])]
        entries = [e.next() for e in entryIterators]
        blocklen = Compiler.pull_uniform_entries(entries, entryIterators)
        self.assertAlmostEqual(blocklen, 120e-9)
        self.assertTrue( all(e.length == blocklen for e in entries) )

    def test_merge_channels(self):
        q1 = Qubit(label='q1')
        q1.pulseParams['length'] = 20e-9
        q2 = Qubit(label='q2')
        q2.pulseParams['length'] = 60e-9
        seqs = [[(X90(q1)+Y90(q1)+X90(q1)) * X(q2)]]
        ll = Compiler.compile_sequences(seqs)

        chLL = Compiler.merge_channels(ll, [q1, q2])
        assert len(chLL[0]) == len(ll[q1][0]) - 2
        assert len(chLL[0]) == len(ll[q2][0]) - 1

if __name__ == "__main__":
    unittest.main()
