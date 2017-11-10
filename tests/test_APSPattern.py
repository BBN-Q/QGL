import h5py
import unittest
import numpy as np

from QGL import *
from QGL.drivers import APSPattern


class APSPatternUtils(unittest.TestCase):
    def setUp(self):
        # self.q1gate = Channels.LogicalMarkerChannel(label='q1-gate')
        # self.q1 = Qubit(label='q1', gate_chan=self.q1gate)
        self.q1 = Qubit(label='q1')
        self.q1.pulse_params['length'] = 30e-9

        ChannelLibrary(library_file=None) # Create a blank ChannelLibrary
        ChannelLibraries.channelLib.channelDict = {'q1': self.q1}
        ChannelLibraries.channelLib.build_connectivity_graph()

    def test_unroll_loops_simple(self):
        q1 = self.q1
        seqs = [repeat(2, [qwait(), X(q1), Id(q1)]), repeat(2, [qwait(), Y(q1),
                                                                Id(q1)])]
        a, b = APSPattern.unroll_loops(seqs)
        assert (a == seqs)
        assert (b == 2)

    def test_unroll_loops(self):
        q1 = self.q1
        seqs = [repeat(2, [qwait(), X(q1), Id(q1)]), repeat(3, [qwait(), Y(q1),
                                                                Id(q1)])]
        a, b = APSPattern.unroll_loops(seqs)

        seqUnrolled = [qwait(), X(q1), Id(q1)] * 2
        assert (a[0] == seqUnrolled)

        seqUnrolled = [qwait(), Y(q1), Id(q1)] * 3
        assert (a[1] == seqUnrolled)

        assert (b == 0)

    def test_unroll_nested_loops(self):
        q1 = self.q1
        # since single qubit pulse will be modified convert to APSWaveform
        Zq1 = APSPattern.APSWaveform(Compiler.Waveform(Z(q1)))
        seqs = [repeat(2, [X(q1), Y(q1)] + repeat(3, [Zq1]) + [Y(q1), X(
            q1)]), [X(q1), Y(q1)]]
        a, b = APSPattern.unroll_loops(seqs)

        Zq1_looped = APSPattern.APSWaveform(Compiler.Waveform(Z(q1)))
        Zq1_looped.repeat = 3
        seqUnrolled = ([X(q1), Y(q1), Zq1_looped, Y(q1), X(q1)]) * 2

        assert (a[0] == seqUnrolled)

        assert (b == 0)

    def test_unroll_single_entry(self):
        q1 = self.q1
        # since single qubit pulse will be modified convert to APSWaveform
        Xq1 = APSPattern.APSWaveform(Compiler.Waveform(Z(q1)))
        seqs = [repeat(5, [Xq1]) + [Y(q1)]]
        a, b = APSPattern.unroll_loops(seqs)

        Xq1_looped = APSPattern.APSWaveform(Compiler.Waveform(Z(q1)))
        Xq1_looped.repeat = 5
        seqUnrolled = [Xq1_looped, Y(q1)]

        assert (a[0] == seqUnrolled)
        assert (b == 0)


if __name__ == "__main__":
    unittest.main()
