import unittest
import numpy as np

from QGL import *
from QGL.drivers import APSPattern
from QGL.PatternUtils import flatten


class APSPatternUtils(unittest.TestCase):
    def setUp(self):
        self.cl = ChannelLibrary(":memory:")
        self.q1gate = Channels.LogicalMarkerChannel(label='q1-gate', channel_db=self.cl.channelDatabase)
        self.q1 = self.cl.new_qubit(label='q1')
        self.q1.gate_chan = self.q1gate
        self.q1.pulse_params['length'] = 30e-9
        self.cl.update_channelDict()

    def test_unroll_loops_simple(self):
        q1 = self.q1

        seqs = [repeat(2, [qwait(), X(q1), Id(q1)]), repeat(2, [qwait(), Y(q1),
                                                                Id(q1)])]
        a, b = APSPattern.unroll_loops(seqs)
        assert (a == seqs)
        assert (b == 2)

    def test_unroll_loops(self):
        q1 = self.q1

        Xw = APSPattern.APSWaveform(Compiler.Waveform(Z(q1)))
        Yw = APSPattern.APSWaveform(Compiler.Waveform(Y(q1)))
        Zw = APSPattern.APSWaveform(Compiler.Waveform(Z(q1)))

        seqs = [repeat(2, [qwait(), Xw, Yw]), repeat(3, [qwait(), Yw,
                                                                Zw])]
        a, b = APSPattern.unroll_loops(seqs)

        seqUnrolled = [qwait(), Xw, Yw] * 2
        assert (a[0] == seqUnrolled)

        seqUnrolled = [qwait(), Yw, Zw] * 3
        assert (a[1] == seqUnrolled)

        assert (b == 0)

    def test_unroll_nested_loops(self):
        q1 = self.q1
        # since single qubit pulse will be modified convert to APSWaveform
        Xw = APSPattern.APSWaveform(Compiler.Waveform(Z(q1)))
        Yw = APSPattern.APSWaveform(Compiler.Waveform(Y(q1)))
        Zw = APSPattern.APSWaveform(Compiler.Waveform(Z(q1)))

        seqs =  [[Xw, Yw] + repeat(2, [Xw, Yw] + repeat(3, [qwait(), Zw, Xw]))]
        unrolled = list(flatten([Xw, Yw] + 2 * [[Xw, Yw] + 3*[qwait(), Zw, Xw]]))

        a, b = APSPattern.unroll_loops(seqs)

        assert (a[0] == unrolled)

        assert (b == 0)

    def test_unroll_single_entry(self):
        q1 = self.q1
        # since single qubit pulse will be modified convert to APSWaveform
        Xw = APSPattern.APSWaveform(Compiler.Waveform(Z(q1)))
        Yw = APSPattern.APSWaveform(Compiler.Waveform(Y(q1)))
        Yw_repeat = APSPattern.APSWaveform(Compiler.Waveform(Y(q1)))
        Yw_repeat.repeat = 5

        seqs = [repeat(5, [Yw]) + [Xw]]
        a, b = APSPattern.unroll_loops(seqs)

        seqUnrolled = [Yw_repeat, Xw]

        assert (a[0] == seqUnrolled)
        assert (b == 0)


if __name__ == "__main__":
    unittest.main()
