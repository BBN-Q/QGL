import h5py
import unittest
import numpy as np
from copy import copy

from QGL import *
from QGL.drivers import APS2Pattern


class APSPatternUtils(unittest.TestCase):
    def setUp(self):
        self.cl = ChannelLibrary(db_resource_name=":memory:")
        self.q1gate = Channels.LogicalMarkerChannel(label='q1-gate', channel_db=self.cl.channelDatabase)
        self.q1 = self.cl.new_qubit(label='q1')
        self.q1.gate_chan = self.q1gate
        self.q1.pulse_params['length'] = 30e-9
        self.cl.update_channelDict()

    def test_synchronize_control_flow(self):
        q1 = self.q1

        pulse = Compiler.Waveform()
        pulse.length = 24
        pulse.key = 12345
        delay = Compiler.Waveform()
        delay.length = 100
        delay.isTimeAmp = True
        blank = Compiler.Waveform(BLANK(q1, pulse.length))

        seq_1 = [qwait(), delay, copy(pulse), qwait(), copy(pulse)]
        seq_2 = [qwait(), copy(blank), qwait(), copy(blank)]
        offsets = {APS2Pattern.wf_sig(pulse): 0}

        instructions = APS2Pattern.create_seq_instructions(
            [seq_1, [], seq_2, [], [], []], offsets)[0]

        instr_types = [
            APS2Pattern.SYNC, APS2Pattern.WAIT, APS2Pattern.WFM,
            APS2Pattern.MARKER, APS2Pattern.WFM, APS2Pattern.WAIT,
            APS2Pattern.WFM, APS2Pattern.MARKER
        ]

        for actual, expected in zip(instructions, instr_types):
            instrOpCode = (actual.header >> 4) & 0xf
            assert (instrOpCode == expected)


if __name__ == "__main__":
    unittest.main()
