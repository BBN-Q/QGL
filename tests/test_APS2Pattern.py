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

    def test_inplace_updates(self):
        q1 = self.q1
        APS2Pattern.SAVE_WF_OFFSETS = True

        mf = RabiAmp(q1, np.linspace(-1, 1, 11))
        aps2_f = os.path.join(os.path.dirname(mf), "Rabi-BBNAPS1.aps2")

        instructions = APS2Pattern.read_instructions(aps2_f)
        waveforms    = APS2Pattern.read_waveforms(aps2_f)

        t_instructions =
        t_waveforms    =

        # Assert that we can read the correct information from file
        for actual, expected in zip(instructions, t_instructions):
            instrOpCode = (actual.header >> 4) & 0xf
            assert (instrOpCode == expected)

        for actual, expected in zip(waveforms, t_waveforms):
            instrOpCode = (actual.header >> 4) & 0xf
            assert (instrOpCode == expected)

        # overwrite the instructions and waveforms
        # here we completely change the experiment from Rabi to
        # SPAM characterization
        spam_mf = SPAM(q1, np.linspace(-1.0, 1.0, 11));

        offset_f = os.path.join(os.path.dirname(spam_mf), "SPAM-BBNAPS1.offsets"))
        with open(offset_f, "rb") as FID:
            offsets = pickle.load(FID)

        #pulses = {l: Utheta(q1, amp=0.5, phase=0) for l in offsets}
        aps2_f = os.path.dirname(mf), "Rabi-BBNAPS1.aps2")
        wfm_f  = os.path.dirname(mf), "SPAM-BBNAPS1.aps2")
        spam_waveforms = APS2Pattern.read_waveforms(wfm_f)
        # The spam_waveforms are raw.  Is this what we want?
        QGL.drivers.APS2Pattern.update_wf_library(aps2_f, wfm_f, offsets)

        spam_instrs = APS2Pattern.read_instructions(wfm_f)
        APS2Pattern.replace_instructions(aps2_f, spam_instrs)

        # assert the data now in the file is what we wrote above
        instructions = APS2Pattern.read_instructions(aps2_f)
        waveforms    = APS2Pattern.read_waveforms(aps2_f)
        for actual, expected in zip(instructions, spam_instrs):
            instrOpCode = (actual.header >> 4) & 0xf
            assert (instrOpCode == expected)

        for actual, expected in zip(waveforms, spam_waveforms):
            instrOpCode = (actual.header >> 4) & 0xf
            assert (instrOpCode == expected)


if __name__ == "__main__":
    unittest.main()
