import unittest
import os
import pickle
import numpy as np
from copy import copy

from QGL import *
from QGL.drivers import APS2Pattern


class APSPatternUtils(unittest.TestCase):
    def setUp(self):
        self.cl = ChannelLibrary(db_resource_name=":memory:")
        #self.q1gate = Channels.LogicalMarkerChannel(label='q1-gate',
        #                                    channel_db=self.cl.channelDatabase)
        self.q1 = self.cl.new_qubit(label='q1')
        #self.q1.gate_chan = self.q1gate
        self.q1.pulse_params['length'] = 30e-9

        ip_addresses = [f"192.168.1.{i}" for i in [23, 24, 25, 28]]
        aps2 = self.cl.new_APS2_rack("Maxwell",
                                     ip_addresses,
                                     tdm_ip="192.168.1.11")
        aps2.px("TDM").trigger_interval = 500e-6
        self.cl.set_master(aps2.px("TDM"))

        # initialize all four APS2 to linear regime
        for i in range(1,4):
            aps2.tx(i).ch(1).I_channel_amp_factor = 0.5
            aps2.tx(i).ch(1).Q_channel_amp_factor = 0.5
            aps2.tx(i).ch(1).amp_factor = 1

        dig_1  = self.cl.new_X6("MyX6", address=0)

        dig_1.record_length = 1024 + 256

        #### QUBIT 1 ######################################################
        #### Qubit 1 Instruments ##########################################
        AM1 = self.cl.new_source("AutodyneM1",
                                 "HolzworthHS9000",
                                 "HS9004A-492-1",
                                 power=16.0,
                                 frequency= 6.74621e9, reference="10MHz")

        q1src = self.cl.new_source("q1source",
                                   "HolzworthHS9000",
                                   "HS9004A-492-2",
                                   power=16.0,
                                   frequency=5.0122e9, reference="10MHz")

        self.cl.set_measure(self.q1, aps2.tx(2), dig_1.channels[1],
                                            gate=False,
                                            trig_channel=aps2.tx(2).ch("m2"),
                                            generator=AM1)
        self.cl.set_control(self.q1, aps2.tx(4), generator=q1src)

        #### Qubit 2 Measure Chan ########################################

        self.cl["q1"].measure_chan.autodyne_freq = 11e6
        self.cl["q1"].measure_chan.pulse_params = {"length": 1.0e-6,
                                                   "amp": 0.5,
                                                   "sigma": 1.0e-8,
                                                   "shape_fun": "tanh"}

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
        #print(q1)
        #print(q1.chan)

        # amps = np.linspace(-1, 1, 11)
        # seqs = [[Utheta(q1, amp=amp), MEAS(q1)] for amp in amps]
        #
        # axis_descriptor = [{
        #     'name': 'amplitude',
        #     'unit': None,
        #     'points': list(amps),
        #     'partition': 1
        # }]
        #
        # mf = compile_to_hardware(seqs, 'Rabi/Rabi',
        #                                axis_descriptor=axis_descriptor)

        mf = RabiAmp(q1, np.linspace(-1, 1, 11))
        aps2_f = os.path.join(os.path.dirname(mf), "Rabi-Maxwell_U4.aps2")

        # instructions = APS2Pattern.read_instructions(aps2_f)
        # waveforms    = APS2Pattern.read_waveforms(aps2_f)
        #
        # t_instructions =
        # t_waveforms    =
        #
        # # Assert that we can read the correct information from file
        # for actual, expected in zip(instructions, t_instructions):
        #      instrOpCode = (actual.header >> 4) & 0xf
        #      assert (instrOpCode == expected)
        #
        # for actual, expected in zip(waveforms, t_waveforms):
        #     instrOpCode = (actual.header >> 4) & 0xf
        #     assert (instrOpCode == expected)

        # overwrite the instructions and waveforms
        # here we completely change the experiment from Rabi to
        # SPAM characterization
        spam_mf = SPAM(q1, np.linspace(-1.0, 1.0, 11));

        offset_f = os.path.join(os.path.dirname(spam_mf), "SPAM-Maxwell_U4.offsets")
        with open(offset_f, "rb") as FID:
            offsets = pickle.load(FID)

        #pulses = {l: Utheta(q1, amp=0.5, phase=0) for l in offsets}
        aps2_f = os.path.join(os.path.dirname(mf), "Rabi-Maxwell_U4.aps2")
        wfm_f  = os.path.join(os.path.dirname(spam_mf), "SPAM-Maxwell_U4.aps2")
        spam_waveforms = APS2Pattern.read_waveforms(wfm_f)
        print(spam_waveforms)
        # The spam_waveforms are raw.  Is this what we want?
        APS2Pattern.update_wf_library(aps2_f, spam_waveforms, offsets)

        spam_instrs = APS2Pattern.read_instructions(wfm_f)
        APS2Pattern.replace_instructions(aps2_f, spam_instrs)

        # assert the data now in the file is what we wrote above
        instructions = APS2Pattern.read_instructions(wfm_f)
        waveforms    = APS2Pattern.read_waveforms(wfm_f)
        for actual, expected in zip(instructions, spam_instrs):
            instrOpCode = (actual.header >> 4) & 0xf
            assert (instrOpCode == expected)

        for actual, expected in zip(waveforms, spam_waveforms):
            instrOpCode = (actual.header >> 4) & 0xf
            assert (instrOpCode == expected)


if __name__ == "__main__":
    unittest.main()
