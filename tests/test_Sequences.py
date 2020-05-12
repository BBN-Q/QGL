import numpy as np
import unittest, time, os, random, sys
import tempfile
import pickle

from QGL import *
from QGL import GSTTools
import QGL

from QGL.Channels import Edge, Measurement, LogicalChannel, LogicalMarkerChannel, PhysicalMarkerChannel, PhysicalQuadratureChannel
from QGL.drivers import APSPattern, APS2Pattern

# Pulled in logger to help debug stand-alone run issue with the config.AWGDir
# (the configuration was NOT gettng loaded when run independently)
#
import logging
logger = logging.getLogger( 'sequences')

import os
istravis = os.environ.get('TRAVIS') == 'true'

#Determine if tests are running in Travis

class AWGTestHelper(object):
    testFileDirectory = './tests/test_data/awg/'

    def __init__(self, translatorModule=None, tolerance=2.0 / 2**13):
        self.channels = {}
        self.assign_channels()
        self.set_awg_dir()
        if translatorModule:
            self.read_function = translatorModule.read_sequence_file
        else:
            self.read_function = None
        self.tolerance = tolerance

    def finalize_map(self, mapping):
        for name, value in mapping.items():
            self.channels[name].phys_chan = self.channels[value]

        self.cl = ChannelLibrary(db_resource_name=":memory:")
        self.cl.clear()
        self.cl.session.add_all(self.channels.values())
        for chan in self.channels.values():
            chan.channel_db = self.cl.channelDatabase
        self.cl.update_channelDict()
        (self.q1, self.q2) = self.get_qubits()

    def assign_channels(self):

        self.qubit_names = ['q1', 'q2']
        self.logical_names = ['digitizerTrig', 'slave_trig']

        self.assign_logical_channels()

    def assign_logical_channels(self):

        for name in self.logical_names:
            self.channels[name] = LogicalMarkerChannel(label=name)

        for name in self.qubit_names:
            mName = 'M-' + name
            mgName = 'M-' + name + '-gate'
            qgName = name + '-gate'

            mg = LogicalMarkerChannel(label=mgName)
            qg = LogicalMarkerChannel(label=qgName)

            m = Measurement(label=mName,
                            gate_chan=mg,
                            trig_chan=self.channels['digitizerTrig'],
                            meas_type='autodyne')

            q = Qubit(label=name, gate_chan=qg)
            q.pulse_params['length'] = 30e-9
            q.pulse_params['phase'] = pi / 2

            self.channels[name] = q
            self.channels[mName] = m
            self.channels[mgName] = mg
            self.channels[qgName] = qg

        # this block depends on the existence of q1 and q2
        self.channels['cr-gate'] = LogicalMarkerChannel(label='cr-gate')

        q1, q2 = self.channels['q1'], self.channels['q2']
        cr = Edge(label="cr",
                  source=q1,
                  target=q2,
                  gate_chan=self.channels['cr-gate'])
        cr.pulse_params['length'] = 30e-9
        cr.pulse_params['phase'] = pi / 4
        self.channels["cr"] = cr

        mq1q2g = LogicalMarkerChannel(label='M-q1q2-gate')
        self.channels['M-q1q2-gate'] = mq1q2g
        self.channels['M-q1q2'] = Measurement(
            label='M-q1q2',
            gate_chan=mq1q2g,
            trig_chan=self.channels['digitizerTrig'],
            meas_type='autodyne')

    def get_qubits(self):
        return [QGL.ChannelLibraries.channelLib[name]
                for name in self.qubit_names]

    def set_awg_dir(self, footer=""):
        if QGL.config.AWGDir is None:
            QGL.config.load_config()

        if QGL.config.AWGDir is None:
            self.temp_dir = tempfile.TemporaryDirectory()
            QGL.config.AWGDir = self.temp_dir.name
            logger.warning(f"Creating temporary AWG dir at {QGL.config.AWGDir}")

        if not hasattr(self, 'original_awg_dir'):
          self.original_awg_dir = QGL.config.AWGDir

        self.awg_dir = os.path.abspath(self.original_awg_dir + os.path.sep + self.__class__.__name__)
        self.truth_dir = os.path.abspath(self.testFileDirectory + os.path.sep + self.__class__.__name__)

        if footer != "":
            self.awg_dir = self.awg_dir + os.path.sep + footer
            self.truth_dir = self.truth_dir + os.path.sep + footer

        if not os.path.isdir(self.awg_dir):
            os.makedirs(self.awg_dir)
        QGL.config.AWGDir = self.awg_dir

    def compare_sequences(self, seqDir):
        if not self.read_function:
            print("AWGTestHelper.read_function is not defined")
            return

        searchDirectory = os.path.join(self.awg_dir, seqDir)
        truthDirectory = os.path.join(self.truth_dir, seqDir)

        filenames = os.listdir(truthDirectory)

        for filename in filenames:
            truthFile = os.path.join(truthDirectory, filename)
            testFile = os.path.join(searchDirectory, filename)

            self.assertTrue(
                os.path.isfile(truthFile),
                "Truth Data File: {0} not found.".format(truthFile))
            self.compare_file_data(testFile, truthFile)

    def compare_file_data(self, testFile, truthFile):
        awgData = self.read_function(testFile)
        truthData = self.read_function(truthFile)

        awgDataLen = len(awgData)
        truthDataLen = len(truthData)

        self.assertTrue(awgDataLen == truthDataLen,
                        "Expected {0} sequences in file. Found {1}.".format(
                            truthDataLen, awgDataLen))

        for name in truthData:
            self.assertTrue(
                name in awgData,
                "Expected channel {0} not found in file {1}".format(name,
                                                                    testFile))

            for x in range(len(truthData[name])):
                seqA = np.array(truthData[name][x])
                seqB = np.array(awgData[name][x])
                self.compare_sequence(
                    seqA, seqB,
                    "\nFile {0} =>\nChannel {1} Sequence {2}".format(testFile,
                                                                     name, x))

    def compare_sequence(self, seqA, seqB, errorHeader):
        #unroll the time amplitude pairs for comparison
        wfA = np.concatenate([ta[1] * np.ones(int(ta[0]))
                              for ta in seqA]) if len(seqA) else np.empty(0)
        wfB = np.concatenate([ta[1] * np.ones(int(ta[0]))
                              for ta in seqB]) if len(seqB) else np.empty(0)
        self.assertTrue(
            len(wfA) == len(wfB), "{} size {} != size {}".format(
                errorHeader, len(wfA), len(wfB)))

        #Check values
        wf_check = np.allclose(wfA, wfB, rtol=0, atol=self.tolerance)

        if not wf_check:
            bad_idx = np.where(wfA != wfB)[0]
            percent_bad = float(len(bad_idx)) / len(wfA)
            bad_level = np.mean(np.abs(wfA - wfB)[bad_idx]) / self.tolerance
            if percent_bad < 0.6:
                msg = "{0}.\nFailed indices: ({1:.1f}% mismatch)\n{2}".format(
                    errorHeader, 100 * percent_bad, bad_idx)
                msg += "\nAvg failure level: {0}".format(bad_level)
            else:
                msg = "{0} ({1:.1f}% mismatch)".format(errorHeader, 100 *
                                                       percent_bad)
        else:
            msg = ""

        self.assertTrue(wf_check, msg=msg)


class TestSequences(object):
    def compare_sequences(self):
        abstract

    #TODO: fix the nutFreq scaling of the arb_axis_drag so that it doesn't overflow
    def test_misc_seqs1(self):
        """ catch all for sequences not otherwise covered """
        self.set_awg_dir()
        seqs = [Ytheta(self.q1, amp=0.5), Z90m(self.q1),
                Ztheta(self.q1, angle=np.pi / 4),
                arb_axis_drag(self.q1, 10.0e6, np.pi / 4, np.pi / 2, np.pi / 8)]

        for ac in range(0, 24):
            seqs.append(AC(self.q1, ac))

        seqs = [seqs]

        filenames = compile_to_hardware(seqs, 'MISC1/MISC1')
        self.compare_sequences('MISC1')

    def test_misc_seqs2(self):
        """ catch all for sequences not otherwise covered """
        self.set_awg_dir()
        seqs = [[ZX90_CR(self.q1, self.q2)]]

        filenames = compile_to_hardware(seqs, 'MISC2/MISC2')
        self.compare_sequences('MISC2')

    def test_misc_seqs3(self):
        """ catch all for sequences not otherwise covered """
        self.set_awg_dir()
        seqs = [[CNOT_CR(self.q1, self.q2)]]

        filenames = compile_to_hardware(seqs, 'MISC3/MISC3')
        self.compare_sequences('MISC3')

    def test_misc_seqs4(self):
        """ catch all for sequences not otherwise covered """
        self.set_awg_dir()
        seqs = [[CNOT_CR(self.q2, self.q1)]]

        filenames = compile_to_hardware(seqs, 'MISC4/MISC4')
        self.compare_sequences('MISC4')

    # TODO: replace with a [MEAS(q1)*MEAS(q2)] sequence where M-q1 and M-q2 share
    # a physical channel.
    # def test_misc_seqs5(self):
    # 	""" catch all for sequences not otherwise covered """
    # 	self.set_awg_dir()
    # 	seqs = [[MEASmux((self.q1, self.q2))]]

    # 	filenames = compile_to_hardware(seqs, 'MISC5/MISC5')
    # 	self.compare_sequences('MISC5')

    def test_AllXY(self):
        self.set_awg_dir()
        AllXY(self.q1)
        self.compare_sequences('AllXY')

    def test_CR_PiRabi(self):
        self.set_awg_dir()
        PiRabi(self.q1, self.q2, np.linspace(0, 4e-6, 11))
        self.compare_sequences('PiRabi')

    def test_CR_EchoCRLen(self):
        self.set_awg_dir('EchoCRLen')
        EchoCRLen(self.q1, self.q2, np.linspace(0, 2e-6, 11))
        self.compare_sequences('EchoCR')

    def test_CR_EchoCRPhase(self):
        self.set_awg_dir('EchoCRPhase')
        EchoCRPhase(self.q1, self.q2, np.linspace(0, pi / 2, 11))
        self.compare_sequences('EchoCR')

    def test_Decoupling_HannEcho(self):
        self.set_awg_dir()
        HahnEcho(self.q1, np.linspace(0, 5e-6, 11))
        self.compare_sequences('Echo')

    def test_Decoupling_CPMG(self):
        self.set_awg_dir()
        CPMG(self.q1, 2 * np.arange(4), 500e-9)
        self.compare_sequences('CPMG')

    def test_FlipFlop(self):
        self.set_awg_dir()
        FlipFlop(self.q1, np.linspace(0, 5e-6, 11))
        self.compare_sequences('FlipFlop')

    def test_T1T2_InversionRecovery(self):
        self.set_awg_dir()
        InversionRecovery(self.q1, np.linspace(0, 5e-6, 11))
        self.compare_sequences('T1')

    def test_T1T2_Ramsey(self):
        self.set_awg_dir()
        Ramsey(self.q1, np.linspace(0, 5e-6, 11))
        self.compare_sequences('Ramsey')

    def test_SPAM(self):
        self.set_awg_dir()
        SPAM(self.q1, np.linspace(0, pi / 2, 11))
        self.compare_sequences('SPAM')

    def test_Rabi_RabiAmp(self):
        self.set_awg_dir('RabiAmp')
        RabiAmp(self.q1, np.linspace(0, 5e-6, 11))
        self.compare_sequences('Rabi')

    def test_Rabi_RabiWidth(self):
        self.set_awg_dir('RabiWidth')
        RabiWidth(self.q1, np.linspace(0, 5e-6, 11))
        self.compare_sequences('Rabi')

    def test_Rabi_RabiAmp_NQubits(self):
        self.set_awg_dir('RabiAmp2')
        RabiAmp_NQubits((self.q1, self.q2), np.linspace(0, 5e-6, 11))
        self.compare_sequences('Rabi')

    def test_Rabi_RabiAmpPi(self):
        self.set_awg_dir('RabiAmpPi')
        RabiAmpPi(self.q1, self.q2, np.linspace(0, 5e-6, 11))
        self.compare_sequences('Rabi')

    def test_Rabi_SingleShot(self):
        self.set_awg_dir()
        SingleShot(self.q1)
        self.compare_sequences('SingleShot')

    def test_Rabi_PulsedSpec(self):
        self.set_awg_dir()
        PulsedSpec(self.q1)
        self.compare_sequences('Spec')

    def test_RB_SingleQubitRB(self):
        self.set_awg_dir('SingleQubitRB')
        np.random.seed(20152606)  # set seed for create_RB_seqs()
        random.seed(20152606)  # set seed for random.choice()
        SingleQubitRB(self.q1, create_RB_seqs(1, 2**np.arange(1, 7)))
        self.compare_sequences('RB')

    @unittest.expectedFailure
    def test_RB_TwoQubitRB(self):
        """  Fails on APS1, APS2, and Tek7000 due to:
		File "QGL\PatternUtils.py", line 129, in add_gate_pulses
    	if has_gate(chan) and not pulse.isZero and not (chan.gate_chan
		AttributeError: 'CompositePulse' object has no attribute 'isZero'
		"""
        self.set_awg_dir('TwoQubitRB')
        np.random.seed(20152606)  # set seed for create_RB_seqs()
        TwoQubitRB(self.q2,
                   self.q1,
                   create_RB_seqs(2, [2, 4, 8, 16, 32],
                                  repeats=16))
        self.compare_sequences('RB')

    # def test_RB_SingleQubitRB_AC(self):
    # 	self.set_awg_dir('SingleQubitRB_AC')
    # 	np.random.seed(20152606) # set seed for create_RB_seqs
    # 	SingleQubitRB_AC(self.q1,create_RB_seqs(1, 2**np.arange(1,7)))
    # 	self.compare_sequences('RB')

    # def test_RB_SingleQubitIRB_AC(self):
    # 	self.set_awg_dir('SingleQubitIRB_AC')
    # 	SingleQubitIRB_AC(self.q1,'')
    # 	self.compare_sequences('RB')

    # def test_RB_SingleQubitRBT(self):
    # 	self.set_awg_dir('SingleQubitRBT')
    # 	SingleQubitRBT(self.q1,'')
    # 	self.compare_sequences('RBT')

    def test_RB_SimultaneousRB_AC(self):
        self.set_awg_dir('SimultaneousRB_AC')
        np.random.seed(20151709)  # set seed for create_RB_seqs
        seqs1 = create_RB_seqs(1, 2**np.arange(1, 7))
        seqs2 = create_RB_seqs(1, 2**np.arange(1, 7))
        SimultaneousRB_AC((self.q1, self.q2), (seqs1, seqs2))
        self.compare_sequences('RB')

    @unittest.skip("Need to update to new pygsti API")
    def test_1Q_GST(self):

        if istravis:
            raise unittest.SkipTest("FIX ME -- Figure out pygsti integration for Travis.")

        self.set_awg_dir('GST')
        # list of GST gate strings
        if GSTTools.PYGSTI_PRESENT:
            # generate the gate gatestrings
            import pygsti
            from pygsti.construction import std1Q_XYI

            #Create a data set
            gs_target = std1Q_XYI.gs_target
            fiducials = std1Q_XYI.fiducials
            germs = std1Q_XYI.germs
            maxLengths = [1,2,4]

            listOfExperiments = pygsti.construction.make_lsgst_experiment_list(gs_target.gates.keys(), fiducials, fiducials, germs, maxLengths)
        else:
            listOfExperiments = list(np.load('tests/test_data/awg/TestAPS2/GST/GST/listOfExperiments.npy'))

        seqs = list(GSTTools.gst_map_1Q(listOfExperiments, self.q1))
        filenames = compile_to_hardware(seqs, 'GST/GST')
        self.compare_sequences('GST')

    @unittest.skip("Need to update to new pygsti API")
    def test_2Q_GST(self):

        if istravis:
            raise unittest.SkipTest("FIX ME -- Figure out pygsti integration for Travis.")

        self.set_awg_dir('GST')
        def gst_2Qgate_map(q1, q2):
            return {"Gxi": X90(q1)*Id(q2),
                     "Gyi": Y90(q1)*Id(q2),
                     "Gii": Id(q1)*Id(q2),
                     "Gix": Id(q1)*X90(q2),
                     "Giy": Id(q1)*Y90(q2),
                     "Gcnot": CNOT_CR(q2,q1)}

        if GSTTools.PYGSTI_PRESENT:
            import pygsti
            from pygsti.construction import std1Q_XYI, std2Q_XYICNOT
            from itertools import product
            from QGL.GSTTools import SingleQubitCliffordGST, gst_map_2Q
            # note the use of the germs_lite!
            gs = std2Q_XYICNOT
            gs_target = std2Q_XYICNOT.gs_target.copy()

            prep_fiducials = std2Q_XYICNOT.prepStrs
            effect_fiducials = std2Q_XYICNOT.effectStrs
            gs_germs = std2Q_XYICNOT.germs_lite
            #maxLengths = [1,2,4,8,16]
            maxLengths = [1,2]

            print('Creating GST sequences...')
            listOfExperiments = pygsti.construction.make_lsgst_experiment_list(gs_target.gates.keys(), prep_fiducials, effect_fiducials, gs_germs, maxLengths)
        else:
            # list of GST gate strings
            listOfExperiments = list(np.load('tests/test_data/awg/TestAPS2/GST/GST2Q/listOfExperiments.npy'))

        seqs = list(GSTTools.gst_map_2Q(listOfExperiments, (self.q1, self.q2), qgl_map = gst_2Qgate_map(self.q1, self.q2), append_meas=True))

        filenames = compile_to_hardware(seqs, 'GST/GST2Q')
        self.compare_sequences('GST2Q')


class APS2Helper(AWGTestHelper):
    def setUp(self):
        AWGTestHelper.__init__(self, APS2Pattern)
        for name in ['APS1', 'APS2', 'APS3', 'APS4', 'APS5', 'APS6']:
            channelName = name + '-1'
            channel = PhysicalQuadratureChannel(label=channelName, channel=0)
            channel.sampling_rate = 1.2e9
            channel.instrument = name
            channel.translator = 'APS2Pattern'
            self.channels[channelName] = channel

            for m in range(1, 5):
                channelName = "{0}-m{1}".format(name, m)
                channel = PhysicalMarkerChannel(label=channelName, channel=m-1)
                channel.sampling_rate = 1.2e9
                channel.instrument = name
                channel.translator = 'APS2Pattern'
                self.channels[channelName] = channel

        mapping = {'digitizerTrig': 'APS1-m1',
                   'slave_trig': 'APS1-m2',
                   'q1': 'APS1-1',
                   'q1-gate': 'APS1-m3',
                   'M-q1': 'APS2-1',
                   'M-q1-gate': 'APS2-m1',
                   'q2': 'APS3-1',
                   'q2-gate': 'APS3-m1',
                   'M-q2': 'APS4-1',
                   'M-q2-gate': 'APS4-m1',
                   'cr': 'APS5-1',
                   'cr-gate': 'APS5-m1',
                   'M-q1q2': 'APS6-1',
                   'M-q1q2-gate': 'APS6-m1'}

        self.finalize_map(mapping)


class TestAPS2(unittest.TestCase, APS2Helper, TestSequences):
    # TestAPS2 is seperated from APS2Helper so the setup method of APS2Helper may be used in test_json.py
    def setUp(self):
        APS2Helper.__init__(self)
        APS2Helper.setUp(self)

    def test_mux_CR(self):
        self.set_awg_dir()
        #control and CR sharing the same chans
        self.channels['cr'].phys_chan = self.channels['q1'].phys_chan
        self.channels['q1'].frequency = 100e6
        self.channels['cr'].frequency = 200e6
        self.cl.update_channelDict()
        seqs = [[CNOT_CR(self.q1, self.q2)]]

        filenames = compile_to_hardware(seqs, 'CNOT_CR_mux/CNOT_CR_mux')
        self.compare_sequences('CNOT_CR_mux')

    def test_update_in_place(self):
        self.set_awg_dir('RabiAmpInPlace')
        APS2Pattern.SAVE_WF_OFFSETS = True
        RabiAmp(self.q1, np.linspace(0, 5e-6, 11))
        with open(os.path.join(self.awg_dir, "Rabi", "Rabi-APS1.offsets"), "rb") as FID:
            offsets = pickle.load(FID)
        print(offsets)
        pulses = {list(offsets.keys())[0]: Utheta(self.q1, amp=0.0, phase=0)}
        APS2Pattern.update_wf_library(os.path.join(self.awg_dir, "Rabi", "Rabi-APS1.aps2"), pulses, offsets)
        print(os.path.join(self.awg_dir, "Rabi", "Rabi-APS1.aps2"))
        APS2Pattern.SAVE_WF_OFFSETS = False


class TestAPS1(unittest.TestCase, AWGTestHelper, TestSequences):
    def setUp(self):
        AWGTestHelper.__init__(self, APSPattern)
        for name in ['APS1', 'APS2', 'APS3']:
            for i, ch in enumerate(['12', '34']):
                channelName = name + '-' + ch
                channel = PhysicalQuadratureChannel(label=channelName, channel=i)
                channel.sampling_rate = 1.2e9
                channel.instrument = name
                channel.translator = 'APSPattern'
                self.channels[channelName] = channel

            for m in range(1, 5):
                channelName = "{0}-{1}m1".format(name, m)
                channel = PhysicalMarkerChannel(label=channelName, channel=m-1)
                channel.sampling_rate = 1.2e9
                channel.instrument = name
                channel.translator = 'APSPattern'
                self.channels[channelName] = channel

        mapping = {'digitizerTrig': 'APS1-1m1',
                   'slave_trig': 'APS1-2m1',
                   'q1': 'APS1-12',
                   'M-q1': 'APS1-34',
                   'M-q1-gate': 'APS1-3m1',
                   'q1-gate': 'APS1-4m1',
                   'q2': 'APS2-12',
                   'M-q2': 'APS2-34',
                   'M-q2-gate': 'APS2-1m1',
                   'q2-gate': 'APS2-2m1',
                   'cr': 'APS3-12',
                   'cr-gate': 'APS3-1m1',
                   'M-q1q2': 'APS3-34',
                   'M-q1q2-gate': 'APS3-2m1'}

        # override trigger lengths on APS1 to get single blips
        self.channels['slave_trig'].pulse_params['length'] = 0.833e-9
        self.channels['digitizerTrig'].pulse_params['length'] = 0.833e-9
        self.finalize_map(mapping)

    def compare_file_data(self, testFile, truthFile):
        '''
		Override the method in AWGTestHelper so that we can special-case marker comparison
		'''
        awgData = self.read_function(testFile)
        truthData = self.read_function(truthFile)

        awgDataLen = len(awgData)
        truthDataLen = len(truthData)

        self.assertTrue(awgDataLen == truthDataLen,
                        "Expected {0} sequences in file. Found {1}.".format(
                            truthDataLen, awgDataLen))

        for name in truthData:
            self.assertTrue(
                name in awgData,
                "Expected channel {0} not found in file {1}".format(name,
                                                                    testFile))
            isMarker = ('m' == name[-2])

            for x in range(len(truthData[name])):
                seqA = np.array(truthData[name][x])
                seqB = np.array(awgData[name][x])
                if isMarker:
                    self.compare_marker_sequence(
                        seqA, seqB,
                        "\nFile {0} =>\nChannel {1} Sequence {2}".format(
                            testFile, name, x))
                else:
                    self.compare_sequence(
                        seqA, seqB,
                        "\nFile {0} =>\nChannel {1} Sequence {2}".format(
                            testFile, name, x))

    def compare_marker_sequence(self, seqA, seqB, errorHeader):
        markerDistanceTolerance = 4
        self.assertTrue(seqA.size == seqB.size,
                        "{0} size {1} != size {2}".format(
                            errorHeader, str(seqA.size), str(seqB.size)))

        # convert sequences to locations of blips
        idxA = np.where(seqA)[0]
        idxB = np.where(seqB)[0]
        self.assertTrue(
            len(idxA) == len(idxB),
            "{0}.\nNumber of blips did not match: {1} != {2}".format(
                errorHeader, len(idxA), len(idxB)))
        # compare the blip locations element-wise
        if len(idxA) > 0:
            diff = np.abs(idxA - idxB)
        else:
            diff = np.array([0])
        self.assertTrue(
            max(diff) <= markerDistanceTolerance,
            "{0}\nMismatches: {1}".format(errorHeader, diff.nonzero()[0]))

    @unittest.expectedFailure
    def test_Rabi_RabiWidth(self):
        """ test_Rabi_RabiWidth is expected to fail on APS1 with the following error:
			AssertionError: Oops! You have exceeded the waveform memory of the APS
		"""
        TestSequences.test_Rabi_RabiWidth(self)

    @unittest.expectedFailure
    def test_RB_SimultaneousRB_AC(self):
        """ test_RB_SimultaneousRB_AC is expected to fail on APS1 with the following error:
			AssertionError: Oops! You have exceeded the waveform memory of the APS
		"""
        TestSequences.test_RB_SimultaneousRB_AC(self)

if __name__ == "__main__":
    unittest.main()
