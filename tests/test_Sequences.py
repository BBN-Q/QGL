import h5py
import numpy as np
import unittest, time, os, random

import config
import Libraries
from QGL import *
import QGL

from QGL.Channels import Edge, Measurement, LogicalChannel, LogicalMarkerChannel, PhysicalMarkerChannel, PhysicalQuadratureChannel, ChannelLibrary
from instruments.AWGs import APS2, APS, Tek5014, Tek7000
from instruments.InstrumentManager import InstrumentLibrary

BASE_AWG_DIR = config.AWGDir

class AWGTestHelper(object):
	testFileDirectory = './tests/test_data/awg/'

	def __init__(self, clsObj = None, tolerance = 1.5/2**13):
		self.channels = {}
		self.instruments = {}
		self.assign_channels()
		self.set_awg_dir()
		if clsObj:
			self.read_function = clsObj().read_sequence_file
		else:
			self.read_function = None
		self.tolerance = tolerance

	def finalize_map(self, mapping):
		for name,value in mapping.iteritems():
			self.channels[name].physChan = self.channels[value]

		Compiler.channelLib = ChannelLibrary()
		Compiler.channelLib.channelDict = self.channels
		Compiler.channelLib.build_connectivity_graph()

		Compiler.instrumentLib = InstrumentLibrary()
		Compiler.instrumentLib.instrDict = self.instruments

		(self.q1, self.q2) = self.get_qubits()

	def assign_channels(self):

		self.qubit_names = ['q1','q2']
		self.logical_names = ['digitizerTrig', 'slaveTrig']

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

			m = Measurement(label=mName, gateChan = mg, trigChan=self.channels['digitizerTrig'])

			q = Qubit(label=name, gateChan=qg)
			q.pulseParams['length'] = 30e-9
			q.pulseParams['phase'] = pi/2

			self.channels[name] = q
			self.channels[mName] = m
			self.channels[mgName]  = mg
			self.channels[qgName]  = qg

		# this block depends on the existence of q1 and q2
		self.channels['cr-gate'] = LogicalMarkerChannel(label='cr-gate')

		q1, q2 = self.channels['q1'], self.channels['q2']
		cr = Edge(label="cr", source = q1, target = q2, gateChan = self.channels['cr-gate'] )
		cr.pulseParams['length'] = 30e-9
		cr.pulseParams['phase'] = pi/4
		self.channels["cr"] = cr

		mq1q2g = LogicalMarkerChannel(label='M-q1q2-gate')
		self.channels['M-q1q2-gate']  = mq1q2g
		self.channels['M-q1q2']       = Measurement(label='M-q1q2', gateChan = mq1q2g, trigChan=self.channels['digitizerTrig'])

	def get_qubits(self):
		return [QGL.Compiler.channelLib[name] for name in self.qubit_names]

	def set_awg_dir(self, footer = ""):
		cn = self.__class__.__name__

		self.awg_dir = os.path.abspath(BASE_AWG_DIR + os.path.sep + cn)
		self.truth_dir = os.path.abspath(self.testFileDirectory + os.path.sep + cn)

		if footer != "":
			self.awg_dir = self.awg_dir + os.path.sep + footer
			self.truth_dir = self.truth_dir + os.path.sep + footer

		if not os.path.isdir(self.awg_dir):
			os.makedirs(self.awg_dir)
		config.AWGDir = self.awg_dir

	def compare_sequences(self, seqDir):
		if not self.read_function:
			print "AWGTestHelper.read_function is not defined"
			return

		searchDirectory = self.awg_dir + os.path.sep + seqDir
		truthDirectory = self.truth_dir + os.path.sep + seqDir

		filenames = os.listdir(searchDirectory)
		for filename in filenames:
			truthFile = truthDirectory + os.path.sep + filename
			testFile = searchDirectory + os.path.sep + filename

			self.assertTrue(os.path.isfile(truthFile), "Truth Data File: {0} not found.".format(truthFile) )
			self.compare_file_data(testFile, truthFile)

	def compare_file_data(self, testFile, truthFile):
		awgData = self.read_function(testFile)
		truthData = self.read_function(truthFile)

		awgDataLen = len(awgData)
		truthDataLen = len(truthData)

		self.assertTrue(awgDataLen == truthDataLen,
			"Expected {0} sequences in file. Found {1}.".format(truthDataLen, awgDataLen))

		for name in truthData:
			self.assertTrue(name in awgData, "Expected channel {0} not found in file {1}".format(name, testFile))

			for x in range(len(truthData[name])):
				seqA = np.array(truthData[name][x])
				seqB = np.array(awgData[name][x])
				self.compare_sequence(seqA,seqB,  "\nFile {0} =>\nChannel {1} Sequence {2}".format(testFile, name, x))

	def compare_sequence(self, seqA, seqB, errorHeader):
		self.assertTrue( seqA.size == seqB.size, "{0} size {1} != size {2}".format(errorHeader, str(seqA.size), str(seqB.size)))
		# np.testing.assert_allclose(seqA, seqB, rtol=0, atol=self.tolerance, err_msg=errorHeader)
		npdiff = np.allclose(seqA, seqB, rtol=0, atol=self.tolerance)
		diff = np.abs(seqA - seqB) < self.tolerance
		test = npdiff or all(diff)
		if not test:
			bad_idx = np.where(diff == False)[0]
			percent_bad = float(len(bad_idx))/len(seqA)
			bad_level = np.mean(np.abs(seqA - seqB)[bad_idx]) / self.tolerance
			if percent_bad < 0.6:
				msg = "{0}.\nFailed indices: ({1:.1f}% mismatch)\n{2}".format(errorHeader, 100*percent_bad, bad_idx)
				msg += "\nAvg failure level: {0}".format(bad_level)
			else:
				msg = "{0} ({1:.1f}% mismatch)".format(errorHeader, 100*percent_bad)
		else:
			msg = ""
		self.assertTrue(npdiff or all(diff), msg)

class TestSequences(object):

	def compare_sequences(self): abstract

	def test_misc_seqs1(self):
		""" catch all for sequences not otherwise covered """
		self.set_awg_dir()
		seqs = [Ytheta(self.q1, amp = 0.5), Z90m(self.q1), Ztheta(self.q1, angle = np.pi/4),
				arb_axis_drag(self.q1, 5.0, np.pi/4, np.pi/2, np.pi/8)
			   ]

		for ac in range(0, 24):
			seqs.append(AC(self.q1, ac))

		seqs = [seqs]

		filenames = compile_to_hardware(seqs, 'MISC1/MISC1')
		self.compare_sequences('MISC1')

	def test_misc_seqs2(self):
		""" catch all for sequences not otherwise covered """
		self.set_awg_dir()
		seqs = [ZX90_CR(self.q1, self.q2)]

		filenames = compile_to_hardware(seqs, 'MISC2/MISC2')
		self.compare_sequences('MISC2')

	def test_misc_seqs3(self):
		""" catch all for sequences not otherwise covered """
		self.set_awg_dir()
		seqs = [CNOT_CR(self.q1, self.q2)]

		filenames = compile_to_hardware(seqs, 'MISC3/MISC3')
		self.compare_sequences('MISC3')

	def test_misc_seqs4(self):
		""" catch all for sequences not otherwise covered """
		self.set_awg_dir()
		seqs = [CNOT_CR(self.q2, self.q1)]

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
		EchoCRPhase(self.q1, self.q2, np.linspace(0, pi/2, 11))
		self.compare_sequences('EchoCR')

	def test_Decoupling_HannEcho(self):
		self.set_awg_dir()
	 	HahnEcho(self.q1, np.linspace(0, 5e-6, 11))
	 	self.compare_sequences('Echo')

	def test_Decoupling_CPMG(self):
		self.set_awg_dir()
		CPMG(self.q1, 2*np.arange(4), 500e-9)
	 	self.compare_sequences('CPMG')

	def test_FlipFlop(self):
		self.set_awg_dir()
	  	FlipFlop(self.q1, np.linspace(0, 5e-6, 11))
	  	self.compare_sequences('FlipFlop')

	def test_T1T2_InversionRecovery(self):
		self.set_awg_dir()
	  	InversionRecovery(self.q1,  np.linspace(0, 5e-6, 11))
	  	self.compare_sequences('T1')

	def test_T1T2_Ramsey(self):
		self.set_awg_dir()
	 	Ramsey(self.q1, np.linspace(0, 5e-6, 11))
	 	self.compare_sequences('Ramsey')

	def test_SPAM(self):
		self.set_awg_dir()
	  	SPAM(self.q1, np.linspace(0, pi/2, 11))
	  	self.compare_sequences('SPAM')

	def test_Rabi_RabiAmp(self):
		self.set_awg_dir('RabiAmp')
		RabiAmp(self.q1,  np.linspace(0, 5e-6, 11))
		self.compare_sequences('Rabi')

	def test_Rabi_RabiWidth(self):
		self.set_awg_dir('RabiWidth')
		RabiWidth(self.q1,  np.linspace(0, 5e-6, 11))
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
		np.random.seed(20152606) # set seed for create_RB_seqs()
		random.seed(20152606) # set seed for random.choice()
		SingleQubitRB(self.q1, create_RB_seqs(1, 2**np.arange(1,7)))
		self.compare_sequences('RB')

	@unittest.expectedFailure
	def test_RB_TwoQubitRB(self):
		"""  Fails on APS1, APS2, and Tek7000 due to:
		File "QGL\PatternUtils.py", line 129, in add_gate_pulses
    	if has_gate(chan) and not pulse.isZero and not (chan.gateChan
		AttributeError: 'CompositePulse' object has no attribute 'isZero'
		"""
		self.set_awg_dir('TwoQubitRB')
		np.random.seed(20152606) # set seed for create_RB_seqs()
		TwoQubitRB(self.q1, self.q2, create_RB_seqs(2, [2, 4, 8, 16, 32], repeats=16))
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

	@unittest.expectedFailure
	def test_RB_SimultaneousRB_AC(self):
		self.set_awg_dir('SimultaneousRB_AC')
	 	np.random.seed(20151709) # set seed for create_RB_seqs
		seqs1 = create_RB_seqs(1, 2**np.arange(1,7))
		seqs2 = create_RB_seqs(1, 2**np.arange(1,7))
	 	SimultaneousRB_AC((self.q1, self.q2), (seqs1, seqs2))
	 	self.compare_sequences('RB')

class APS2Helper(AWGTestHelper):
	def setUp(self):
		AWGTestHelper.__init__(self, APS2)
		for name in ['APS1', 'APS2', 'APS3', 'APS4', 'APS5', 'APS6']:
			self.instruments[name] = APS2(label=name)

			channelName = name + '-12'
			channel = PhysicalQuadratureChannel(label=channelName)
			channel.AWG = self.instruments[name]
			self.channels[channelName] = channel

			for m in range(1,5):
				channelName = "{0}-12m{1}".format(name,m)
				channel = PhysicalMarkerChannel(label=channelName)
				channel.AWG = self.instruments[name]
				self.channels[channelName] = channel

		mapping = {	'digitizerTrig' : 'APS1-12m1',
				    'slaveTrig'     : 'APS1-12m2',
			        'q1'            : 'APS1-12',
					'q1-gate'       : 'APS1-12m3',
					'M-q1'          : 'APS2-12',
					'M-q1-gate'     : 'APS2-12m1',
					'q2'            : 'APS3-12',
					'q2-gate'       : 'APS3-12m1',
					'M-q2'          : 'APS4-12',
					'M-q2-gate'     : 'APS4-12m1',
					'cr'            : 'APS5-12',
					'cr-gate'       : 'APS5-12m1',
					'M-q1q2'        : 'APS6-12',
					'M-q1q2-gate'   : 'APS6-12m1'}

		self.finalize_map(mapping)

class TestAPS2(unittest.TestCase, APS2Helper, TestSequences):
	# TestAPS2 is seperated from APS2Helper so the setup method of APS2Helper may be used in test_json.py
	def setUp(self):
		APS2Helper.__init__(self)
		APS2Helper.setUp(self)

class TestAPS1(unittest.TestCase, AWGTestHelper, TestSequences):

	def setUp(self):
		AWGTestHelper.__init__(self, APS)
		for name in ['APS1', 'APS2', 'APS3']:
			self.instruments[name] = APS(label=name)

			for ch in ['12', '34']:
				channelName = name + '-' + ch
				channel = PhysicalQuadratureChannel(label=channelName)
				channel.AWG = self.instruments[name]
				self.channels[channelName] = channel

			for m in range(1,5):
				channelName = "{0}-{1}m1".format(name,m)
				channel = PhysicalMarkerChannel(label=channelName)
				channel.AWG = self.instruments[name]
				self.channels[channelName] = channel

		mapping = {	'digitizerTrig' : 'APS1-1m1',
					'slaveTrig'     : 'APS1-2m1',
					'q1'            : 'APS1-12',
					'M-q1'          : 'APS1-34',
					'M-q1-gate'     : 'APS1-3m1',
					'q1-gate'       : 'APS1-4m1',
					'q2'            : 'APS2-12',
					'M-q2'          : 'APS2-34',
					'M-q2-gate'     : 'APS2-1m1',
					'q2-gate'       : 'APS2-2m1',
					'cr'            : 'APS3-12',
					'cr-gate'       : 'APS3-1m1',
					'M-q1q2'        : 'APS3-34',
					'M-q1q2-gate'   : 'APS3-2m1'}

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
			"Expected {0} sequences in file. Found {1}.".format(truthDataLen, awgDataLen))

		for name in truthData:
			self.assertTrue(name in awgData, "Expected channel {0} not found in file {1}".format(name, testFile))
			isMarker = ('m' == name[-2])

			for x in range(len(truthData[name])):
				seqA = np.array(truthData[name][x])
				seqB = np.array(awgData[name][x])
				if isMarker:
					self.compare_marker_sequence(seqA, seqB, "\nFile {0} =>\nChannel {1} Sequence {2}".format(testFile, name, x))
				else:
					self.compare_sequence(seqA, seqB, "\nFile {0} =>\nChannel {1} Sequence {2}".format(testFile, name, x))

	def compare_marker_sequence(self, seqA, seqB, errorHeader):
		markerDistanceTolerance = 4
		self.assertTrue( seqA.size == seqB.size, "{0} size {1} != size {2}".format(errorHeader, str(seqA.size), str(seqB.size)))

		# convert sequences to locations of blips
		idxA = np.where(seqA)[0]
		idxB = np.where(seqB)[0]
		self.assertTrue(len(idxA) == len(idxB), "{0}.\nNumber of blips did not match: {1} != {2}".format(errorHeader, len(idxA), len(idxB)))
		# compare the blip locations element-wise
		if len(idxA) > 0:
			diff = np.abs(idxA - idxB)
		else:
			diff = np.array([0])
		self.assertTrue(max(diff) <= markerDistanceTolerance, "{0}\nMismatches: {1}".format(errorHeader, diff.nonzero()[0]))

	@unittest.expectedFailure
	def test_Rabi_RabiWidth(self):
		""" test_Rabi_RabiWidth is expected to fail on APS1 with the following error:
			AssertionError: Oops! You have exceeded the waveform memory of the APS
		"""
		TestSequences.test_Rabi_RabiWidth(self)

class TestTek5014(unittest.TestCase, AWGTestHelper, TestSequences):

	def setUp(self):
		AWGTestHelper.__init__(self, Tek5014)
		for name in ['TEK1', 'TEK2']:
			self.instruments[name] = Tek5014(label=name)

			for ch in ['12', '34']:
				channelName = name + '-' + ch
				channel = PhysicalQuadratureChannel(label=channelName)
				channel.AWG = self.instruments[name]
				self.channels[channelName] = channel

			for m in ['1m1', '1m2', '2m1', '2m2', '3m1', '3m2', '4m1', '4m2']:
				channelName = "{0}-{1}".format(name,m)
				channel = PhysicalMarkerChannel(label=channelName)
				channel.AWG = self.instruments[name]
				self.channels[channelName] = channel

		mapping = { 'digitizerTrig'	: 'TEK1-1m2',
					'slaveTrig'   	: 'TEK1-2m2',
					'q1'			: 'TEK1-12',
					'M-q1'			: 'TEK1-12',
					'M-q1-gate'		: 'TEK1-1m1',
					'q1-gate'		: 'TEK1-2m1',
					'q2'			: 'TEK1-34',
					'M-q2'			: 'TEK1-34',
					'M-q2-gate'		: 'TEK1-3m1',
					'q2-gate'		: 'TEK1-4m1',
					'cr'            : 'TEK2-12',
					'cr-gate'       : 'TEK2-1m1',
					'M-q1q2'        : 'TEK2-34',
					'M-q1q2-gate'   : 'TEK2-2m1'}

		self.finalize_map(mapping)

	@unittest.expectedFailure
	def test_misc_seqs2(self):
		""" Fails due to a divide by zero
		    File "C:\Projects\Q\lib\PyQLab\QGL\TekPattern.py", line 77, in merge_waveform
   			  for entry in chAB['linkList'][n % len(chAB['linkList'])]:
			ZeroDivisionError: integer division or modulo by zero
		"""
		TestSequences.test_misc_seqs2(self)

	@unittest.expectedFailure
	def test_misc_seqs5(self):
		""" Fails due to a divide by zero
		    File "C:\Projects\Q\lib\PyQLab\QGL\TekPattern.py", line 77, in merge_waveform
   			  for entry in chAB['linkList'][n % len(chAB['linkList'])]:
			ZeroDivisionError: integer division or modulo by zero
		"""
		TestSequences.test_misc_seqs5(self)

	# multiple tests will fail with an attribute error:
	# AttributeError: 'Wait' object has no attribute 'isTimeAmp' at line 78
	# in TekPattern.py in merge_waveform

	@unittest.expectedFailure
	def test_misc_seqs1(self):
		TestSequences.test_misc_seqs1(self)

	@unittest.expectedFailure
	def test_misc_seqs3(self):
		TestSequences.test_misc_seqs3(self)

	@unittest.expectedFailure
	def test_misc_seqs4(self):
		TestSequences.test_misc_seqs4(self)

	@unittest.expectedFailure
	def test_AllXY(self):
		TestSequences.test_AllXY(self)

	@unittest.expectedFailure
	def test_CR_PiRabi(self):
		TestSequences.test_CR_PiRabi(self)

	@unittest.expectedFailure
	def test_CR_EchoCRLen(self):
		TestSequences.test_CR_EchoCRLen(self)

	@unittest.expectedFailure
	def test_CR_EchoCRPhase(self):
		TestSequences.test_CR_EchoCRPhase(self)

	@unittest.expectedFailure
	def test_Decoupling_HannEcho(self):
		TestSequences.test_Decoupling_HannEcho(self)

	@unittest.expectedFailure
	def test_Decoupling_CPMG(self):
		TestSequences.test_Decoupling_CPMG(self)

	@unittest.expectedFailure
	def test_FlipFlop(self):
		TestSequences.test_FlipFlop(self)

	@unittest.expectedFailure
	def test_T1T2_InversionRecovery(self):
		TestSequences.test_T1T2_InversionRecovery(self)

	@unittest.expectedFailure
	def test_T1T2_Ramsey(self):
		TestSequences.test_T1T2_Ramsey(self)

	@unittest.expectedFailure
	def test_SPAM(self):
		TestSequences.test_SPAM(self)

	@unittest.expectedFailure
	def test_Rabi_RabiAmp(self):
		TestSequences.test_Rabi_RabiAmp(self)

	@unittest.expectedFailure
	def test_Rabi_RabiWidth(self):
		TestSequences.test_Rabi_RabiWidth(self)

	@unittest.expectedFailure
	def test_Rabi_RabiAmp_NQubits(self):
		TestSequences.test_Rabi_RabiAmp_NQubits(self)

	@unittest.expectedFailure
	def test_Rabi_RabiAmpPi(self):
		TestSequences.test_Rabi_RabiAmpPi(self)

	@unittest.expectedFailure
	def test_Rabi_SingleShot(self):
		TestSequences.test_Rabi_SingleShot(self)

	@unittest.expectedFailure
	def test_Rabi_PulsedSpec(self):
		TestSequences.test_Rabi_PulsedSpec(self)

	@unittest.expectedFailure
	def test_RB_SingleQubitRB(self):
		TestSequences.test_RB_SingleQubitRB(self)

# class TestTek7000(unittest.TestCase, AWGTestHelper, TestSequences):

# 	def setUp(self):
# 		AWGTestHelper.__init__(self, TekPattern.read_Tek_file)
# 		for name in ['TEK1', 'TEK2', 'TEK3', 'TEK4', 'TEK5']:
# 			self.instruments[name] = Tek7000(label=name)

# 			for ch in ['12']:
# 				channelName = name + '-' + ch
# 				channel = PhysicalQuadratureChannel(label=channelName)
# 				channel.AWG = self.instruments[name]
# 				self.channels[channelName] = channel

# 			for m in ['1m1', '1m2', '2m1', '2m2']:
# 				channelName = "{0}-{1}".format(name,m)
# 				channel = PhysicalMarkerChannel(label=channelName)
# 				channel.AWG = self.instruments[name]
# 				self.channels[channelName] = channel

# 		mapping = { 'digitizerTrig'	:'TEK1-1m2',
# 					'slaveTrig'   	:'TEK1-2m2',
# 					'q1'			:'TEK1-12',
# 					'M-q1'			:'TEK2-12',
# 					'M-q1-gate'		:'TEK1-1m1',
# 					'q1-gate'		:'TEK1-2m1',
# 					'q2'			:'TEK3-12',
# 					'M-q2'			:'TEK4-12',
# 					'M-q2-gate'		:'TEK2-1m1',
# 					'q2-gate'		:'TEK2-2m1',
# 					'cr'            :'TEK5-12',
# 					'cr-gate'       :'TEK5-1m1'}
# 		self.finalize_map(mapping)

if __name__ == "__main__":
    unittest.main()
