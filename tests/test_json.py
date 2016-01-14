import unittest
import os
import json
import config
import numpy as np
from atom.api import Atom, atomlist

import Libraries
from MeasFilters import RawStream, DigitalDemod, KernelIntegration, Correlator, StateComparator, StreamSelector
from Sweeps import *
from test_Sequences import APS2Helper
import ExpSettingsVal


testsRunning = 0

class JSONTestHelper(object):

	
	testFileNames = [Libraries.channelLib.libFile, 
					 Libraries.instrumentLib.libFile,
					 Libraries.measLib.libFile,
					 Libraries.sweepLib.libFile]

	@classmethod
	def backupFiles(cls):		
		global testsRunning
		testsRunning = testsRunning + 1
		if testsRunning != 1:
			return
		for fileName in cls.testFileNames:
			print "Testing for file {0}".format(fileName)
			if os.path.exists(fileName):
				print "Backing up file {0}".format(fileName)
				if os.path.exists(fileName + '.org'):
					os.remove(fileName + '.org')
		 		os.rename(fileName, fileName + '.org')

	@classmethod
	def restoreFiles(cls):
		global testsRunning
		testsRunning = testsRunning - 1
		if testsRunning > 0:
			return
		for fileName in cls.testFileNames:
			print "Testing for file {0}".format(fileName)
			if os.path.exists(fileName):
				os.remove(fileName)
			if os.path.exists(fileName+ '.org'):
				print "Restoring file {0}".format(fileName)
	 			os.rename(fileName + '.org', fileName)

	def validate_json_dictionary(self, testDict, validDict):

		# top level
		for key in validDict.keys():
			print key
			assertIn(key, testDict)
			if isinstance(testDict[key], dict):
				print testDict[key]

	def validate_atomlist(self, test, truth):
		for i in range(len(truth)):
			if isinstance(truth[i], Atom):
				self.validate_atom(test[i], truth[i])
			elif isinstance(truth[i], atomlist):
				self.validate_atomlist(test[i], truth[i])
			else:
				self.assertEqual(test[i], truth[i])

	def validate_atom(self, test, truth):
		self.assertIsInstance(test, truth.__class__)
		testState = test.__getstate__()
		truthState = truth.__getstate__()
		for state in truthState:
			if isinstance(truthState[state], Atom):
				self.validate_atom(testState[state], truthState[state])
			elif isinstance(truthState[state], atomlist):
				self.validate_atomlist(testState[state], truthState[state])
			elif isinstance(truthState[state], np.ndarray):
				np.testing.assert_allclose(testState[state], truthState[state], rtol=1e-5, atol=0)
			else:
				self.assertEqual(testState[state], truthState[state])

	def validate_library(self, test, truth):
		for element in truth:
			self.assertIn(element, test)
			self.assertIsInstance(test[element], truth[element].__class__)
			if isinstance(truth[element], Atom):
				self.validate_atom(test[element], truth[element])
			else:
				print "Did not test: ", element


class TestAWGJSON(unittest.TestCase, APS2Helper, JSONTestHelper):

	
	@classmethod
	def setUpClass(cls):
		cls.backupFiles()

	@classmethod
	def tearDownClass(cls):
		cls.restoreFiles()

	def setUp(self):
		APS2Helper.__init__(self)
		APS2Helper.setUp(self)

	def test_channels_instruments_library(self):
		# test channels and instruments together as they are coupled

		Libraries.channelLib.channelDict = self.channels
		Libraries.instrumentLib.instrDict = self.instruments

		ExpSettingsVal.list_config()
		ExpSettingsVal.validate_lib()

		Libraries.channelLib.write_to_file()
		Libraries.instrumentLib.write_to_file()

		Libraries.channelLib.channelDict = {}
		Libraries.channelLib.load_from_library()

		Libraries.instrumentLib.instrDict = {}
		Libraries.instrumentLib.load_from_library()

		self.validate_library(Libraries.channelLib.channelDict, self.channels)
		self.validate_library(Libraries.instrumentLib.instrDict, self.instruments)

class TestMeasJSON(unittest.TestCase, JSONTestHelper):

	

	@classmethod
	def setUpClass(cls):
		cls.backupFiles()

	@classmethod
	def tearDownClass(cls):
		cls.restoreFiles()

	def setUp(self):
		
		self.measurements = {}
		self.measurements['R1'] = RawStream(label='R1', saveRecords = True,recordsFilePath = '/tmp/records', channel='1')
		self.measurements['M1'] = DigitalDemod(label='M1',  saveRecords = True,recordsFilePath = '/tmp/records', IFfreq=10e6, samplingRate=250e6)
		#self.measurements['M12'] = Correlator(label='M12')
		self.measurements['KI1'] = KernelIntegration(label='KI1', boxCarStart=100, boxCarStop=500, IFfreq=10e6, samplingRate=250e6)
		self.measurements['SC1'] = StateComparator(label='SC1', threshold = 0.5, integrationTime = 100)
		self.measurements['SS1'] = StreamSelector(label='SS1', stream = 'test', saveRecords = True, recordsFilePath = '/tmp/records')
	
	def test_measurements_library(self):
		Libraries.measLib.filterDict = self.measurements

		ExpSettingsVal.list_config()
		ExpSettingsVal.validate_lib()

		Libraries.measLib.write_to_file()
		Libraries.measLib.filterDict = {}
		Libraries.measLib.load_from_library()

		# can not test instrumentLib directly here must test enclosed dictionary
		self.validate_library(Libraries.measLib.filterDict, self.measurements)

class TestSweepJSON(unittest.TestCase, JSONTestHelper):

	

	@classmethod
	def setUpClass(cls):
		cls.backupFiles()

	@classmethod
	def tearDownClass(cls):
		cls.restoreFiles()


	def setUp(self):

		self.sweeps = {}
		self.sweeps['PS'] = PointsSweep(label='PS', start = 1.0 , stop = 10.0, numPoints = 10 )
		self.sweeps['PW'] = Power(label='PW', instr = 'CSSRC',
										units = 'dBm', start = 1.0 ,
										stop = 10.0, numPoints = 10 )

		self.sweeps['F'] = Frequency(label='F', instr = 'CSSRC',
									   start = 1.0 , stop = 10.0, numPoints = 10)

		self.sweeps['HF'] = HeterodyneFrequency(label='HF', instr1 = 'CSSRC1', instr2 = 'CSSRC2',
									   start = 1.0 , stop = 10.0, numPoints = 10)

		self.sweeps['SN'] = SegmentNum(label='SN', start = 1.0 , stop = 10.0, numPoints = 10)

		self.sweeps['SNWC1'] = SegmentNumWithCals(label='SNWC', numCals = 2, start = 1.0 , stop = 10.0, numPoints = 10)

		self.sweeps['R'] = Repeat(label='R', numRepeats = 2)

		self.sweeps['AWGC'] = AWGChannel(channel = '1', mode = 'Amp.',
										start = 1.0 , stop = 10.0, numPoints = 10 )

		self.sweeps['AWGS'] = AWGSequence(sequenceFile = '/tmp/sf',
										start = 1 , stop = 10, step = 2 )

		self.sweeps['A'] = Attenuation(channel= 1 , instr = 'AInst', start = 1.0 , stop = 10.0, numPoints = 10 )
		self.sweeps['D'] = DC( instr = 'DCInst', start = 1.0 , stop = 10.0, numPoints = 10 )
		self.sweeps['T'] = Threshold( instr = 'TInst', stream = '(1,1)', start = 1.0 , stop = 10.0, numPoints = 10 )

	def test_sweeps_library(self):
		Libraries.sweepLib.sweepDict = self.sweeps

		ExpSettingsVal.list_config()
		ExpSettingsVal.validate_lib()

		Libraries.sweepLib.write_to_file()
		Libraries.sweepLib.sweepDict = {}
		Libraries.sweepLib.load_from_library()

		# can not test instrumentLib directly here must test enclosed dictionary
		self.validate_library(Libraries.sweepLib.sweepDict, self.sweeps)

if __name__ == "__main__":
    unittest.main()
