'''
Module for writing hdf5 APS2 files from LL's and patterns

Copyright 2014 Raytheon BBN Technologies

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

import h5py
import os
import numpy as np
from warnings import warn
from itertools import chain, izip_longest
from copy import copy, deepcopy
import Compiler, ControlFlow
import APSPattern

#Some constants
ADDRESS_UNIT = 4 #everything is done in units of 4 timesteps
MIN_ENTRY_LENGTH = 8
MAX_WAVEFORM_PTS = 2**28 #maximum size of waveform memory
MAX_WAVEFORM_VALUE = 2**13-1 #maximum waveform value i.e. 14bit DAC
MAX_NUM_INSTRUCTIONS = 2**26
MAX_REPEAT_COUNT = 2**16-1;
MAX_TRIGGER_COUNT = 2**32-1

#APS command bits
WAIT_TRIG_BIT = 15
# CMP_OP        = 11 # 11-13
# OP_CODE       = 8 # 8-10
# CMP_MASK      = 0 # 0-7

# instruction encodings
WFM    = 0x0
MARKER = 0x1
WAIT   = 0x2
LOAD   = 0x3
REPEAT = 0x4
CMP    = 0x5
GOTO   = 0x6
CALL   = 0x7
RET    = 0x8
SYNC   = 0x9
PFETCH = 0xA

# WFM/MARKER op codes
PLAY      = 0x0
WAIT_TRIG = 0x1
WAIT_SYNC = 0x2
TA_PAIR_BIT = 45

# CMP op encodings
EQUAL       = 0x0
NOTEQUAL    = 0x1
GREATERTHAN = 0x2
LESSTHAN    = 0x3

def create_wf_vector(wfLib):
	'''
	Helper function to create the wf vector and offsets into it.
	'''
	max_pts_needed = 0
	for wf in wfLib.values():
		if len(wf) == 1:
			max_pts_needed += ADDRESS_UNIT
		else:
			max_pts_needed += len(wf)

	wfVec = np.zeros(max_pts_needed, dtype=np.int16)
	offsets = {}
	idx = 0
	for key, wf in wfLib.items():
		#Clip the wf
		wf[wf>1] = 1.0
		wf[wf<-1] = -1.0
		#TA pairs need to be repeated ADDRESS_UNIT times
		if wf.size == 1:
			wf = wf.repeat(ADDRESS_UNIT)
		#Ensure the wf is an integer number of ADDRESS_UNIT's 
		trim = wf.size%ADDRESS_UNIT
		if trim:
			wf = wf[:-trim]
		assert idx + wf.size < MAX_WAVEFORM_PTS, 'Oops! You have exceeded the waveform memory of the APS'
		wfVec[idx:idx+wf.size] = np.uint16(np.round(MAX_WAVEFORM_VALUE*wf))
		offsets[key] = idx
		idx += wf.size 
					
	#Trim the waveform 
	# wfVec = wfVec[0:idx]
	wfVec.resize(idx+1)

	return wfVec, offsets

def calc_marker_delay(entry):
	#The firmware cannot handle 0 delay markers so push out one clock cycle
	if entry.markerDelay1 is not None:
		if entry.markerDelay1 < ADDRESS_UNIT:
			entry.markerDelay1 = ADDRESS_UNIT
		markerDelay1 = entry.markerDelay1//ADDRESS_UNIT
	else:
		markerDelay1 = 0

	if entry.markerDelay2 is not None:
		if entry.markerDelay2 < ADDRESS_UNIT:
			entry.markerDelay2 = ADDRESS_UNIT
		markerDelay2 = entry.markerDelay2//ADDRESS_UNIT
	else:
		markerDelay2 = 0

	return markerDelay1, markerDelay2

class Instruction:
	def __init__(self, header, payload, label=None, target=None):
		self.header = header
		self.payload = payload
		self.label = label
		self.target = target

	@property
	def address(self):
		return self.payload & 0xffffffff # bottom 32-bits of payload

	@address.setter
	def address(self, value):
		self.payload |= value & 0xffffffff

	def flatten(self):
		return (self.header << 56) | self.payload

def Waveform(addr, count, isTA, write=False, label=None):
	header = (WFM << 4) | (write & 0x1)
	count = int(count)
	count = ((count / ADDRESS_UNIT)-1) & 0x00ffffff
	addr = (addr / ADDRESS_UNIT) & 0x00ffffff
	payload = PLAY << 46 | (isTA & 0x1) | (count << 24) | addr
	return Instruction(header, payload, label)

def Marker(sel, state, count, write=False, label=None):
	header = (MARKER << 4) | ((sel & 0x3) << 2) | (write & 0x1)
	count = int(count)
	four_count = ((count / ADDRESS_UNIT)-1) & 0xffffffff
	count_rem = count % ADDRESS_UNIT
	# do something with count_rem....
	transition = 0

	payload = PLAY << 46 | (transition << 33) | ((state & 0x1) << 32) | four_count
	return Instruction(header, payload, label)

def Command(cmd, payload, label=None):
	header = (cmd << 4)
	if isinstance(payload, int):
		return Instruction(header, payload, label)
	else:
		return Instruction(header, 0, label, target=payload)

def Sync(label=None):
	return Command(SYNC, 0, label=label)

def Wait(label=None):
	return Command(WAIT, 0, label=label)

def Cmp(op, mask, label=None):
	return Command(CMP, (op << 8) | (mask & 0xff), label=label)

def Goto(addr, label=None):
	return Command(GOTO, addr, label=label)

def Call(addr, label=None):
	return Command(CALL, addr, label)

def Return(label=None):
	return Command(RET, 0, label=label)

def Load(count, label=None):
	return Command(LOAD, count, label=label)

def Repeat(addr, label=None):
	return Command(REPEAT, 0, label=label)

def create_instr_data(seq, offsets):
	'''
	Helper function to create LL data vectors from a list of miniLL's and an offset dictionary
	keyed on the wf keys.
	'''

	#Preallocate the bank data and do some checking for miniLL lengths
	seqLengths = np.array([len(miniLL) for miniLL in seq])
	numEntries = sum(seqLengths)
	assert numEntries < MAX_NUM_INSTRUCTIONS, 'Oops! too many instructions: {0}'.format(numEntries)
	assert seq[-1][-1].instruction == 'GOTO', 'Link list must end with a GOTO instruction'

	cmpTable = {'==': EQUAL, '!=': NOTEQUAL, '>': GREATERTHAN, '<': LESSTHAN}

	# always start with SYNC
	instructions = [Sync(label='start')]

	for entry in chain.from_iterable(seq):
		# switch on the IR instruction type
		if isinstance(entry, Compiler.LLWaveform):
			instructions.append(Waveform(offsets[entry.key],
				                         entry.length,
				                         entry.isTimeAmp,
				                         True,
				                         label=entry.label))

		elif isinstance(entry, ControlFlow.ComparisonInstruction):
			instructions.append(Cmp(cmpTable[entry.operator], entry.mask, label=entry.label))

		elif isinstance(entry, ControlFlow.ControlInstruction):
			# zero argument commands
			if entry.instruction == 'WAIT':
				instructions.append(Wait(label=entry.label))
			elif entry.instruction == 'RETURN':
				instructions.append(Return(label=entry.label))
			# target argument commands
			elif entry.instruction == 'GOTO':
				instructions.append(Goto(entry.target, label=entry.label))
			elif entry.instruction == 'CALL':
				instructions.append(Call(entry.target, label=entry.label))
			elif entry.instruction == 'REPEAT':
				instructions.append(Call(entry.target, label=entry.label))
			# value argument commands
			elif entry.instruction == 'LOAD':
				instructions.append(Load(entry.value, label=entry.label))

	instructions.append(Goto('start'))
	print(instructions)

	resolve_symbols(instructions)

	data = [instr.flatten() for instr in instructions]
	return data

def resolve_symbols(seq):
	symbols = {}
	# create symbol look-up table
	for ct, entry in enumerate(seq):
		if entry.label:
			symbols[entry.label] = ct
	# then update
	for entry in seq:
		if entry.target and entry.target in symbols:
			entry.address = symbols[entry.target]

def merge_APS_markerData(IQLL, markerLL, markerNum):
	'''
	Helper function to merge two marker channels into an IQ channel.
	'''

	markerAttr = 'markerDelay' + str(markerNum)

	# expand link lists to the same length (copying first element of shorter one)
	for miniLL_IQ, miniLL_m in izip_longest(IQLL, markerLL):
		if not miniLL_IQ:
			IQLL.append(deepcopy(IQLL[0]))
		if not miniLL_m:
			markerLL.append(deepcopy(markerLL[0]))

	#Step through the all the miniLL's together
	for miniLL_IQ, miniLL_m in zip(IQLL, markerLL):
		#Find the cummulative length for each entry of IQ channel
		timePts = np.cumsum([0] + [entry.totLength for entry in miniLL_IQ])

		#Find the switching points of the marker channels
		switchPts = []
		curIndex = 0
		for curEntry, nextEntry in zip(miniLL_m[:-1], miniLL_m[1:]):
			curIndex += curEntry.totLength
			if curEntry.key != nextEntry.key:
				switchPts.append(curIndex)
				
		#Assume switch pts seperated by 1 point are single trigger blips
		blipPts = (np.diff(switchPts) == 1).nonzero()[0]
		for pt in blipPts[::-1]:
			del switchPts[pt+1]
		#Ensure the IQ LL is long enough to support the blips
		if switchPts:
			if max(switchPts) > timePts[-1]:
				assert miniLL_IQ[-1].isTimeAmp
				miniLL_IQ[-1].length += max(switchPts) - timePts[-1] + 4 

		#Now map onto linklist elements
		curIQIdx = 0
		trigQueue = []
		for switchPt in switchPts:
			#If the trigger count is too long we need to move to the next IQ entry
			while ((switchPt - timePts[curIQIdx]) > ADDRESS_UNIT * MAX_TRIGGER_COUNT) or (len(trigQueue) > 1):
				# update the trigger queue, dropping triggers that have played
				trigQueue = [t - miniLL_IQ[curIQIdx].length for t in trigQueue]
				trigQueue = [t for t in trigQueue if t >= 0]
				curIQIdx += 1
			#Push on the trigger count
			if switchPt - timePts[curIQIdx] <= 0:
				setattr(miniLL_IQ[curIQIdx], markerAttr, 0)
				print("Had to push marker blip out to start of next entry.")
			else:
				setattr(miniLL_IQ[curIQIdx], markerAttr, switchPt - timePts[curIQIdx])
				trigQueue.insert(0, switchPt - timePts[curIQIdx])
			# update the trigger queue
			trigQueue = [t - miniLL_IQ[curIQIdx].length for t in trigQueue]
			trigQueue = [t for t in trigQueue if t >= 0]
			curIQIdx += 1

	#Replace any remaining empty entries with None
	for miniLL_IQ in IQLL:
		for entry in miniLL_IQ:
			if not hasattr(entry, markerAttr):
				setattr(entry, markerAttr, None)

def write_APS2_file(awgData, fileName):
	'''
	Main function to pack channel LLs into an APS h5 file.
	'''

	#Preprocess the LL data to handle APS restrictions
	seqs = [APSPattern.preprocess_APS(seq, awgData['ch12']['wfLib']) for seq in awgData['ch12']['linkList']]

	#Merge the the marker data into the IQ linklists
	# merge_APS_markerData(seqs, awgData['ch1m1']['linkList'], 1)
	# merge_APS_markerData(seqs, awgData['ch1m2']['linkList'], 2)
	# merge_APS_markerData(seqs, awgData['ch2m1']['linkList'], 3)
	# merge_APS_markerData(seqs, awgData['ch2m2']['linkList'], 4)
	
	#Open the HDF5 file
	if os.path.isfile(fileName):
		os.remove(fileName)
	with h5py.File(fileName, 'w') as FID:  
		FID['/'].attrs['Version'] = 3.0
		FID['/'].attrs['channelDataFor'] = np.uint16([1,2])
   
		#Create the waveform vectors
		wfInfo = []
		wfInfo.append(create_wf_vector({key:wf.real for key,wf in awgData['ch12']['wfLib'].items()}))
		wfInfo.append(create_wf_vector({key:wf.imag for key,wf in awgData['ch12']['wfLib'].items()}))

		#Create the groups and datasets
		for chanct in range(2):
			chanStr = '/chan_{0}'.format(chanct+1)
			chanGroup = FID.create_group(chanStr)
			#Write the waveformLib to file
			FID.create_dataset(chanStr+'/waveforms', data=wfInfo[chanct][0])

			#Write the instructions to channel 1
			if np.mod(chanct,2) == 0:
				FID.create_dataset(chanStr+'/instructions', data=create_instr_data(seqs, wfInfo[chanct][1]))

def read_APS2_file(fileName):
	return {'ch1': np.zeros(1), 'ch2': np.zeros(1)}
