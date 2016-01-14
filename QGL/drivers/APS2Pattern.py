'''
Module for writing hdf5 APS2 files from sequences and patterns

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
from copy import copy
from QGL import Compiler, ControlFlow, BlockLabel, PatternUtils
from QGL.PatternUtils import hash_pulse, flatten

#Some constants
SAMPLING_RATE = 1.2e9
ADDRESS_UNIT = 4 #everything is done in units of 4 timesteps
MIN_ENTRY_LENGTH = 8
MAX_WAVEFORM_PTS = 2**28 #maximum size of waveform memory
WAVEFORM_CACHE_SIZE = 2**17
MAX_WAVEFORM_VALUE = 2**13-1 #maximum waveform value i.e. 14bit DAC
MAX_NUM_INSTRUCTIONS = 2**26
MAX_REPEAT_COUNT = 2**16-1;
MAX_TRIGGER_COUNT = 2**32-1

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
LOADCMP = 0XB

# WFM/MARKER op codes
PLAY      = 0x0
WAIT_TRIG = 0x1
WAIT_SYNC = 0x2
WFM_OP_OFFSET = 46
TA_PAIR_BIT   = 45

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
		#For now assert we fit in a single waveform cache until we get PREFETCH working.
		#assert idx + wf.size < MAX_WAVEFORM_PTS, 'Oops! You have exceeded the waveform memory of the APS'
		assert idx + wf.size < WAVEFORM_CACHE_SIZE, 'Oops! You have exceeded the waveform cache of the APS'
		wfVec[idx:idx+wf.size] = np.uint16(np.round(MAX_WAVEFORM_VALUE*wf))
		offsets[key] = idx
		idx += wf.size

	#Trim the waveform
	wfVec.resize(idx)

	return wfVec, offsets

class Instruction(object):
	def __init__(self, header, payload, label=None, target=None):
		self.header = header
		self.payload = payload
		self.label = label
		self.target = target

	@classmethod
	def unflatten(cls, instr):
		return cls(header = int(long(instr) >> 56) & 0xff, payload = long(instr) & 0xffffffffffffff)

	def __repr__(self):
		return self.__str__()

	def __str__(self):

		opCodes = ["WFM", "MARKER", "WAIT", "LOAD", "REPEAT", "CMP", "GOTO", "CALL", "RET", "SYNC", "PFETCH", "LOADCMP"]


		labelPart = "{0} ".format(self.label) if self.label else ""

		instrOpCode = (self.header >> 4) & 0xf
		out = labelPart + "Instruction(" + opCodes[instrOpCode] + '|'
		if self.header & 0x1:
			out += "WRITEFLAG=1"
		else:
			out += "WRITEFLAG=0"

		if instrOpCode == 0x1:
			out += ", ENGINESELECT={}".format((self.header >> 2) & 0x3)

		out += "; "

		if self.target:
			out += str(self.target) + "/"

		if instrOpCode == 0x0:
			wfOpCode = (self.payload >> 46) & 0x3
			wfOpCodes = ["PLAY", "TRIG", "SYNC"]
			out += wfOpCodes[wfOpCode]
			out += ", TA bit={}".format((self.payload >> 45) & 0x1)
			out += ", count = {}".format((self.payload >> 24) & 2**21-1)
			out += ", addr = {}".format(self.payload & 2**24-1)

		elif instrOpCode == 0x1:
			mrkOpCode = (self.payload >> 46) & 0x3
			mrkOpCodes = ["PLAY", "TRIG", "SYNC"]
			out += mrkOpCodes[mrkOpCode]
			out += ", state={}".format((self.payload >> 32) & 0x1)
			out += ", count = {}".format(self.payload & 2**32-1)

		elif instrOpCode == 0x5:
			cmpCodes = ["EQUAL", "NOTEQUAL", "GREATERTHAN", "LESSTHAN"]
			cmpCode = (self.payload >> 8) & 0x3
			out += ", " + cmpCodes[cmpCode]
			out += ", mask = {}".format(self.payload & 0xff)

		elif instrOpCode == 0x6 or instrOpCode == 0x7 or instrOpCode == 0x8:
			out += ", target addr = {}".format(self.payload & 2**26-1)

		out += ')'

		return out

	def __eq__(self, other):
		return self.header == other.header and self.payload == other.payload and self.label == other.label

	def __ne__(self, other):
		return not self == other

	def __hash__(self):
		return hash((self.header, self.payload, self.label))

	@property
	def address(self):
		return self.payload & 0xffffffff # bottom 32-bits of payload

	@address.setter
	def address(self, value):
		self.payload |= value & 0xffffffff

	@property
	def writeFlag(self):
		return self.header & 0x1

	@writeFlag.setter
	def writeFlag(self, value):
		self.header |= value & 0x1

	@property
	def opcode(self):
		return self.header >> 4

	def flatten(self):
		return long((self.header << 56) | self.payload)

def Waveform(addr, count, isTA, write=False, label=None):
	header = (WFM << 4) | (write & 0x1)
	count = int(count)
	count = ((count // ADDRESS_UNIT)-1) & 0x000fffff # 20 bit count
	addr = (addr // ADDRESS_UNIT) & 0x00ffffff # 24 bit addr
	payload = (PLAY << WFM_OP_OFFSET) | ((int(isTA) & 0x1) << TA_PAIR_BIT) | (count << 24) | addr
	return Instruction(header, payload, label)

def Marker(sel, state, count, write=False, label=None):
	header = (MARKER << 4) | ((sel & 0x3) << 2) | (write & 0x1)
	count = int(count)
	four_count = ((count // ADDRESS_UNIT)-1) & 0xffffffff # 32 bit count
	count_rem = count % ADDRESS_UNIT
	if state == 0:
		transitionWords = {0: 0b0000, 1: 0b1000, 2: 0b1100, 3: 0b1110}
		transition = transitionWords[count_rem]
	else:
		transitionWords = {0: 0b1111, 1: 0b0111, 2: 0b0011, 3: 0b0001}
		transition = transitionWords[count_rem]
	payload = (PLAY << WFM_OP_OFFSET) | (transition << 33) | ((state & 0x1) << 32) | four_count
	return Instruction(header, payload, label)

def Command(cmd, payload, write=False, label=None):
	header = (cmd << 4)
	if isinstance(payload, (int, long)):
		instr = Instruction(header, payload, label)
	else:
		instr = Instruction(header, 0, label, target=payload)
	instr.writeFlag = write
	return instr

def Sync(label=None):
	return Command(SYNC, WAIT_SYNC << WFM_OP_OFFSET, write=True, label=label)

def Wait(label=None):
	return Command(WAIT, WAIT_TRIG << WFM_OP_OFFSET, write=True, label=label)

def LoadCmp(label=None):
	return Command(LOADCMP, 0, label=label)

def Cmp(op, mask, label=None):
	return Command(CMP, (op << 8) | (mask & 0xff), label=label)

def Goto(addr, label=None):
	return Command(GOTO, addr, label=label)

def Call(addr, label=None):
	return Command(CALL, addr, label=label)

def Return(label=None):
	return Command(RET, 0, label=label)

def Load(count, label=None):
	return Command(LOAD, count, label=label)

def Repeat(addr, label=None):
	return Command(REPEAT, 0, label=label)

def preprocess(seqs, shapeLib, T):
	for seq in seqs:
		PatternUtils.propagate_frame_changes(seq)
	seqs = PatternUtils.convert_lengths_to_samples(seqs, SAMPLING_RATE, ADDRESS_UNIT)
	PatternUtils.quantize_phase(seqs, 1.0/2**13)
	wfLib = build_waveforms(seqs, shapeLib)
	PatternUtils.correct_mixers(wfLib, T)
	return seqs, wfLib

def wf_sig(wf):
	'''
	Compute a signature of a Compiler.Waveform that identifies the relevant properties for
	two Waveforms to be considered "equal" in the waveform library. For example, we ignore
	length of TA waveforms.
	'''
	if wf.isZero or (wf.isTimeAmp and wf.frequency == 0): # 2nd condition necessary until we support RT SSB
		return (wf.amp, wf.phase)
	else:
		return (wf.key, wf.amp, round(wf.phase * 2**13), wf.length, wf.frequency)

def build_waveforms(seqs, shapeLib):
	# apply amplitude, phase, and modulation and add the resulting waveforms to the library
	wfLib = {}
	for wf in flatten(seqs):
		if isinstance(wf, Compiler.Waveform) and wf_sig(wf) not in wfLib:
			shape = np.exp(1j*wf.phase) * wf.amp * shapeLib[wf.key]
			if wf.frequency != 0 and wf.amp != 0:
				shape *= np.exp(-1j*2*np.pi*wf.frequency*np.arange(wf.length)/SAMPLING_RATE) #minus from negative frequency qubits
			wfLib[wf_sig(wf)] = shape
	return wfLib

def timestamp_entries(seq):
	t = 0
	for ct in range(len(seq)):
		seq[ct].startTime = t
		t += seq[ct].length

def synchronize_clocks(seqs):
	# Control-flow instructions (CFIs) must occur at the same time on all channels.
	# Therefore, we need to "reset the clock" by synchronizing the accumulated
	# time at each CFI to the largest value on any channel
	syncInstructions = [filter(lambda s: isinstance(s, ControlFlow.ControlInstruction), seq) for seq in seqs if seq]

	# Add length to control-flow instructions to make accumulated time match at end of CFI.
	# Keep running tally of how much each channel has been shifted so far.
	localShift = [0 for _ in syncInstructions]
	for ct in range(len(syncInstructions[0])):
		step = [seq[ct] for seq in syncInstructions]
		endTime = max((s.startTime + shift for s, shift in zip(step, localShift)))
		for ct, s in enumerate(step):
			s.length = endTime - (s.startTime + localShift[ct])
			# localShift[ct] += endTime - (s.startTime + localShift[ct])
			# the += and the last term cancel, therefore:
			localShift[ct] = endTime - s.startTime
	# re-timestamp to propagate changes across the sequences
	for seq in seqs:
		timestamp_entries(seq)
	# then transfer the control flow "lengths" back into start times
	for seq in syncInstructions:
		for instr in seq:
			instr.startTime += instr.length
			instr.length = 0

def create_seq_instructions(seqs, offsets):
	'''
	Helper function to create instruction vector from an IR sequence and an offset dictionary
	keyed on the wf keys.

	Seqs is a list of lists containing waveform and marker data, e.g.
	[wfSeq, m1Seq, m2Seq, m3Seq, m4Seq]

	We take the strategy of greedily grabbing the next instruction that occurs in time, accross
	all	waveform and marker channels.
	'''

	# timestamp all entries before filtering (where we lose time information on control flow)
	for seq in seqs:
		timestamp_entries(seq)
	synchronize_clocks(seqs)

	# filter out sequencing instructions from the waveform and marker lists, so that seqs becomes:
	# [control-flow, wfs, m1, m2, m3, m4]
	# control instructions get broadcast so pull them from the first non-empty sequence
	try:
		ct = next(i for i,j in enumerate(seqs) if j)
	except StopIteration:
		print("No non-empty sequences to create!")
		raise
	controlInstrs = filter(lambda s: isinstance(s, (ControlFlow.ControlInstruction, BlockLabel.BlockLabel)),
		                   seqs[ct])
	for ct in range(len(seqs)):
		if seqs[ct]:
			seqs[ct] = filter(lambda s: isinstance(s, Compiler.Waveform), seqs[ct])

	seqs.insert(0, controlInstrs)

	# create (seq, startTime) pairs over all sequences
	timeTuples = []
	for ct, seq in enumerate(seqs):
		timeTuples += [(entry.startTime, ct) for entry in seq]
	timeTuples.sort()

	# keep track of where we are in each sequence
	curIdx = np.zeros(len(seqs), dtype=np.int64)

	cmpTable = {'==': EQUAL, '!=': NOTEQUAL, '>': GREATERTHAN, '<': LESSTHAN}

	# always start with SYNC (stealing label from beginning of sequence)
	if isinstance(seqs[0][0], BlockLabel.BlockLabel):
		label = seqs[0][0]
		timeTuples.pop(0)
		curIdx[0] += 1
	else:
		label = None
	instructions = [Sync(label=label)]
	label = None

	while len(timeTuples) > 0:
		startTime, curSeq = timeTuples.pop(0)
		entry = seqs[curSeq][curIdx[curSeq]]
		nextStartTime = timeTuples[0][0] if len(timeTuples) > 0 else -1
		writeFlag = (startTime != nextStartTime)
		curIdx[curSeq] += 1

		# poor man's way of deciding waveform or marker is to use curSeq
		if curSeq == 1: # waveform channel
			if entry.length < MIN_ENTRY_LENGTH:
				continue
			instructions.append(Waveform(offsets[wf_sig(entry)],
				                         entry.length,
				                         entry.isTimeAmp or entry.isZero,
				                         write=writeFlag,
				                         label=label))
		elif curSeq > 1: # a marker channel
			if entry.length < MIN_ENTRY_LENGTH:
				continue
			markerSel = curSeq - 2
			state = not entry.isZero
			instructions.append(Marker(markerSel,
				                       state,
				                       entry.length,
				                       write=writeFlag,
				                       label=label))

		else: # otherwise we are dealing with labels and control-flow
			if isinstance(entry, BlockLabel.BlockLabel):
				# carry label forward to next entry
				label = entry
				continue
			# zero argument commands
			elif isinstance(entry, ControlFlow.Wait):
				instructions.append(Wait(label=label))
			elif isinstance(entry, ControlFlow.LoadCmp):
				instructions.append(LoadCmp(label=label))
			elif isinstance(entry, ControlFlow.Sync):
				instructions.append(Sync(label=label))
			elif isinstance(entry, ControlFlow.Return):
				instructions.append(Return(label=label))
			# target argument commands
			elif isinstance(entry, ControlFlow.Goto):
				instructions.append(Goto(entry.target, label=label))
			elif isinstance(entry, ControlFlow.Call):
				instructions.append(Call(entry.target, label=label))
			elif isinstance(entry, ControlFlow.Repeat):
				instructions.append(Repeat(entry.target, label=label))
			# value argument commands
			elif isinstance(entry, ControlFlow.Load):
				instructions.append(Load(entry.value-1, label=label))
			elif isinstance(entry, ControlFlow.ComparisonInstruction):
				instructions.append(Cmp(cmpTable[entry.operator], entry.mask, label=label))
		label = None

	return instructions

def create_instr_data(seqs, offsets):
	'''
	Constructs the complete instruction data vector, and does basic checks for validity.
	'''
	maxlen = max([len(s) for s in seqs])
	instructions = []
	for ct in range(maxlen):
		instructions += create_seq_instructions([s[ct] if ct < len(s) else [] for s in seqs], offsets)

	resolve_symbols(instructions)

	if instructions[-1] != Goto(0):
		instructions.append(Goto(0))

	assert len(instructions) < MAX_NUM_INSTRUCTIONS, 'Oops! too many instructions: {0}'.format(len(instructions))
	data = np.array([instr.flatten() for instr in instructions], dtype=np.uint64)
	return data

def resolve_symbols(seq):
	symbols = {}
	# create symbol look-up table
	for ct, entry in enumerate(seq):
		if entry.label and entry.label not in symbols:
			symbols[entry.label] = ct
	# then update
	for entry in seq:
		if entry.target:
			entry.address = symbols[entry.target]

def compress_marker(markerLL):
	'''
	Compresses adjacent entries of the same state into single entries
	'''
	for seq in markerLL:
		idx = 0
		while idx+1 < len(seq):
			if (isinstance(seq[idx], Compiler.Waveform)
				and isinstance(seq[idx+1], Compiler.Waveform)
				and seq[idx].isZero == seq[idx+1].isZero):

				seq[idx].length += seq[idx+1].length
				del seq[idx+1]
			else:
				idx += 1

def write_APS2_file(awgData, fileName):
	'''
	Main function to pack channel sequences into an APS2 h5 file.
	'''
	# Convert QGL IR into a representation that is closer to the hardware.
	awgData['ch12']['linkList'], wfLib = preprocess(awgData['ch12']['linkList'],
	                                                awgData['ch12']['wfLib'],
	                                                awgData['ch12']['correctionT'])

	# compress marker data
	for field in ['ch12m1', 'ch12m2', 'ch12m3', 'ch12m4']:
		if 'linkList' in awgData[field].keys():
			PatternUtils.convert_lengths_to_samples(awgData[field]['linkList'], SAMPLING_RATE)
			compress_marker(awgData[field]['linkList'])
		else:
			awgData[field]['linkList'] = []

	#Create the waveform vectors
	wfInfo = []
	wfInfo.append(create_wf_vector({key:wf.real for key,wf in wfLib.items()}))
	wfInfo.append(create_wf_vector({key:wf.imag for key,wf in wfLib.items()}))

	# build instruction vector
	instructions = create_instr_data([awgData[s]['linkList'] for s in ['ch12', 'ch12m1', 'ch12m2', 'ch12m3', 'ch12m4']],
		wfInfo[0][1])

	#Open the HDF5 file
	if os.path.isfile(fileName):
		os.remove(fileName)
	with h5py.File(fileName, 'w') as FID:
		FID['/'].attrs['Version'] = 3.0
		FID['/'].attrs['channelDataFor'] = np.uint16([1,2])

		#Create the groups and datasets
		for chanct in range(2):
			chanStr = '/chan_{0}'.format(chanct+1)
			chanGroup = FID.create_group(chanStr)
			#Write the waveformLib to file
			FID.create_dataset(chanStr+'/waveforms', data=wfInfo[chanct][0])

			#Write the instructions to channel 1
			if np.mod(chanct,2) == 0:
				FID.create_dataset(chanStr+'/instructions', data=instructions)

def read_APS2_file(fileName):
	chanStrs = ['ch1', 'ch2', 'ch12m1', 'ch12m2', 'ch12m3', 'ch12m4']
	seqs = {ch: [] for ch in chanStrs}
	with h5py.File(fileName, 'r') as FID:
		ch1wf = (1.0/MAX_WAVEFORM_VALUE)*FID['/chan_1/waveforms'].value.flatten()
		ch2wf = (1.0/MAX_WAVEFORM_VALUE)*FID['/chan_2/waveforms'].value.flatten()
		instructions = FID['/chan_1/instructions'].value.flatten()

		for data in instructions:
			instr = Instruction.unflatten(data)
			if instr.opcode == WAIT:
				for ch in chanStrs:
					seqs[ch].append(np.array([], dtype=np.float64))
			elif instr.opcode == WFM:
				addr = (instr.payload & 0x00ffffff) * ADDRESS_UNIT
				count = (instr.payload >> 24) & 0xfffff
				count = (count + 1) * ADDRESS_UNIT
				isTA = (instr.payload >> 45) & 0x1
				if not isTA:
					seqs['ch1'][-1] = np.append( seqs['ch1'][-1], ch1wf[addr:addr + count] )
					seqs['ch2'][-1] = np.append( seqs['ch2'][-1], ch2wf[addr:addr + count] )
				else:
					seqs['ch1'][-1] = np.append( seqs['ch1'][-1], np.array([ch1wf[addr]] * count) )
					seqs['ch2'][-1] = np.append( seqs['ch2'][-1], np.array([ch2wf[addr]] * count) )
			elif instr.opcode == MARKER:
				chan = 'ch12m' + str(((instr.header >> 2) & 0x3) + 1)
				count = instr.payload & 0xffffffff
				count = (count + 1) * ADDRESS_UNIT
				state = (instr.payload >> 32) & 0x1
				seqs[chan][-1] = np.append( seqs[chan][-1], np.array([state] * count) )

	return seqs
