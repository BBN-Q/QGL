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

# Python 2/3 compatibility: use 'int' that subclasses 'long'
from builtins import int

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
WFM        = 0x0
MARKER     = 0x1
WAIT       = 0x2
LOAD       = 0x3
REPEAT     = 0x4
CMP        = 0x5
GOTO       = 0x6
CALL       = 0x7
RET        = 0x8
SYNC       = 0x9
MODULATION = 0xA
LOADCMP    = 0xB

# WFM/MARKER op codes
PLAY          = 0x0
WAIT_TRIG     = 0x1
WAIT_SYNC     = 0x2
WFM_OP_OFFSET = 46
TA_PAIR_BIT   = 45

# CMP op encodings
EQUAL       = 0x0
NOTEQUAL    = 0x1
GREATERTHAN = 0x2
LESSTHAN    = 0x3

def get_empty_channel_set():
	return {'ch12':{}, 'ch12m1':{}, 'ch12m2':{}, 'ch12m3':{}, 'ch12m4':{}}

def get_seq_file_extension():
	return '.h5'

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
		self.payload = int(payload)
		self.label = label
		self.target = target

	@classmethod
	def unflatten(cls, instr):
		return cls(header = (int(instr) >> 56) & 0xff, payload = int(instr) & 0xffffffffffffff)

	def __repr__(self):
		return self.__str__()

	def __str__(self):

		opCodes = ["WFM", "MARKER", "WAIT", "LOAD", "REPEAT", "CMP", "GOTO", "CALL", "RET", "SYNC", "MODULATION", "LOADCMP"]

		out = "{0}: ".format(self.label) if self.label else ""

		instrOpCode = (self.header >> 4) & 0xf
		out += opCodes[instrOpCode]

		if (instrOpCode == MARKER) or (instrOpCode == WFM) or (instrOpCode == MODULATION):
			if (instrOpCode == MARKER) or (instrOpCode == WFM):
				out += "; engine={}, ".format((self.header >> 2) & 0x3)
			else:
				out += "; "
			if self.header & 0x1:
				out += "write=1 | "
			else:
				out += "write=0 | "

		if self.target:
			out += str(self.target) + "/"

		if instrOpCode == WFM:
			wfOpCode = (self.payload >> 46) & 0x3
			wfOpCodes = ["PLAY", "TRIG", "SYNC"]
			out += wfOpCodes[wfOpCode]
			out += "; TA bit={}".format((self.payload >> 45) & 0x1)
			out += ", count = {}".format((self.payload >> 24) & 2**21-1)
			out += ", addr = {}".format(self.payload & 2**24-1)

		elif instrOpCode == MARKER:
			mrkOpCode = (self.payload >> 46) & 0x3
			mrkOpCodes = ["PLAY", "TRIG", "SYNC"]
			out += mrkOpCodes[mrkOpCode]
			out += "; state={}".format((self.payload >> 32) & 0x1)
			out += ", count = {}".format(self.payload & 2**32-1)

		elif instrOpCode == MODULATION:
			modulatorOpCode = (self.payload >> 45) & 0x7
			modulatorOpCodes = ["MODULATE", "RESET_PHASE", "TRIG", "SET_FREQ", "SYNC", "SET_PHASE", "", "UPDATE_FRAME"]
			out += modulatorOpCodes[modulatorOpCode]
			out += "; nco_select=0x{:x}".format((self.payload >> 40) & 0xf)
			if modulatorOpCode == 0x0:
				out += ", count={:d}".format(self.payload & 0xffffffff)
			elif modulatorOpCode == 0x3:
				out += ", increment=0x{:08x}".format(self.payload & 0xffffffff)
			elif modulatorOpCode == 0x5:
				out += ", phase=0x{:08x}".format(self.payload & 0xffffffff)
			elif modulatorOpCode == 0x7:
				out += ", frame_change=0x{:08x}".format(self.payload & 0xffffffff)

		elif instrOpCode == CMP:
			cmpCodes = ["EQUAL", "NOTEQUAL", "GREATERTHAN", "LESSTHAN"]
			cmpCode = (self.payload >> 8) & 0x3
			out += " | " + cmpCodes[cmpCode]
			out += ", mask = {}".format(self.payload & 0xff)

		elif (instrOpCode == GOTO) or (instrOpCode == CALL) or (instrOpCode == RET) or (instrOpCode == REPEAT):
			out += " | target addr = {}".format(self.payload & 2**26-1)

		elif instrOpCode == LOAD:
			out += " | count = {}".format(self.payload)

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
		return int((self.header << 56) | (self.payload & 0xffffffffffff))

def Waveform(addr, count, isTA, write=False, label=None):
	header = (WFM << 4) | (0x3 << 2) | (write & 0x1) #broadcast to both engines
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
	if isinstance(payload, int):
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
	return Command(REPEAT, addr, label=label)

def preprocess(seqs, shapeLib):
	seqs = PatternUtils.convert_lengths_to_samples(seqs, SAMPLING_RATE, ADDRESS_UNIT)
	wfLib = build_waveforms(seqs, shapeLib)
	modulator_seqs = extract_modulation_seqs(seqs)
	return seqs, modulator_seqs, wfLib

def wf_sig(wf):
	'''
	Compute a signature of a Compiler.Waveform that identifies the relevant properties for
	two Waveforms to be considered "equal" in the waveform library. For example, we ignore
	length of TA waveforms.
	'''
	if wf.isZero or (wf.isTimeAmp and wf.frequency == 0): # 2nd condition necessary until we support RT SSB
		return (wf.amp)
	else:
		#TODO: why do we need the length?
		return (wf.key, wf.amp, wf.length)

class ModulationCommand(object):
	"""docstring for ModulationCommand"""
	def __init__(self, instruction, nco_select, frequency=0, phase=0, length=0):
		super(ModulationCommand, self).__init__()
		self.instruction = instruction
		self.nco_select = nco_select
		self.frequency = frequency
		self.phase = phase
		self.length = length

	def __str__(self):
		out = "Modulation({}, nco_select=0x{:x}".format(self.instruction, self.nco_select)
		if self.instruction == "MODULATE":
			out += ", length={})".format(self.length)
		elif self.instruction == "SET_FREQ":
			out += ", frequency={})".format(self.frequency)
		elif self.instruction == "SET_PHASE" or self.instruction == "UPDATE_FRAME":
			out += ", phase={})".format(self.phase)
		else:
			out += ")"
		return out

	def _repr_pretty_(self, p, cycle):
		p.text(str(self))

	def __repr__(self):
		return str(self)

	def to_instruction(self, write_flag=True, label=None):
		#Modulator op codes
		MODULATOR_OP_OFFSET  = 44
		NCO_SELECT_OP_OFFSET = 40
		MODULATION_CLOCK = 300e6

		op_code_map = {"MODULATE":0x0, "RESET_PHASE":0x2, "SET_FREQ":0x6, "SET_PHASE":0xa, "UPDATE_FRAME":0xe}
		payload = (op_code_map[self.instruction] << MODULATOR_OP_OFFSET) | (self.nco_select << NCO_SELECT_OP_OFFSET);
		if self.instruction == "MODULATE":
			payload |= (int(self.length/ADDRESS_UNIT -1) & 0xffffffff) #zero-indexed quad count
		elif self.instruction == "SET_FREQ":
			payload |= (int(self.frequency * (1/MODULATION_CLOCK) * 2**28) & 0xffffffff)
		elif (self.instruction == "SET_PHASE") | (self.instruction == "UPDATE_FRAME"):
			payload |= (int((np.mod(1/(2*np.pi) * self.phase + 0.5, 1) - 0.5) * 2**28) & 0xffffffff)

		instr = Instruction(MODULATION << 4, payload, label)
		instr.writeFlag = write_flag
		return instr

def extract_modulation_seqs(seqs):
	"""
	Extract modulation commands from phase, frequency and frameChange of waveforms
	in an IQ waveform sequence. Assume single NCO for now.
	"""
	modulator_seqs = []
	cur_freq = 0
	cur_phase = 0
	for seq in seqs:
		modulator_seq = []
		#check whether we have modulation commands
		freqs = np.unique([entry.frequency for entry in filter(lambda s: isinstance(s,Compiler.Waveform), seq)])
		no_freq_cmds = np.allclose(freqs, 0)
		phases = [entry.phase for entry in filter(lambda s: isinstance(s,Compiler.Waveform), seq)]
		no_phase_cmds = np.allclose(phases, 0)
		frame_changes = [entry.frameChange for entry in filter(lambda s: isinstance(s,Compiler.Waveform), seq)]
		no_frame_cmds = np.allclose(frame_changes, 0)
		no_modulation_cmds = no_freq_cmds and no_phase_cmds and no_frame_cmds
		for entry in seq:
			#copies to avoid same object having different timestamps later
			#copy through BlockLabel
			if isinstance(entry, BlockLabel.BlockLabel):
				modulator_seq.append(copy(entry))
			#mostly copy through control-flow
			elif isinstance(entry, ControlFlow.ControlInstruction):
				#heuristic to insert phase reset before trigger if we have modulation commands
				if isinstance(entry, ControlFlow.Wait):
					if not ( no_modulation_cmds and (cur_freq == 0) and (cur_phase == 0)):
						modulator_seq.append(ModulationCommand("RESET_PHASE", 0x1))
				elif isinstance(entry, ControlFlow.Return):
					cur_freq = 0 #makes sure that the frequency is set in the first sequence after the definition of subroutines
				modulator_seq.append(copy(entry))
			elif isinstance(entry, Compiler.Waveform):
				if not no_modulation_cmds:
					#insert phase and frequency commands before modulation command so they have the same start time
					phase_freq_cmds = []
					if cur_freq != entry.frequency:
						phase_freq_cmds.append( ModulationCommand("SET_FREQ", 0x1, frequency=entry.frequency) )
						cur_freq = entry.frequency
					if cur_phase != entry.phase:
						phase_freq_cmds.append( ModulationCommand("SET_PHASE", 0x1, phase=entry.phase) )
						cur_phase = entry.phase
					for cmd in phase_freq_cmds:
						modulator_seq.append(cmd)
					#now apply modulation for count command
					if (len(modulator_seq) > 0) and (isinstance(modulator_seq[-1], ModulationCommand)) and (modulator_seq[-1].instruction == "MODULATE"):
						modulator_seq[-1].length += entry.length
					else:
						modulator_seq.append( ModulationCommand("MODULATE", 0x1, length=entry.length))
					#now apply non-zero frame changes after so it is applied at end
					if entry.frameChange != 0:
						modulator_seq.append( ModulationCommand("UPDATE_FRAME", 0x1, phase=entry.frameChange) )
		modulator_seqs.append(modulator_seq)
	return modulator_seqs

def build_waveforms(seqs, shapeLib):
	# apply amplitude, and add the resulting waveforms to the library
	wfLib = {}
	for wf in flatten(seqs):
		if isinstance(wf, Compiler.Waveform) and wf_sig(wf) not in wfLib:
			shape = wf.amp * shapeLib[wf.key]
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
	syncInstructions = [list(filter(lambda s: isinstance(s, ControlFlow.ControlInstruction), seq)) for seq in seqs if seq]

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
	[wfSeq, modulationSeq, m1Seq, m2Seq, m3Seq, m4Seq]

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
	controlInstrs = list(filter(lambda s: isinstance(s, (ControlFlow.ControlInstruction, BlockLabel.BlockLabel)),
		                        seqs[ct]))
	for ct in range(len(seqs)):
		if seqs[ct]:
			seqs[ct] = list(filter(lambda s: isinstance(s, (Compiler.Waveform, ModulationCommand)), seqs[ct]))

	seqs.insert(0, controlInstrs)

	# create (seq, startTime) pairs over all sequences
	timeTuples = []
	for ct, seq in enumerate(seqs):
		timeTuples += [(entry.startTime, ct) for entry in seq]
	timeTuples.sort()

	# keep track of where we are in each sequence
	indexes = np.zeros(len(seqs), dtype=np.int64)

	cmpTable = {'==': EQUAL, '!=': NOTEQUAL, '>': GREATERTHAN, '<': LESSTHAN}

	# always start with SYNC (stealing label from beginning of sequence)
	if isinstance(seqs[0][0], BlockLabel.BlockLabel):
		label = seqs[0][0]
		timeTuples.pop(0)
		indexes[0] += 1
	else:
		label = None
	instructions = [Sync(label=label)]
	label = None

	while len(timeTuples) > 0:
		#pop off all entries that have the same time
		entries = []
		start_time = 0
		while True:
			start_time, seq_idx = timeTuples.pop(0)
			entries.append( (seqs[seq_idx][indexes[seq_idx]], seq_idx) )
			indexes[seq_idx] += 1
			next_start_time = timeTuples[0][0] if len(timeTuples) > 0 else -1
			if start_time != next_start_time:
				break

		# potentially reorder
		# heuristics:
		# 1. modulation phase updates should happen before SYNC/TRIG but in between SYNC and TRIG if we have both
		# 2. modulation phase updates should happen before NCO select commands
		# 3. SET_PHASE should happen after RESET_PHASE
		# 4. instructions to different engines should have single write flag
		def find_and_pop_entries(predicate):
			matched = []
			for ct, entry in enumerate(entries):
				if predicate(entry):
					matched.append(entries.pop(ct))
			return matched

		if len(entries) > 1:
			sync_entry = find_and_pop_entries(lambda e: isinstance(e[0], ControlFlow.Sync))
			trig_entry = find_and_pop_entries(lambda e: isinstance(e[0], ControlFlow.Wait))
			control_flow_entries = find_and_pop_entries(lambda e: isinstance(e[0], ControlFlow.ControlInstruction))
			reset_entry = find_and_pop_entries(lambda e: isinstance(e[0], ModulationCommand) and e[0].instruction == "RESET_PHASE")
			frame_entry = find_and_pop_entries(lambda e: isinstance(e[0], ModulationCommand) and e[0].instruction == "UPDATE_FRAME")
			phase_entry = find_and_pop_entries(lambda e: isinstance(e[0], ModulationCommand) and e[0].instruction == "SET_PHASE")
			freq_entry = find_and_pop_entries(lambda e: isinstance(e[0], ModulationCommand) and e[0].instruction == "SET_FREQ")
			#SYNC and TRIG:
			reordered_entries = sync_entry + control_flow_entries + reset_entry + phase_entry + freq_entry + frame_entry + trig_entry
			write_flags = [True]*len(reordered_entries)
			for entry in entries:
				reordered_entries.append(entry)
				write_flags.append(False)
			write_flags[-1] = True
			entries = reordered_entries
		else:
			write_flags = [True]

		for ct,(entry,seq_idx) in enumerate(entries):
			if seq_idx == 0:
				#labels and control-flow
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
				elif isinstance(entry, ControlFlow.LoadRepeat):
					instructions.append(Load(entry.value-1, label=label))
				elif isinstance(entry, ControlFlow.ComparisonInstruction):
					instructions.append(Cmp(cmpTable[entry.operator], entry.mask, label=label))

			elif seq_idx == 1: # waveform engine
				if entry.length < MIN_ENTRY_LENGTH:
					continue
				instructions.append(Waveform(offsets[wf_sig(entry)],
					                         entry.length,
					                         entry.isTimeAmp or entry.isZero,
					                         write=write_flags[ct],
					                         label=label))

			elif seq_idx == 2: # modulation engine
				if entry.instruction == "MODULATE" and entry.length < MIN_ENTRY_LENGTH:
					continue
				instructions.append(entry.to_instruction(write_flag=write_flags[ct], label=label))

			else: # a marker engine
				if entry.length < MIN_ENTRY_LENGTH:
					continue
				markerSel = seq_idx - 3
				state = not entry.isZero
				instructions.append(Marker(markerSel,
					                       state,
					                       entry.length,
					                       write=write_flags[ct],
					                       label=label))

			#clear label
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

def write_sequence_file(awgData, fileName):
	'''
	Main function to pack channel sequences into an APS2 h5 file.
	'''
	# Convert QGL IR into a representation that is closer to the hardware.
	awgData['ch12']['linkList'], modulation_seqs, wfLib = preprocess(awgData['ch12']['linkList'],
	                                                awgData['ch12']['wfLib'])

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
	seq_data = [awgData[s]['linkList'] for s in ['ch12', 'ch12m1', 'ch12m2', 'ch12m3', 'ch12m4']]
	seq_data.insert(1, modulation_seqs)
	instructions = create_instr_data(seq_data, wfInfo[0][1])

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

def read_sequence_file(fileName):
	chanStrs = ['ch1', 'ch2', 'ch12m1', 'ch12m2', 'ch12m3', 'ch12m4', 'mod_phase']
	seqs = {ch: [] for ch in chanStrs}
	with h5py.File(fileName, 'r') as FID:
		ch1wf = (1.0/MAX_WAVEFORM_VALUE)*FID['/chan_1/waveforms'].value.flatten()
		ch2wf = (1.0/MAX_WAVEFORM_VALUE)*FID['/chan_2/waveforms'].value.flatten()
		instructions = FID['/chan_1/instructions'].value.flatten()
		NUM_NCO = 2
		freq = np.zeros(NUM_NCO) #radians per timestep
		phase = np.zeros(NUM_NCO)
		frame = np.zeros(NUM_NCO)
		next_freq = np.zeros(NUM_NCO)
		next_phase = np.zeros(NUM_NCO)
		next_frame = np.zeros(NUM_NCO)
		accumulated_phase = np.zeros(NUM_NCO)
		reset_flag = [False, False]

		for data in instructions:
			instr = Instruction.unflatten(data)

			modulator_opcode = instr.payload >> 44

			#update phases at these boundaries
			if (instr.opcode == WAIT) | (instr.opcode == SYNC) | ((instr.opcode) == MODULATION and (modulator_opcode == 0x0) ):
				for ct in range(NUM_NCO):
					if reset_flag[ct]:
						accumulated_phase[ct] = 0
						reset_flag[ct] = False
				freq = next_freq
				phase = next_phase
				frame = next_frame

			#Assume new sequence at every WAIT
			if instr.opcode == WAIT:
				for ch in chanStrs:
					seqs[ch].append(np.array([], dtype=np.float64))
			elif instr.opcode == WFM:
				addr = (instr.payload & 0x00ffffff) * ADDRESS_UNIT
				count = (instr.payload >> 24) & 0xfffff
				count = (count + 1) * ADDRESS_UNIT
				isTA = (instr.payload >> 45) & 0x1
				chan = 'ch12m' + str(((instr.header >> 2) & 0x3) + 1)
				ch1_select_bit = (instr.header >> 2) & 0x1
				ch2_select_bit = (instr.header >> 3) & 0x1
				if not isTA:
					if ch1_select_bit:
						seqs['ch1'][-1] = np.append( seqs['ch1'][-1], ch1wf[addr:addr + count] )
					if ch2_select_bit:
						seqs['ch2'][-1] = np.append( seqs['ch2'][-1], ch2wf[addr:addr + count] )
				else:
					if ch1_select_bit:
						seqs['ch1'][-1] = np.append( seqs['ch1'][-1], np.array([ch1wf[addr]] * count) )
					if ch2_select_bit:
						seqs['ch2'][-1] = np.append( seqs['ch2'][-1], np.array([ch2wf[addr]] * count) )
			elif instr.opcode == MARKER:
				chan = 'ch12m' + str(((instr.header >> 2) & 0x3) + 1)
				count = instr.payload & 0xffffffff
				count = (count + 1) * ADDRESS_UNIT
				state = (instr.payload >> 32) & 0x1
				seqs[chan][-1] = np.append( seqs[chan][-1], np.array([state] * count) )
			elif instr.opcode == MODULATION:
				# modulator_op_code_map = {"MODULATE":0x0, "RESET_PHASE":0x2, "SET_FREQ":0x6, "SET_PHASE":0xa, "UPDATE_FRAME":0xe}
				nco_select_bits = (instr.payload >> 40) & 0xf
				if modulator_opcode == 0x0:
					#modulate
					count = ((instr.payload & 0xffffffff) + 1) * ADDRESS_UNIT
					nco_select = {0b0001:0, 0b0010:1, 0b0100:2, 0b1000:3}[nco_select_bits]
					seqs['mod_phase'][-1] = np.append(seqs['mod_phase'][-1], freq[nco_select]*np.arange(count) + accumulated_phase[nco_select] + phase[nco_select] + frame[nco_select])
					accumulated_phase += count*freq
				else:
					phase_rad = 2*np.pi * (instr.payload & 0xffffffff) / 2**28
					for ct in range(NUM_NCO):
						if (nco_select_bits >> ct) & 0x1:
							if modulator_opcode == 0x2:
								#reset
								next_phase[ct] = 0
								next_frame[ct] = 0
								reset_flag[ct] = True
							elif modulator_opcode == 0x6:
								#set frequency
								next_freq[ct] = phase_rad / ADDRESS_UNIT
							elif modulator_opcode == 0xa:
								#set phase
								next_phase[ct] = phase_rad
							elif modulator_opcode == 0xe:
								#update frame
								next_frame[ct] += phase_rad

		#Apply modulation if we have any
		for ct, (ch1,ch2,mod_phase) in enumerate(zip(seqs['ch1'], seqs['ch2'], seqs['mod_phase'])):
			if mod_phase.size:
				modulated = np.exp(1j*mod_phase) * (ch1 + 1j*ch2)
				seqs['ch1'][ct] = np.real(modulated)
				seqs['ch2'][ct] = np.imag(modulated)
		del seqs['mod_phase']

	return seqs
