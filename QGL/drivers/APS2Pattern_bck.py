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

import os
import logging
from warnings import warn
from copy import copy
from future.moves.itertools import zip_longest
import pickle

import h5py
import numpy as np

from QGL import Compiler, ControlFlow, BlockLabel, PatternUtils
from QGL.PatternUtils import hash_pulse, flatten

# Python 2/3 compatibility: use 'int' that subclasses 'long'
from builtins import int

#Some constants
SAMPLING_RATE = 1.2e9
ADDRESS_UNIT = 4  #everything is done in units of 4 timesteps
MIN_ENTRY_LENGTH = 8
MAX_WAVEFORM_PTS = 2**28  #maximum size of waveform memory
WAVEFORM_CACHE_SIZE = 2**17
MAX_WAVEFORM_VALUE = 2**13 - 1  #maximum waveform value i.e. 14bit DAC
MAX_NUM_INSTRUCTIONS = 2**26
MAX_REPEAT_COUNT = 2**16 - 1
MAX_TRIGGER_COUNT = 2**32 - 1

# instruction encodings
WFM = 0x0
MARKER = 0x1
WAIT = 0x2
LOAD = 0x3
REPEAT = 0x4
CMP = 0x5
GOTO = 0x6
CALL = 0x7
RET = 0x8
SYNC = 0x9
MODULATION = 0xA
LOADCMP = 0xB
PREFETCH = 0xC
NOP = 0XF

# WFM/MARKER op codes
PLAY = 0x0
WAIT_TRIG = 0x1
WAIT_SYNC = 0x2
WFM_PREFETCH = 0x3
WFM_OP_OFFSET = 46
TA_PAIR_BIT = 45

# CMP op encodings
EQUAL = 0x0
NOTEQUAL = 0x1
GREATERTHAN = 0x2
LESSTHAN = 0x3

# Whether we use PHASE_OFFSET modulation commands or bake it into waveform
# Default to false as we usually don't have many variants
USE_PHASE_OFFSET_INSTRUCTION = False

# Whether to save the waveform offsets for partial compilation
SAVE_WF_OFFSETS = False

# Do we want a pulse file per instrument or per channel
SEQFILE_PER_CHANNEL = False

def get_empty_channel_set():
    return {'ch12': {}, 'ch12m1': {}, 'ch12m2': {}, 'ch12m3': {}, 'ch12m4': {}}


def get_seq_file_extension():
    return '.h5'


def is_compatible_file(filename):
    with h5py.File(filename, 'r') as FID:
        target = FID['/'].attrs['target hardware']
        if isinstance(target, str):
            target = target.encode('utf-8')
        if target == b'APS2':
            return True
    return False

def create_wf_vector(wfLib, seqs):
    '''
	Helper function to create the wf vector and offsets into it.
	'''
    max_pts_needed = 0
    for wf in wfLib.values():
        if len(wf) == 1:
            max_pts_needed += ADDRESS_UNIT
        else:
            max_pts_needed += len(wf)

    #If we have more than fits in cache we'll need to align and prefetch
    need_prefetch = max_pts_needed > WAVEFORM_CACHE_SIZE

    idx = 0

    if not need_prefetch:
        offsets = [{}]
        cache_lines = []
        #if we can fit them all in just pack
        wfVec = np.zeros(max_pts_needed, dtype=np.int16)
        for key, wf in wfLib.items():
            #Clip the wf
            wf[wf > 1] = 1.0
            wf[wf < -1] = -1.0
            #TA pairs need to be repeated ADDRESS_UNIT times
            if wf.size == 1:
                wf = wf.repeat(ADDRESS_UNIT)
            #Ensure the wf is an integer number of ADDRESS_UNIT's
            trim = wf.size % ADDRESS_UNIT
            if trim:
                wf = wf[:-trim]
            wfVec[idx:idx + wf.size] = np.int16(MAX_WAVEFORM_VALUE * wf)
            offsets[-1][key] = idx
            idx += wf.size

        #Trim the waveform
        wfVec.resize(idx)

    else:
        #otherwise fill in one cache line at a time
        CACHE_LINE_LENGTH = WAVEFORM_CACHE_SIZE / 2
        wfVec = np.zeros(CACHE_LINE_LENGTH, dtype=np.int16)
        offsets = [{}]
        cache_lines = []
        for seq in seqs:
            #go through sequence and see what we need to add
            pts_to_add = 0
            for entry in seq:
                if isinstance(entry, Compiler.Waveform):
                    sig = wf_sig(entry)
                    if sig not in offsets[-1]:
                        pts_to_add += entry.length

            #If what we need to add spills over then add a line and start again
            if (idx % CACHE_LINE_LENGTH) + pts_to_add > CACHE_LINE_LENGTH:
                idx = int(CACHE_LINE_LENGTH * (
                    (idx + CACHE_LINE_LENGTH) // CACHE_LINE_LENGTH))
                wfVec = np.append(wfVec,
                                  np.zeros(CACHE_LINE_LENGTH,
                                           dtype=np.int16))
                offsets.append({})

            for entry in seq:
                if isinstance(entry, Compiler.Waveform):
                    sig = wf_sig(entry)
                    if sig not in offsets[-1]:
                        wf = wfLib[sig]
                        wf[wf > 1] = 1.0
                        wf[wf < -1] = -1.0
                        #TA pairs need to be repeated ADDRESS_UNIT times
                        if wf.size == 1:
                            wf = wf.repeat(ADDRESS_UNIT)
                        #Ensure the wf is an integer number of ADDRESS_UNIT's
                        trim = wf.size % ADDRESS_UNIT
                        if trim:
                            wf = wf[:-trim]
                        wfVec[idx:idx + wf.size] = np.int16(
                            MAX_WAVEFORM_VALUE * wf)
                        offsets[-1][sig] = idx
                        idx += wf.size

            cache_lines.append(int(idx // CACHE_LINE_LENGTH))

    return wfVec, offsets, cache_lines


class Instruction(object):
    def __init__(self, header, payload, label=None, target=None):
        self.header = header
        self.payload = int(payload)
        self.label = label
        self.target = target

    @classmethod
    def unflatten(cls, instr):
        return cls(header=(int(instr) >> 56) & 0xff,
                   payload=int(instr) & 0xffffffffffffff)

    def __repr__(self):
        return self.__str__()

    def __str__(self):

        opCodes = ["WFM", "MARKER", "WAIT", "LOAD", "REPEAT", "CMP", "GOTO",
                   "CALL", "RET", "SYNC", "MODULATION", "LOADCMP", "PREFETCH",
                   "NOP", "NOP", "NOP"]

        out = "{0} ".format(self.label) if self.label else ""

        instrOpCode = (self.header >> 4) & 0xf
        out += opCodes[instrOpCode]

        if (instrOpCode == MARKER) or (instrOpCode == WFM) or (
                instrOpCode == MODULATION):
            if (instrOpCode == MARKER) or (instrOpCode == WFM):
                out += "; engine={}, ".format((self.header >> 2) & 0x3)
            else:
                out += "; "
            if self.header & 0x1:
                out += "write=1 | "
            else:
                out += "write=0 | "

        if self.target:
            out += " {}".format(self.target)

        if instrOpCode == WFM:
            wfOpCode = (self.payload >> 46) & 0x3
            wfOpCodes = ["PLAY", "TRIG", "SYNC", "PREFETCH"]
            out += wfOpCodes[wfOpCode]
            out += "; TA bit={}".format((self.payload >> 45) & 0x1)
            out += ", count = {}".format((self.payload >> 24) & 2**21 - 1)
            out += ", addr = {}".format(self.payload & 2**24 - 1)

        elif instrOpCode == MARKER:
            mrkOpCode = (self.payload >> 46) & 0x3
            mrkOpCodes = ["PLAY", "TRIG", "SYNC"]
            out += mrkOpCodes[mrkOpCode]
            out += "; state={}".format((self.payload >> 32) & 0x1)
            out += ", count = {}".format(self.payload & 2**32 - 1)

        elif instrOpCode == MODULATION:
            modulatorOpCode = (self.payload >> 45) & 0x7
            modulatorOpCodes = ["MODULATE", "RESET_PHASE", "TRIG", "SET_FREQ",
                                "SYNC", "SET_PHASE", "", "UPDATE_FRAME"]
            out += modulatorOpCodes[modulatorOpCode]
            out += "; nco_select=0x{:x}".format((self.payload >> 40) & 0xf)
            if modulatorOpCode == 0x0:
                out += ", count={:d}".format(self.payload & 0xffffffff)
            elif modulatorOpCode == 0x3:
                out += ", increment=0x{:08x}".format(self.payload & 0xffffffff)
            elif modulatorOpCode == 0x5:
                out += ", phase=0x{:08x}".format(self.payload & 0xffffffff)
            elif modulatorOpCode == 0x7:
                out += ", frame_change=0x{:08x}".format(self.payload &
                                                        0xffffffff)

        elif instrOpCode == CMP:
            cmpCodes = ["EQUAL", "NOTEQUAL", "GREATERTHAN", "LESSTHAN"]
            cmpCode = (self.payload >> 8) & 0x3
            out += " | " + cmpCodes[cmpCode]
            out += ", value = {}".format(self.payload & 0xff)

        elif any(
            [instrOpCode == op for op in [GOTO, CALL, RET, REPEAT, PREFETCH]]):
            out += " | target addr = {}".format(self.payload & 2**26 - 1)

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
        return self.payload & 0xffffffff  # bottom 32-bits of payload

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
    header = (WFM << 4) | (0x3 << 2) | (write &
                                        0x1)  #broadcast to both engines
    count = int(count)
    count = ((count // ADDRESS_UNIT) - 1) & 0x000fffff  # 20 bit count
    addr = (addr // ADDRESS_UNIT) & 0x00ffffff  # 24 bit addr
    payload = (PLAY << WFM_OP_OFFSET) | ((int(isTA) & 0x1)
                                         << TA_PAIR_BIT) | (count << 24) | addr
    return Instruction(header, payload, label)


def WaveformPrefetch(addr):
    header = (WFM << 4) | (0x3 << 2) | (0x1)
    payload = (WFM_PREFETCH << WFM_OP_OFFSET) | addr
    return Instruction(header, payload, None)


def Marker(sel, state, count, write=False, label=None):
    header = (MARKER << 4) | ((sel & 0x3) << 2) | (write & 0x1)
    count = int(count)
    four_count = ((count // ADDRESS_UNIT) - 1) & 0xffffffff  # 32 bit count
    count_rem = count % ADDRESS_UNIT
    if state == 0:
        transitionWords = {0: 0b0000, 1: 0b1000, 2: 0b1100, 3: 0b1110}
        transition = transitionWords[count_rem]
    else:
        transitionWords = {0: 0b1111, 1: 0b0111, 2: 0b0011, 3: 0b0001}
        transition = transitionWords[count_rem]
    payload = (PLAY << WFM_OP_OFFSET) | (transition << 33) | (
        (state & 0x1) << 32) | four_count
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


def Cmp(op, value, label=None):
    return Command(CMP, (op << 8) | (value & 0xff), label=label)


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


def Prefetch(addr, label=None):
    return Command(PREFETCH, addr)


def NoOp():
    return Instruction.unflatten(0xffffffffffffffff)


def preprocess(seqs, shapeLib):
    seqs = PatternUtils.convert_lengths_to_samples(
        seqs, SAMPLING_RATE, ADDRESS_UNIT, Compiler.Waveform)
    wfLib = build_waveforms(seqs, shapeLib)
    inject_modulation_cmds(seqs)
    return seqs, wfLib


def wf_sig(wf):
    '''
	Compute a signature of a Compiler.Waveform that identifies the relevant properties for
	two Waveforms to be considered "equal" in the waveform library. For example, we ignore
	length of TA waveforms.
	'''
    if wf.isZero or wf.isTimeAmp:  # 2nd condition necessary until we support RT SSB
        if USE_PHASE_OFFSET_INSTRUCTION:
            return (wf.amp)
        else:
            return (wf.amp, round(wf.phase * 2**13))
    else:
        #TODO: why do we need the length?
        if USE_PHASE_OFFSET_INSTRUCTION:
            return (wf.key, wf.amp, wf.length)
        else:
            return (wf.key, round(wf.phase * 2**13), wf.amp, wf.length)


class ModulationCommand(object):
    """docstring for ModulationCommand"""

    def __init__(self,
                 instruction,
                 nco_select,
                 frequency=0,
                 phase=0,
                 length=0):
        super(ModulationCommand, self).__init__()
        self.instruction = instruction
        self.nco_select = nco_select
        self.frequency = frequency
        self.phase = phase
        self.length = length

    def __str__(self):
        out = "Modulation({}, nco_select=0x{:x}".format(self.instruction,
                                                        self.nco_select)
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
        MODULATOR_OP_OFFSET = 44
        NCO_SELECT_OP_OFFSET = 40
        MODULATION_CLOCK = 300e6

        op_code_map = {"MODULATE": 0x0,
                       "RESET_PHASE": 0x2,
                       "SET_FREQ": 0x6,
                       "SET_PHASE": 0xa,
                       "UPDATE_FRAME": 0xe}
        payload = (op_code_map[self.instruction] << MODULATOR_OP_OFFSET) | (
            self.nco_select << NCO_SELECT_OP_OFFSET)
        if self.instruction == "MODULATE":
            #zero-indexed quad count
            payload |= np.uint32(self.length / ADDRESS_UNIT - 1)
        elif self.instruction == "SET_FREQ":
            # frequencies can span -2 to 2 or 0 to 4 in unsigned
            payload |= np.uint32(
                (self.frequency / MODULATION_CLOCK if self.frequency > 0 else
                 self.frequency / MODULATION_CLOCK + 4) * 2**28)
        elif (self.instruction == "SET_PHASE") | (
                self.instruction == "UPDATE_FRAME"):
            #phases can span -0.5 to 0.5 or 0 to 1 in unsigned
            payload |= np.uint32(np.mod(self.phase / (2 * np.pi), 1) * 2**28)

        instr = Instruction(MODULATION << 4, payload, label)
        instr.writeFlag = write_flag
        return instr

def inject_modulation_cmds(seqs):
	"""
	Inject modulation commands from phase, frequency and frameChange of waveforms
	in an IQ waveform sequence. Assume up to 2 NCOs for now.
	"""
	cur_freq = 0
	cur_phase = 0
	for ct,seq in enumerate(seqs):
		#check whether we have modulation commands
		freqs = np.unique([entry.frequency for entry in filter(lambda s: isinstance(s,Compiler.Waveform), seq)])
		if len(freqs) > 2:
			raise Exception("Max 2 frequencies on the same channel allowed.")
		no_freq_cmds = np.allclose(freqs, 0)
		phases = [entry.phase for entry in filter(lambda s: isinstance(s,Compiler.Waveform), seq)]
		no_phase_cmds = np.allclose(phases, 0)
		frame_changes = [entry.frameChange for entry in filter(lambda s: isinstance(s,Compiler.Waveform), seq)]
		no_frame_cmds = np.allclose(frame_changes, 0)
		no_modulation_cmds = no_freq_cmds and no_phase_cmds and no_frame_cmds

		if no_modulation_cmds:
			continue

		mod_seq = []
		pending_frame_update = False

		for entry in seq:

			#copies to avoid same object having different timestamps later
			#copy through BlockLabel
			if isinstance(entry, BlockLabel.BlockLabel):
				mod_seq.append(copy(entry))
			#mostly copy through control-flow
			elif isinstance(entry, ControlFlow.ControlInstruction):
				#heuristic to insert phase reset before trigger if we have modulation commands
				if isinstance(entry, ControlFlow.Wait):
					if not ( no_modulation_cmds and (cur_freq == 0) and (cur_phase == 0)):
						mod_seq.append(ModulationCommand("RESET_PHASE", 0x3))
						for nco_ind, freq in enumerate(freqs):
							mod_seq.append( ModulationCommand("SET_FREQ", nco_ind + 1, frequency = -freq) )
				elif isinstance(entry, ControlFlow.Return):
					cur_freq = 0 #makes sure that the frequency is set in the first sequence after the definition of subroutines
				mod_seq.append(copy(entry))
			elif isinstance(entry, Compiler.Waveform):
				if not no_modulation_cmds:
					#select nco
					nco_select = (list(freqs)).index(entry.frequency) + 1
					cur_freq = entry.frequency
					if USE_PHASE_OFFSET_INSTRUCTION and (entry.length > 0) and (cur_phase != entry.phase):
						mod_seq.append( ModulationCommand("SET_PHASE", nco_select, phase=entry.phase) )
						cur_phase = entry.phase
					#now apply modulation for count command and waveform command, if non-zero length
					if entry.length > 0:
						mod_seq.append(entry)
						# if we have a modulate waveform modulate pattern and there is no pending frame update we can append length to previous modulation command
						if (len(mod_seq) > 1) and (isinstance(mod_seq[-1], Compiler.Waveform)) and (isinstance(mod_seq[-2], ModulationCommand)) and (mod_seq[-2].instruction == "MODULATE") \
						and mod_seq[-1].frequency == freqs[mod_seq[-2].nco_select - 1] and not pending_frame_update:
							mod_seq[-2].length += entry.length
						else:
							mod_seq.append( ModulationCommand("MODULATE", nco_select, length = entry.length))
							pending_frame_update = False
					#now apply non-zero frame changes after so it is applied at end
					if entry.frameChange != 0:
						pending_frame_update = True
						#zero length frame changes (Z pulses) need to be combined with the previous frame change or injected where possible
						if entry.length == 0:
							#if the last is a frame change then we can add to the frame change
							if isinstance(mod_seq[-1], ModulationCommand) and mod_seq[-1].instruction == "UPDATE_FRAME":
								mod_seq[-1].phase += entry.frameChange
							#if last entry was pulse without frame change we add frame change
							elif (isinstance(mod_seq[-1], Compiler.Waveform)) or (mod_seq[-1].instruction == "MODULATE"):
								mod_seq.append( ModulationCommand("UPDATE_FRAME", nco_select, phase=entry.frameChange) )
							#if this is the first entry with a wait for trigger then we can inject a frame change
							#before the wait for trigger but after the RESET_PHASE
							elif isinstance(mod_seq[-1], ControlFlow.Wait):
								mod_seq.insert(-1, ModulationCommand("UPDATE_FRAME", nco_select, phase=entry.frameChange) )
							elif isinstance(mod_seq[-2], ControlFlow.Wait) and isinstance(mod_seq[-1], ModulationCommand) and mod_seq[-1].instruction == "SET_FREQ":
								mod_seq.insert(-2, ModulationCommand("UPDATE_FRAME", nco_select, phase=entry.frameChange) )
							#otherwise drop and error if frame has been defined
							else:
								raise Exception("Unable to implement zero time Z pulse")
						else:
							mod_seq.append( ModulationCommand("UPDATE_FRAME", nco_select, phase=entry.frameChange) )

		seqs[ct] = mod_seq

def build_waveforms(seqs, shapeLib):
    # apply amplitude (and optionally phase) and add the resulting waveforms to the library
    wfLib = {}
    for wf in flatten(seqs):
        if isinstance(wf, Compiler.Waveform) and wf_sig(wf) not in wfLib:
            shape = wf.amp * shapeLib[wf.key]
            if not USE_PHASE_OFFSET_INSTRUCTION:
                shape *= np.exp(1j * wf.phase)
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
    syncInstructions = [list(filter(
        lambda s: isinstance(s, ControlFlow.ControlInstruction), seq))
                        for seq in seqs if seq]

    # Add length to control-flow instructions to make accumulated time match at end of CFI.
    # Keep running tally of how much each channel has been shifted so far.
    localShift = [0 for _ in syncInstructions]
    for ct in range(len(syncInstructions[0])):
        step = [seq[ct] for seq in syncInstructions]
        endTime = max((s.startTime + shift
                       for s, shift in zip(step, localShift)))
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
	[wfSeq & modulationSeq, m1Seq, m2Seq, m3Seq, m4Seq]

	We take the strategy of greedily grabbing the next instruction that occurs in time, accross
	all	waveform and marker channels.
	'''

    # timestamp all entries before filtering (where we lose time information on control flow)
    for seq in seqs:
        timestamp_entries(seq)

    synchronize_clocks(seqs)

    # create (seq, startTime) pairs over all sequences
    timeTuples = []
    for ct, seq in enumerate(seqs):
        timeTuples += [(entry.startTime, ct) for entry in seq]
    timeTuples.sort()

    # keep track of where we are in each sequence
    indexes = np.zeros(len(seqs), dtype=np.int64)

    cmpTable = {'==': EQUAL, '!=': NOTEQUAL, '>': GREATERTHAN, '<': LESSTHAN}

    # always start with SYNC (stealing label from beginning of sequence)
    # unless it is a subroutine (using last entry as return as tell)
    label = None
    instructions = []
    for ct, seq in enumerate(seqs):
        if len(seq):
            first_non_empty = ct
            break
    if not isinstance(seqs[first_non_empty][-1], ControlFlow.Return):
        if isinstance(seqs[first_non_empty][0], BlockLabel.BlockLabel):
            label = seqs[first_non_empty][0]
            timeTuples.pop(0)
            indexes[first_non_empty] += 1
        instructions.append(Sync(label=label))
        label = None

    while len(timeTuples) > 0:
        #pop off all entries that have the same time
        entries = []
        start_time = 0
        while True:
            start_time, seq_idx = timeTuples.pop(0)
            entries.append((seqs[seq_idx][indexes[seq_idx]], seq_idx))
            indexes[seq_idx] += 1
            next_start_time = timeTuples[0][0] if len(timeTuples) > 0 else -1
            if start_time != next_start_time:
                break

        write_flags = [True] * len(entries)
        for ct, (entry, seq_idx) in enumerate(entries):
            #use first non empty sequence for control flow
            if seq_idx == first_non_empty and (
                    isinstance(entry, ControlFlow.ControlInstruction) or
                    isinstance(entry, BlockLabel.BlockLabel)):
                if isinstance(entry, BlockLabel.BlockLabel):
                    # carry label forward to next entry
                    label = entry
                    continue
                # control flow instructions
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
                    instructions.append(Load(entry.value - 1, label=label))
                elif isinstance(entry, ControlFlow.ComparisonInstruction):
                    # TODO modify Cmp operator to load from specified address
                    instructions.append(Cmp(cmpTable[entry.operator],
                                            entry.value,
                                            label=label))

                continue

            if seq_idx == 0:
                #analog - waveforms or modulation
                if isinstance(entry, Compiler.Waveform):
                    if entry.length < MIN_ENTRY_LENGTH:
                        warn("Dropping Waveform entry of length %s!" % entry.length)
                        continue
                    instructions.append(Waveform(offsets[wf_sig(entry)],
                                                 entry.length,
                                                 entry.isTimeAmp or
                                                 entry.isZero,
                                                 write=write_flags[ct],
                                                 label=label))
                elif isinstance(entry, ModulationCommand):
                    instructions.append(entry.to_instruction(
                        write_flag=write_flags[ct],
                        label=label))

            else:  # a marker engine
                if isinstance(entry, Compiler.Waveform):
                    if entry.length < MIN_ENTRY_LENGTH:
                        warn("Dropping entry!")
                        continue
                    markerSel = seq_idx - 1
                    state = not entry.isZero
                    instructions.append(Marker(markerSel,
                                               state,
                                               entry.length,
                                               write=write_flags[ct],
                                               label=label))

            #clear label
            label = None

    return instructions


def create_instr_data(seqs, offsets, cache_lines):
    '''
	Constructs the complete instruction data vector, and does basic checks for validity.

	Subroutines will be placed at least 8 cache lines away from sequences and aligned to cache line
	'''
    logger = logging.getLogger(__name__)
    logger.debug('')

    seq_instrs = []
    need_prefetch = len(cache_lines) > 0
    num_cache_lines = len(set(cache_lines))
    cache_line_changes = np.concatenate(
        ([0], np.where(np.diff(cache_lines))[0] + 1))
    for ct, seq in enumerate(zip_longest(*seqs, fillvalue=[])):
        seq_instrs.append(create_seq_instructions(
            list(seq), offsets[cache_lines[ct]]
            if need_prefetch else offsets[0]))
        #if we need wf prefetching and have moved waveform cache lines then inject prefetch for the next line
        if need_prefetch and (ct in cache_line_changes):
            next_cache_line = cache_lines[cache_line_changes[(np.where(
                ct == cache_line_changes)[0][0] + 1) % len(
                    cache_line_changes)]]
            seq_instrs[-1].insert(0, WaveformPrefetch(int(
                next_cache_line * WAVEFORM_CACHE_SIZE / 2)))
            #steal label
            seq_instrs[-1][0].label = seq_instrs[-1][1].label
            seq_instrs[-1][1].label = None

    #concatenate instructions
    instructions = []
    subroutines_start = -1
    for ct, seq in enumerate(seq_instrs):
        #Use last instruction being return as mark of start of subroutines
        if (seq[-1].header >> 4) == RET:
            subroutines_start = ct
            break
        instructions += seq

    #if we have any subroutines then group in cache lines
    if subroutines_start >= 0:
        subroutine_instrs = []
        subroutine_cache_line = {}
        CACHE_LINE_LENGTH = 128
        offset = 0
        for sub in seq_instrs[subroutines_start:]:
            #TODO for now we don't properly handle prefetching mulitple cache lines
            if len(sub) > CACHE_LINE_LENGTH:
                warnings.warn(
                    "Subroutines longer than {} instructions may not be prefetched correctly")
            #Don't unecessarily split across a cache line
            if (len(sub) + offset > CACHE_LINE_LENGTH) and (
                    len(sub) < CACHE_LINE_LENGTH):
                pad_instrs = 128 - ((offset + 128) % 128)
                subroutine_instrs += [NoOp()] * pad_instrs
                offset = 0
            if offset == 0:
                line_label = sub[0].label
            subroutine_cache_line[sub[0].label] = line_label
            subroutine_instrs += sub
            offset += len(sub) % CACHE_LINE_LENGTH
        logger.debug("Placed {} subroutines into {} cache lines".format(
            len(seq_instrs[subroutines_start:]), len(subroutine_instrs) //
            CACHE_LINE_LENGTH))
        #inject prefetch commands before waits
        wait_idx = [idx for idx, instr in enumerate(instructions)
                    if (instr.header >> 4) == WAIT] + [len(instructions)]
        instructions_with_prefetch = instructions[:wait_idx[0]]
        last_prefetch = None
        for start, stop in zip(wait_idx[:-1], wait_idx[1:]):
            call_targets = [instr.target for instr in instructions[start:stop]
                            if (instr.header >> 4) == CALL]
            needed_lines = set()
            for target in call_targets:
                needed_lines.add(subroutine_cache_line[target])
            if len(needed_lines) > 8:
                raise RuntimeError(
                    "Unable to prefetch more than 8 cache lines")
            for needed_line in needed_lines:
                if needed_line != last_prefetch:
                    instructions_with_prefetch.append(Prefetch(needed_line))
                    last_prefetch = needed_line
            instructions_with_prefetch += instructions[start:stop]

        instructions = instructions_with_prefetch
        #pad out instruction vector to ensure circular cache never loads a subroutine
        pad_instrs = 7 * 128 + (128 - ((len(instructions) + 128) % 128))
        instructions += [NoOp()] * pad_instrs

        instructions += subroutine_instrs

    #turn symbols into integers addresses
    resolve_symbols(instructions)

    assert len(instructions) < MAX_NUM_INSTRUCTIONS, \
    'Oops! too many instructions: {0}'.format(len(instructions))

    return np.fromiter((instr.flatten() for instr in instructions), np.uint64,
                       len(instructions))


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
        while idx + 1 < len(seq):
            if (isinstance(seq[idx], Compiler.Waveform) and
                    isinstance(seq[idx + 1], Compiler.Waveform) and
                    seq[idx].isZero == seq[idx + 1].isZero):

                seq[idx].length += seq[idx + 1].length
                del seq[idx + 1]
            else:
                idx += 1


def write_sequence_file(awgData, fileName):
    '''
	Main function to pack channel sequences into an APS2 h5 file.
	'''
    # Convert QGL IR into a representation that is closer to the hardware.
    awgData['ch12']['linkList'], wfLib = preprocess(
        awgData['ch12']['linkList'], awgData['ch12']['wfLib'])

    # compress marker data
    for field in ['ch12m1', 'ch12m2', 'ch12m3', 'ch12m4']:
        if 'linkList' in awgData[field].keys():
            PatternUtils.convert_lengths_to_samples(awgData[field]['linkList'],
                                                    SAMPLING_RATE, 1,
                                                    Compiler.Waveform)
            compress_marker(awgData[field]['linkList'])
        else:
            awgData[field]['linkList'] = []

    #Create the waveform vectors
    wfInfo = []
    wfInfo.append(create_wf_vector({key: wf.real
                                    for key, wf in wfLib.items()}, awgData[
                                        'ch12']['linkList']))
    wfInfo.append(create_wf_vector({key: wf.imag
                                    for key, wf in wfLib.items()}, awgData[
                                        'ch12']['linkList']))

    if SAVE_WF_OFFSETS:
        #create a set of all waveform signatures in offset dictionaries
        #we could have multiple offsets for the same pulse becuase it could
        #be repeated in multiple cache lines
        wf_sigs = set()
        for offset_dict in wfInfo[0][1]:
            wf_sigs |= set(offset_dict.keys())
        #create dictionary linking entry labels (that's what we'll have later) with offsets
        offsets = {}
        for seq in awgData['ch12']['linkList']:
            for entry in seq:
                if len(wf_sigs) == 0:
                    break
                if isinstance(entry, Compiler.Waveform):
                    sig = wf_sig(entry)
                    if sig in wf_sigs:
                        #store offsets and wavefor lib length
                        #time ampltidue entries are clamped to ADDRESS_UNIT
                        wf_length = ADDRESS_UNIT if entry.isTimeAmp else entry.length
                        offsets[entry.label] = ([_[sig] for _ in wfInfo[0][1]],
                                                wf_length)
                        wf_sigs.discard(sig)

            #break out of outer loop too
            if len(wf_sigs) == 0:
                break

        #now pickle the label=>offsets
        with open(os.path.splitext(fileName)[0] + ".offsets", "wb") as FID:
            pickle.dump(offsets, FID)

    # build instruction vector
    seq_data = [awgData[s]['linkList']
                for s in ['ch12', 'ch12m1', 'ch12m2', 'ch12m3', 'ch12m4']]
    instructions = create_instr_data(seq_data, wfInfo[0][1], wfInfo[0][2])

    #Open the HDF5 file
    if os.path.isfile(fileName):
        os.remove(fileName)
    with h5py.File(fileName, 'w') as FID:
        FID['/'].attrs['Version'] = 4.0
        FID['/'].attrs['target hardware'] = 'APS2'
        FID['/'].attrs['minimum firmware version'] = 4.0
        FID['/'].attrs['channelDataFor'] = np.uint16([1, 2])

        #Create the groups and datasets
        for chanct in range(2):
            chanStr = '/chan_{0}'.format(chanct + 1)
            chanGroup = FID.create_group(chanStr)
            #Write the waveformLib to file
            FID.create_dataset(chanStr + '/waveforms', data=wfInfo[chanct][0])

            #Write the instructions to channel 1
            if np.mod(chanct, 2) == 0:
                FID.create_dataset(chanStr + '/instructions',
                                   data=instructions)


def read_sequence_file(fileName):
    """
	Reads a HDF5 sequence file and returns a dictionary of lists.
	Dictionary keys are channel strings such as ch1, ch12m1
	Lists are or tuples of time-amplitude pairs (time, output)
	"""
    chanStrs = ['ch1', 'ch2', 'ch12m1', 'ch12m2', 'ch12m3', 'ch12m4',
                'mod_phase']
    seqs = {ch: [] for ch in chanStrs}

    def start_new_seq():
        for ct, ch in enumerate(chanStrs):
            if (ct < 2) or (ct == 6):
                #analog or modulation channel
                seqs[ch].append([])
            else:
                #marker channel
                seqs[ch].append([])

    with h5py.File(fileName, 'r') as FID:
        file_version = FID["/"].attrs["Version"]
        wf_lib = {}
        wf_lib['ch1'] = (
            1.0 /
            MAX_WAVEFORM_VALUE) * FID['/chan_1/waveforms'].value.flatten()
        wf_lib['ch2'] = (
            1.0 /
            MAX_WAVEFORM_VALUE) * FID['/chan_2/waveforms'].value.flatten()
        instructions = FID['/chan_1/instructions'].value.flatten()
        NUM_NCO = 2
        freq = np.zeros(NUM_NCO)  #radians per timestep
        phase = np.zeros(NUM_NCO)
        frame = np.zeros(NUM_NCO)
        next_freq = np.zeros(NUM_NCO)
        next_phase = np.zeros(NUM_NCO)
        next_frame = np.zeros(NUM_NCO)
        accumulated_phase = np.zeros(NUM_NCO)
        reset_flag = [False]*NUM_NCO

        for data in instructions:
            instr = Instruction.unflatten(data)

            modulator_opcode = instr.payload >> 44

            #update phases at these boundaries
            if (instr.opcode == WAIT) | (instr.opcode == SYNC) | (
                (instr.opcode) == MODULATION and (modulator_opcode == 0x0)):
                for ct in range(NUM_NCO):
                    if reset_flag[ct]:
                        #would expect this to be zero but this is first non-zero point
                        accumulated_phase[ct] = next_freq[ct] * ADDRESS_UNIT
                        reset_flag[ct] = False
                freq[:] = next_freq[:]
                phase[:] = next_phase[:]
                frame[:] = next_frame[:]

            #Assume new sequence at every WAIT
            if instr.opcode == WAIT:
                start_new_seq()

            elif instr.opcode == WFM and ((
                (instr.payload >> WFM_OP_OFFSET) & 0x3) == PLAY):
                addr = (instr.payload & 0x00ffffff) * ADDRESS_UNIT
                count = (instr.payload >> 24) & 0xfffff
                count = (count + 1) * ADDRESS_UNIT
                isTA = (instr.payload >> 45) & 0x1
                chan_select_bits = ((instr.header >> 2) & 0x1,
                                    (instr.header >> 3) & 0x1)
                #On older firmware we broadcast by default whereas on newer we respect the engine select
                for chan, select_bit in zip(('ch1', 'ch2'), chan_select_bits):
                    if (file_version < 4) or select_bit:
                        if isTA:
                            seqs[chan][-1].append((count, wf_lib[chan][addr]))
                        else:
                            for sample in wf_lib[chan][addr:addr + count]:
                                seqs[chan][-1].append((1, sample))

            elif instr.opcode == MARKER:
                chan = 'ch12m' + str(((instr.header >> 2) & 0x3) + 1)
                count = instr.payload & 0xffffffff
                count = (count + 1) * ADDRESS_UNIT
                state = (instr.payload >> 32) & 0x1
                seqs[chan][-1].append((count, state))

            elif instr.opcode == MODULATION:
                # modulator_op_code_map = {"MODULATE":0x0, "RESET_PHASE":0x2, "SET_FREQ":0x6, "SET_PHASE":0xa, "UPDATE_FRAME":0xe}
                nco_select_bits = (instr.payload >> 40) & 0xf
                if modulator_opcode == 0x0:
                    #modulate
                    count = ((instr.payload & 0xffffffff) + 1) * ADDRESS_UNIT
                    nco_select = {0b0001: 0,
                                  0b0010: 1,
                                  0b0100: 2,
                                  0b1000: 3}[nco_select_bits]
                    seqs['mod_phase'][-1] = np.append(
                        seqs['mod_phase'][-1], freq[nco_select] *
                        np.arange(count) + accumulated_phase[nco_select] +
                        phase[nco_select] + frame[nco_select])
                    accumulated_phase += count * freq
                else:
                    phase_rad = 2 * np.pi * (instr.payload &
                                             0xffffffff) / 2**28
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
        for ct, (
                ch1, ch2, mod_phase
        ) in enumerate(zip(seqs['ch1'], seqs['ch2'], seqs['mod_phase'])):
            if len(mod_phase):
                #only really works if ch1, ch2 are broadcast together
                mod_ch1 = []
                mod_ch2 = []
                cum_time = 0
                for ((time_ch1, amp_ch1),
                     (time_ch2, amp_ch2)) in zip(ch1, ch2):
                    if (amp_ch1 != 0) or (amp_ch2 != 0):
                        assert time_ch1 == time_ch2
                        if time_ch1 == 1:
                            #single timestep
                            modulated = np.exp(1j * mod_phase[cum_time]) * (
                                amp_ch1 + 1j * amp_ch2)
                            mod_ch1.append((1, modulated.real))
                            mod_ch2.append((1, modulated.imag))
                        else:
                            #expand TA
                            modulated = np.exp(
                                1j *
                                mod_phase[cum_time:cum_time + time_ch1]) * (
                                    amp_ch1 + 1j * amp_ch2)
                            for val in modulated:
                                mod_ch1.append((1, val.real))
                                mod_ch2.append((1, val.imag))
                    else:
                        mod_ch1.append((time_ch1, amp_ch1))
                        mod_ch2.append((time_ch2, amp_ch2))

                    cum_time += time_ch1
                seqs['ch1'][ct] = mod_ch1
                seqs['ch2'][ct] = mod_ch2

        del seqs['mod_phase']

    return seqs


def update_wf_library(filename, pulses, offsets):
    """
	Update a H5 waveform library in place give an iterable of (pulseName, pulse)
	tuples and offsets into the waveform library.
	"""
    assert USE_PHASE_OFFSET_INSTRUCTION == False
    #load the h5 file
    with h5py.File(filename) as FID:
        for label, pulse in pulses.items():
            #create a new waveform
            if pulse.isTimeAmp:
                shape = np.repeat(pulse.amp * np.exp(1j * pulse.phase), 4)
            else:
                shape = pulse.amp * np.exp(1j * pulse.phase) * pulse.shape
            try:
                length = offsets[label][1]
            except KeyError:
                print("\t{} not found in offsets so skipping".format(pulse))
                continue
            for offset in offsets[label][0]:
                print("\tUpdating {} at offset {}".format(pulse, offset))
                FID['/chan_1/waveforms'][offset:offset + length] = np.int16(
                    MAX_WAVEFORM_VALUE * shape.real)
                FID['/chan_2/waveforms'][offset:offset + length] = np.int16(
                    MAX_WAVEFORM_VALUE * shape.imag)
