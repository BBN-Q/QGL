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
from itertools import zip_longest
import pickle

import struct
import sys
import numpy as np

from QGL import Compiler, ControlFlow, BlockLabel, PatternUtils
from QGL import PulseSequencer
from QGL.PatternUtils import hash_pulse, flatten
from QGL import TdmInstructions

from .APS2Pattern import *

# Python 2/3 compatibility: use 'int' that subclasses 'long'
from builtins import int

logger = logging.getLogger(__name__)

#Some constants
SAMPLING_RATE = 5e9
MAX_WAVEFORM_VALUE = 2**15 - 1  #maximum waveform value i.e. 16bit DAC
MODULATION_CLOCK = 312.5e6

def get_seq_file_extension():
    return '.aps3'

def is_compatible_file(filename):
    with open(filename, 'rb') as FID:
        byte = FID.read(4)
        if byte == b'APS3':
            return True
    return False

def preprocess(seqs, shapeLib):
    seqs = PatternUtils.convert_lengths_to_samples(
        seqs, SAMPLING_RATE, ADDRESS_UNIT, Compiler.Waveform)
    wfLib = build_waveforms(seqs, shapeLib)
    inject_modulation_cmds(seqs)
    return seqs, wfLib

def read_sequence_file(fileName):
    """
    Reads a HDF5 sequence file and returns a dictionary of lists.
    Dictionary keys are channel strings such as ch1, m1
    Lists are or tuples of time-amplitude pairs (time, output)
    """
    chanStrs = ['ch1', 'ch2', 'm1', 'mod_phase']
    seqs = {ch: [] for ch in chanStrs}

    def start_new_seq():
        for ct, ch in enumerate(chanStrs):
            if (ct < 2) or (ct == 6):
                #analog or modulation channel
                seqs[ch].append([])
            else:
                #marker channel
                seqs[ch].append([])

    with open(fileName, 'rb') as FID:
        target_hw    = FID.read(4).decode('utf-8')
        file_version = struct.unpack('<f', FID.read(4))[0]
        min_fw       = struct.unpack('<f', FID.read(4))[0]
        num_chans    = struct.unpack('<H', FID.read(2))[0]

        inst_len     = struct.unpack('<Q', FID.read(8))[0]
        instructions = np.frombuffer(FID.read(8*inst_len), dtype=np.uint64)

        wf_lib = {}
        for i in range(num_chans):
            wf_len  = struct.unpack('<Q', FID.read(8))[0]
            wf_dat  = np.frombuffer(FID.read(2*wf_len), dtype=np.int16)
            wf_lib[f'ch{i+1}'] = ( 1.0 / MAX_WAVEFORM_VALUE) * wf_dat.flatten()

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
                chan = 'm' + str(((instr.header >> 2) & 0x3) + 1)
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

def read_waveforms(filename):
    with open(filename, 'rb') as FID:
        target_hw    = FID.read(4).decode('utf-8')
        file_version = struct.unpack('<f', FID.read(4))[0]
        min_fw       = struct.unpack('<f', FID.read(4))[0]
        num_chans    = struct.unpack('<H', FID.read(2))[0]

        inst_len     = struct.unpack('<Q', FID.read(8))[0]
        instructions = np.frombuffer(FID.read(8*inst_len), dtype=np.uint64)

        wf_dat = []
        for i in range(num_chans):
            wf_len  = struct.unpack('<Q', FID.read(8))[0]
            dat = ( 1.0 / MAX_WAVEFORM_VALUE) * np.frombuffer(FID.read(2*wf_len), dtype=np.int16).flatten()
            wf_dat.append(dat)
        return wf_dat

def write_sequence_file(awgData, fileName):
    '''
    Main function to pack channel sequences into an APS2 h5 file.
    '''
    # Convert QGL IR into a representation that is closer to the hardware.
    awgData['ch1']['linkList'], wfLib = preprocess(
        awgData['ch1']['linkList'], awgData['ch1']['wfLib'])

    # compress marker data
    for field in ['m1']:
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
                                        'ch1']['linkList']))
    wfInfo.append(create_wf_vector({key: wf.imag
                                    for key, wf in wfLib.items()}, awgData[
                                        'ch1']['linkList']))

    if SAVE_WF_OFFSETS:
        #create a set of all waveform signatures in offset dictionaries
        #we could have multiple offsets for the same pulse becuase it could
        #be repeated in multiple cache lines
        wf_sigs = set()
        for offset_dict in wfInfo[0][1]:
            wf_sigs |= set(offset_dict.keys())
        #create dictionary linking entry labels (that's what we'll have later) with offsets
        offsets = {}
        for seq in awgData['ch1']['linkList']:
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
                for s in ['ch1', 'm1']]
    instructions = create_instr_data(seq_data, wfInfo[0][1], wfInfo[0][2])

    #Open the binary file
    if os.path.isfile(fileName):
        os.remove(fileName)

    with open(fileName, 'wb') as FID:
        FID.write(b'APS3')                     # target hardware
        FID.write(np.float32(4.0).tobytes())   # Version
        FID.write(np.float32(4.0).tobytes())   # minimum firmware version
        FID.write(np.uint16(2).tobytes())      # number of channels
        # FID.write(np.uint16([1, 2]).tobytes()) # channelDataFor
        FID.write(np.uint64(instructions.size).tobytes()) # instructions length
        FID.write(instructions.tobytes()) # instructions in uint64 form

        #Create the groups and datasets
        for chanct in range(2):
            #Write the waveformLib to file
            print('Channel ' + str(chanct) + ' with size ' + str(wfInfo[chanct][0].size))
            if wfInfo[chanct][0].size == 0:
                #If there are no waveforms, ensure that there is some element
                #so that the waveform group gets written to file.
                #TODO: Fix this in libaps2
                data = np.array([0], dtype=np.int16)
            else:
                data = wfInfo[chanct][0]
            FID.write(np.uint64(data.size).tobytes()) # waveform data length for channel
            FID.write(data.tobytes())

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
        CACHE_LINE_LENGTH = int(np.round(WAVEFORM_CACHE_SIZE / 2)) - 1
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
                                  np.zeros(int(CACHE_LINE_LENGTH),
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
