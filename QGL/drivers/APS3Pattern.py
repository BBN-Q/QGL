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
            if wfInfo[chanct][0].size == 0:
                #If there are no waveforms, ensure that there is some element
                #so that the waveform group gets written to file.
                #TODO: Fix this in libaps2
                data = np.array([0], dtype=np.int16)
            else:
                data = wfInfo[chanct][0]
            FID.write(np.uint64(data.size).tobytes()) # waveform data length for channel
            FID.write(data.tobytes())
