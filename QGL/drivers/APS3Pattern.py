'''
Module for writing hdf5 APS3 files from sequences and patterns

Copyright 2018 Raytheon BBN Technologies

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

import struct
import io
import os
import h5py
import numpy as np
import re

from .APSPattern import unroll_loops
from QGL.Compiler import Waveform
from QGL import Compiler, ControlFlow, BlockLabel, PatternUtils
from QGL.PatternUtils import hash_pulse, flatten

MAX_WAVEFORM_VALUE = 2**15 - 1  #maximum waveform value i.e. 16bit DAC
SAMPLING_RATE = 2.5e9

# Do we want a pulse file per instrument or per channel
SEQFILE_PER_CHANNEL = False
ADDRESS_UNIT = 1

MAX_SEQS = 1024

def get_seq_file_extension():
    return '.h5'

def is_compatible_file(filename):
    with h5py.File(filename, 'r') as FID:
        target = FID['/'].attrs['target hardware']
        if isinstance(target, str):
            target = target.encode('utf-8')
        if target == b'APS3':
            return True
    return False

def get_empty_channel_set():
    return {'ch1': {}, 'ch1m1': {}}

def flatten_waveforms(n, linkList, wfLib):
    data = np.array([], dtype=np.complex128)
    for entry in linkList[n % len(linkList)]:
        if not isinstance(entry, APS3Waveform):
            continue
        if not entry.isTimeAmp:
            data = np.append(data, wfLib[wf_sig(entry)])
        else:
            data = np.append(data, wfLib[wf_sig(entry)][0] *
                                        np.ones(entry.length * entry.repeat))

    pad_width = ((len(data)-1 | 15)+1) - len(data)
    data = np.pad(data, (0,pad_width), 'constant')

    return data


def add_slave_trigger():
    pass

def wf_sig(wf):
    if wf.isZero or (wf.isTimeAmp and wf.frequency == 0):
        return (wf.amp, wf.phase)
    else:
        return (wf.key, wf.amp, round(wf.phase * 2**15), wf.length, wf.frequency)

class APS3Waveform(object):
    """
    More specific APS3 version of a waveform
    """
    def __init__(self, waveform):
        self.label        = waveform.label
        self.key          = waveform.key
        self.amp          = waveform.amp
        self.length       = PatternUtils.convert_length_to_samples(waveform.length, SAMPLING_RATE, ADDRESS_UNIT)
        self.phase        = waveform.phase
        self.frameChange  = waveform.frameChange
        self.isTimeAmp    = waveform.isTimeAmp
        self.frequency    = waveform.frequency
        self.repeat       = 1

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self.isTimeAmp:
            TA = 'HIGH' if self.amp != 0 else 'LOW'
            return "APS3Waveform-TA(" + TA + ", " + str(self.length) + ")"
        else:
            return "APS3Waveform(" + self.label + ", " + str(
                self.key)[:6] + ", " + str(self.length) + ")"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self == other

    @property
    def isZero(self):
        return self.amp == 0

TAZShape = np.zeros(1, dtype=np.complex)
TAZKey = hash_pulse(TAZShape)


def padding_entry(length):
    entry = Compiler.Waveform()
    entry.length = length / SAMPLING_RATE
    entry.key = TAZKey
    entry.isTimeAmp = True
    return APS3Waveform(entry)

def build_waveforms(seqs, shapeLib):
    # apply amplitude, phase, and modulation and add the resulting waveforms to the library
    wfLib = {wf_sig(padding_entry(0)): TAZShape}
    for wf in flatten(seqs):
        if isinstance(wf, APS3Waveform) and wf_sig(wf) not in wfLib:
            shape = np.exp(1j * wf.phase) * wf.amp * shapeLib[wf.key]
            if wf.frequency != 0 and wf.amp != 0:
                shape *= np.exp(
                    1j * 2 * np.pi * wf.frequency * np.arange(wf.length) /
                    SAMPLING_RATE)  #minus from negative frequency qubits
            wfLib[wf_sig(wf)] = shape
    return wfLib

def compress_sequences(seqs):
    '''
	Drop zero-length pulses and combine adjacent TA pairs into single entries
	'''
    for seq in seqs:
        ct = 1
        while ct < len(seq):
            prevEntry = seq[ct - 1]
            curEntry = seq[ct]
            if isinstance(curEntry, Waveform) and curEntry.length == 0:
                del seq[ct]
            elif isinstance(prevEntry, Waveform) and isinstance(curEntry, Waveform) and \
               prevEntry.isTimeAmp and curEntry.isTimeAmp and \
               prevEntry.amp == curEntry.amp and \
               prevEntry.phase == curEntry.phase:
                prevEntry.length += curEntry.length
                prevEntry.frameChange += curEntry.frameChange
                del seq[ct]
            ct += 1

def get_marker_delay(markerLL, time=False):
    '''Get the measurement marker delay. Assumptions:
        - Only one digitizer trigger.
        - The marker is at the start of the trigger pulse (length is set in APS3 firmware).
        - Marker must be aligned between sequences.
    '''
    #find the index of the first trigger pulse in each linked list
    tidx = [idx for LL in markerLL for idx, p in enumerate(LL) if getattr(p, "label", None) == "TRIG"]
    assert len(tidx) == len(markerLL), "More digitizer triggers than sequences!"
    #get the time to the end of the sequence for each linklist
    delays = [sum([p.length for p in LL[idx:]]) for LL, idx in zip(markerLL, tidx)]
    #check that all delays are equal
    assert all(x == delays[0] for x in delays), "Not all markers have the same position!"

    #convert to samples
    if time:
        return delays[0]
    else:
        #round up to nearest mutiple of 16
        return (int(delays[0] * SAMPLING_RATE)| 15 ) + 1



def preprocess(seqs, shapeLib):
    '''Pre-process waveforms for APS3 compatibility:
        - Unroll loops in miniLLs.
        - Apply frame changes to pulses.
    '''
    for seq in seqs:
        for ct,e in enumerate(seq):
            if isinstance(e, Compiler.Waveform):
                seq[ct] = APS3Waveform(e)
    seqs, miniLLrepeat = unroll_loops(seqs) #unroll_loops
    for seq in seqs:
        PatternUtils.propagate_frame_changes(seq, wf_type=Waveform)
    PatternUtils.quantize_phase(seqs, 1.0 / 2**15, wf_type=Waveform) #I think this is right...
    compress_sequences(seqs)
    wfLib = build_waveforms(seqs, shapeLib)
    #For now don't build waveforms that are too long :)
    #for ct in range(len(seqs)):
    #    seqs[ct] = apply_min_pulse_constraints(seqs[ct], wfLib)
    return seqs, miniLLrepeat, wfLib


def write_sequence_file(awgData, fileName, miniLLRepeat=1):
    '''
	Main function to pack channel LLs into an APS h5 file.
	'''
    #Preprocess the sequence data to handle APS restrictions
    LL, repeat, wfLib = preprocess(awgData['ch1']['linkList'],
                                          awgData['ch1']['wfLib'])

    if repeat != 0:
        miniLLRepeat *= repeat

    #Calculate the marker delay
    mdelay = get_marker_delay(awgData['ch1m1']['linkList'])

    nseq = len(awgData['ch1']['linkList'])

    if (nseq > MAX_SEQS):
        raise ValuError("Too many sequences! APS3 can only support {} sequences".format(MAX_SEQS))

    seq_data = []
    for n in range(nseq):
        seq_data.append(flatten_waveforms(n, LL, wfLib))

    seq_lens = np.asarray([len(data) for data in seq_data])

    #Open the HDF5 file
    if os.path.isfile(fileName):
        os.remove(fileName)
    with h5py.File(fileName, 'w') as FID:
         FID['/'].attrs['target hardware'] = 'APS3'
         FID['/'].attrs['num sequences'] = len(awgData['ch1']['linkList'])
         FID['/'].attrs['marker delay'] = np.uint32(mdelay)
         FID.create_dataset('seq_lens', shape=seq_lens.shape, data = seq_lens)

         for ct in range(len(seq_data)):
             FID.create_dataset('seq_data_{:d}'.format(ct), shape=seq_data[ct].shape,
                data=seq_data[ct], dtype=np.complex)

def read_sequence_file(fileName):
    AWGData = {'ch1I': [], 'ch1Q': [], 'ch1m1': []}

    with h5py.File(fileName, 'r') as FID:

        seq_lens = FID['seq_lens'][:]
        marker_delay = FID['/'].attrs['marker delay']

        for ct in range(len(seq_lens)):
            Idata = np.real(FID['seq_data_{:d}'.format(ct)][:])
            Qdata = np.imag(FID['seq_data_{:d}'.format(ct)][:])

            Idata = [(1,s) for s in Idata]
            Qdata = [(1,s) for s in Qdata]

            AWGData['ch1I'].append(Idata)
            AWGData['ch1Q'].append(Qdata)


            marker = np.zeros(seq_lens[ct])
            marker[seq_lens[ct]-marker_delay:seq_lens[ct]-marker_delay + 50] = 1.0
            marker = [(1,m) for m in marker]

            AWGData['ch1m1'].append(marker)

    return AWGData
