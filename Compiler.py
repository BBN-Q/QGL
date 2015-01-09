'''
functions for compiling lists of pulses/pulseBlocks down to the hardware level.

Copyright 2013 Raytheon BBN Technologies

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
import numpy as np
import os
import collections
import itertools
from warnings import warn
from copy import copy

import config
import PatternUtils
import Channels
from PulsePrimitives import Id
import PulseSequencer
import ControlFlow
from BlockLabel import label
import instruments
from instruments.AWGs import get_empty_channel_set

from APSPattern import write_APS_file
from APS2Pattern import write_APS2_file
from TekPattern import write_Tek_file
from mm import multimethod

# global parameter libraries
channelLib = {}
instrumentLib = {}

markerWFLib = {PatternUtils.TAZKey:np.zeros(1, dtype=np.bool), PatternUtils.markerHighKey:np.ones(1, dtype=np.bool) }

def get_channel_label(chanKey):
    ''' Takes in a channel key and returns a channel label '''
    if type(chanKey) != tuple:
        return chanKey.label
    else:
        return "".join([chan.label for chan in chanKey])


def setup_awg_channels(logicalChannels):
    awgs = set([])
    for chan in logicalChannels:
        awgs.add(chan.AWG)

    data = {awg.label:get_empty_channel_set(awg) for awg in awgs}
    for awgdata in data.values():
        for chan in awgdata.keys():
            awgdata[chan] = {'linkList': [], 'wfLib': {PatternUtils.TAZKey: np.zeros(1, dtype=np.complex)}}
    return data

def map_logical_to_physical(linkLists, wfLib):
    physicalChannels = {chan: channelLib[get_channel_label(chan)].physChan.label for chan in linkLists.keys()}
    awgData = setup_awg_channels(linkLists.keys())

    for chan in linkLists.keys():
        awgName, awgChan = physicalChannels[chan].split('-')
        awgData[awgName]['ch'+awgChan] = {'linkList': linkLists[chan], 'wfLib': wfLib[chan]}

    return awgData

def channel_delay_map(awgData):
    chanDelays = {}
    # loop through all used IQkeys
    for IQkey in [awgName + '-' + chanName[2:] for awgName, awg in awgData.items() for chanName in awg.keys()]:
        chan = channelLib[IQkey]
        chanDelays[IQkey] = chan.delay + chan.AWG.delay
    return PatternUtils.normalize_delays(chanDelays)

def compile_to_hardware(seqs, fileName, suffix='', alignMode="right"):
    '''
    Compiles 'seqs' to a hardware description and saves it to 'fileName'. Other inputs:
        suffix : string to append to end of fileName (e.g. with fileNames = 'test' and suffix = 'foo' might save to test-APSfoo.h5)
        alignMode : 'left' or 'right' (default 'left')
    '''

    #Add the digitizer trigger to measurements
    PatternUtils.add_digitizer_trigger(seqs, channelLib['digitizerTrig'])

    # Add gating/blanking pulses
    PatternUtils.add_gate_pulses(seqs)

    # Add the slave trigger
    PatternUtils.add_slave_trigger(seqs, channelLib['slaveTrig'])

    # find channel set at top level to account for individual sequence channel variability
    channels = set([])
    for seq in seqs:
        channels |= find_unique_channels(seq)

    #Compile all the pulses/pulseblocks to linklists and waveform libraries
    linkLists, wfLib = compile_sequences(seqs, channels)

    # apply gating constraints
    for chan, LL in linkLists.items():
        if isinstance(chan, Channels.LogicalMarkerChannel):
            linkLists[chan] = PatternUtils.apply_gating_constraints(chan.physChan, LL)

    # map logical to physical channels
    awgData = map_logical_to_physical(linkLists, wfLib)

    # construct channel delay map
    chanDelays = channel_delay_map(awgData)

    # for each physical channel need to:
    # 1) delay
    # 2) apply SSB if necessary
    # 3) mixer correct
    for awgName, data in awgData.items():
        for chanName, chanData in data.items():
            if chanData:
                # construct IQkey using existing convention
                IQkey = awgName + '-' + chanName[2:]
                chanObj = channelLib[IQkey]

                # apply channel delay
                PatternUtils.delay(chanData['linkList'], chanDelays[IQkey], chanObj.samplingRate)

                # For quadrature channels, apply SSB and mixer correction
                if isinstance(chanObj, Channels.PhysicalQuadratureChannel):

                    #At this point we finally have the timing of all the pulses so we can apply SSB
                    if hasattr(chanObj, 'SSBFreq') and abs(chanObj.SSBFreq) > 0:
                        PatternUtils.apply_SSB(chanData['linkList'], chanData['wfLib'], chanObj.SSBFreq, chanObj.samplingRate)

                    PatternUtils.correctMixer(chanData['wfLib'], chanObj.correctionT)

                #Remove unused waveforms
                compress_wfLib(chanData['linkList'], chanData['wfLib'])

    #Loop back through to fill empty channels and write to file
    fileList = []
    for awgName, data in awgData.items():
        #If all the channels are empty then do not bother writing the file
        if all([chan is None for chan in data.values()]):
            continue

        # convert to hardware formats
        # create the target folder if it does not exist
        targetFolder = os.path.split(os.path.normpath(os.path.join(config.AWGDir, fileName)))[0]
        if not os.path.exists(targetFolder):
            os.mkdir(targetFolder)
        fullFileName = os.path.normpath(os.path.join(config.AWGDir, fileName + '-' + awgName + suffix + instrumentLib[awgName].seqFileExt))

        write_sequence_file(instrumentLib[awgName], data, fullFileName)

        fileList.append(fullFileName)

    #Return the filenames we wrote
    return fileList

def compile_sequences(seqs, channels=None):
    '''
    Main function to convert sequences to miniLL's and waveform libraries.
    '''
    # all sequences should start with a WAIT
    for seq in seqs:
        if seq[0] != ControlFlow.Wait():
            seq.insert(0, ControlFlow.Wait())
    # last sequence should end with a GOTO back to the first sequence
    if not (hasattr(seqs[-1][-1], 'instruction') and seqs[-1][-1].instruction == 'GOTO'):
        seqs[-1].append(ControlFlow.Goto(label(seqs[0])))

    resolve_offsets(seqs)

    wfLib = {}
    # use seqs[0] as prototype in case we were not given a set of channels
    miniLL, wfLib = compile_sequence(seqs[0], wfLib, channels)
    linkLists = {chan: [LL] for chan, LL in miniLL.items()}
    for seq in seqs[1:]:
        miniLL, wfLib = compile_sequence(seq, wfLib, channels)
        for chan in linkLists.keys():
            linkLists[chan].append(miniLL[chan])


    #Print a message so for the experiment we know how many sequences there are
    print('Compiled {} sequences.'.format(len(seqs)))
    return linkLists, wfLib

def compile_sequence(seq, wfLib={}, channels=None):
    '''
    Converts a single sequence into a miniLL and waveform library.
    Returns a single-entry list of a miniLL and the updated wfLib
    '''

    #Find the set of logical channels used here and initialize them
    if not channels:
        channels = find_unique_channels(seq)

    logicalLLs = {}
    for chan in channels:
        logicalLLs[chan] = []
        if chan not in wfLib:
            if isinstance(chan, Channels.LogicalMarkerChannel):
                wfLib[chan] = markerWFLib
            else:
                wfLib[chan] = {PatternUtils.TAZKey: np.zeros(1, dtype=np.complex)}
    carriedPhase = {ch: 0 for ch in channels}
    for block in normalize(flatten(seq), channels):
        # control flow instructions just need to broadcast to all channels
        if isinstance(block, ControlFlow.ControlInstruction):
            for chan in channels:
                logicalLLs[chan] += [copy(block)]
            continue
        #Align the block
        blockLength = block.maxPts
        # drop length 0 blocks but push frame change onto next non-zero entry
        if blockLength == 0:
            carriedPhase = {ch: carriedPhase[ch]+block.pulses[ch].frameChange for ch in channels}
            continue
        for chan in channels:
            # add aligned LL entry(ies) (if the block contains a composite pulse, may get back multiple waveforms and LL entries)
            wfs, LLentries = align(block.label, block.pulses[chan], blockLength, block.alignment)
            for wf in wfs:
                if isinstance(chan, Channels.LogicalMarkerChannel):
                    wf = wf.astype(np.bool)
                if PatternUtils.hash_pulse(wf) not in wfLib:
                    wfLib[chan][PatternUtils.hash_pulse(wf)] = wf
            # Frame changes are then propagated through
            logicalLLs[chan] += propagate_frame(LLentries, carriedPhase[chan])
        carriedPhase = {ch: 0 for ch in channels}

    # loop through again to find phases, frame changes, and SSB modulation for quadrature channels
    for chan in channels:
        if not isinstance(chan, (Channels.Qubit, Channels.Measurement)):
            continue
        curFrame = 0
        for entry in logicalLLs[chan]:
            if isinstance(entry, ControlFlow.ControlInstruction):
                continue
            # frame update
            shape = np.copy(wfLib[chan][entry.key])

            # See if we can turn into a TA pair
            # fragile: if you buffer a square pulse it will not be constant valued
            if np.all(shape == shape[0]):
                entry.isTimeAmp = True
                shape = shape[:1]
                # convert near zeros to PatternUtils.TAZKey
                if abs(shape[0]) < 1e-6:
                    entry.key = PatternUtils.TAZKey

            #Rotate for phase and frame change (don't rotate zeros...)
            if entry.key != PatternUtils.TAZKey:
                shape *= np.exp(1j*(entry.phase+curFrame))
                shapeHash = PatternUtils.hash_pulse(shape)
                if shapeHash not in wfLib[chan]:
                    wfLib[chan][shapeHash] = shape
                entry.key = shapeHash
            curFrame += entry.frameChange

    return logicalLLs, wfLib

def propagate_frame(entries, frame):
    '''
    Propagates a frame change through a list of LL entries, dropping zero length entries
    '''
    # The first LL entry picks up the carried phase.
    entries[0].phase += frame
    entries[0].frameChange += frame
    # then push frame changes from zero length entries forward
    for prevEntry, thisEntry in zip(entries, entries[1:]):
        if prevEntry.length == 0:
            thisEntry.phase += prevEntry.frameChange
            thisEntry.frameChange += prevEntry.frameChange
    # then drop zero length entries
    for ct in reversed(range(len(entries))):
        if entries[ct].length == 0:
            del entries[ct]
    return entries

def compress_wfLib(seqs, wfLib):
    '''
    Helper function to remove unused waveforms from the library.
    '''
    usedKeys = set([PatternUtils.TAZKey, PatternUtils.markerHighKey])
    for miniLL in seqs:
        for entry in miniLL:
            if isinstance(entry, LLWaveform):
                usedKeys.add(entry.key)

    unusedKeys = set(wfLib.keys()) - usedKeys
    for key in unusedKeys:
        del wfLib[key]

def find_unique_channels(seq):
    channels = set([])
    for step in flatten(seq):
        if isinstance(step, PulseSequencer.Pulse):
            if isinstance(step.qubits, Channels.Channel):
                channels |= set([step.qubits])
            else:
                channels |= set(step.qubits)
        elif hasattr(step, 'pulses'):
            channels |= set(step.pulses.keys())
    return channels

def normalize(seq, channels=None):
    '''
    For mixed lists of Pulses and PulseBlocks, converts to list of PulseBlocks
    with uniform channels on each PulseBlock. We inject Id's where necessary.
    '''
    # promote to PulseBlocks
    seq = [p.promote() for p in seq]

    if not channels:
        channels = find_unique_channels(seq)

    # inject Id's for PulseBlocks not containing every channel
    for block in filter(lambda x: isinstance(x, PulseSequencer.PulseBlock), seq):
        emptyChannels = channels - set(block.pulses.keys())
        for ch in emptyChannels:
            block.pulses[ch] = Id(ch, length=0)
    return seq

@multimethod(instruments.AWGs.APS, dict, unicode)
def write_sequence_file(awg, data, filename):
    write_APS_file(data, filename)

@multimethod(instruments.AWGs.Tek5014, dict, unicode)
def write_sequence_file(awg, data, filename):
    write_Tek_file(data, filename, 1)

@multimethod(instruments.AWGs.APS2, dict, unicode)
def write_sequence_file(awg, data, filename):
    write_APS2_file(data, filename)

class LLWaveform(object):
    '''
    IQ LL elements for quadrature mod channels.
    '''
    def __init__(self, pulse=None, label=None):
        self.repeat = 1
        self.label = label
        if pulse is None:
            self.key = None
            self.length = 0
            self.phase = 0
            self.frameChange = 0
            self.isTimeAmp = False
            self.repeat = 1
        else:
            self.key = PatternUtils.hash_pulse(pulse.shape)
            self.length = pulse.length
            self.phase = pulse.phase
            self.frameChange = pulse.frameChange
            self.isTimeAmp = pulse.isTimeAmp

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        labelPart = "{0}: ".format(self.label) if self.label else ""
        if self.isTimeAmp:
            TA = 'HIGH' if self.key != PatternUtils.TAZKey else 'LOW'
            return labelPart + "LLWaveform-TA(" + TA + ", " + str(self.length) + ")"
        else:
            return labelPart + "LLWaveform(" + self.key[:6] + ", " + str(self.length) + ")"

    @property
    def isZero(self):
        return self.key == PatternUtils.TAZKey

    @property
    def totLength(self):
        return self.length*self.repeat

def create_padding_LL(length=0):
    tmpLL = LLWaveform()
    tmpLL.isTimeAmp = True
    tmpLL.key = PatternUtils.TAZKey
    tmpLL.length = length
    return tmpLL

def align(label, pulse, blockLength, alignment, cutoff=12):
    # check for composite pulses
    if hasattr(pulse, 'pulses'):
        entries = [LLWaveform(p) for p in pulse.pulses]
        entries[0].label = label
        shapes = [p.shape for p in pulse.pulses]
    else:
        entries = [LLWaveform(pulse, label)]
        shapes = [pulse.shape]
    padLength = blockLength - pulse.length
    if padLength == 0:
        # no padding element required
        return shapes, entries
    if (padLength < cutoff) and (alignment == "left" or alignment == "right"):
        # pad the first/last shape on one side
        if alignment == "left":
            shapes[-1] = np.hstack((shapes[-1], np.zeros(padLength)))
            entries[-1].key = PatternUtils.hash_pulse(shapes[-1])
        else: #right alignment
            shapes[0] = np.hstack((np.zeros(padLength), shapes[0]))
            entries[0].key = PatternUtils.hash_pulse(shapes[0])
    elif (padLength < 2*cutoff and alignment == "center"):
        # pad the both sides of the shape(s)
        if len(shapes) == 1:
            shapes[0] = np.hstack(( np.zeros(np.floor(padLength/2)), shapes[0], np.zeros(np.ceil(padLength/2)) ))
            entries[0].key = PatternUtils.hash_pulse(shapes[0])
        else:
            shapes[0] = np.hstack(( np.zeros(np.floor(padLength/2)), shapes[0]))
            shapes[-1] = np.hstack(( np.zeroes(np.ceil(padLength/2)), shapes[-1]))
            entries[0].key = PatternUtils.hash_pulse(shapes[0])
            entries[-1].key = PatternUtils.hash_pulse(shapes[-1])
    elif padLength == blockLength:
        #Here we have a zero-length sequence which just needs to be expanded
        entries[0].key = PatternUtils.TAZKey
        entries[0].length = blockLength
        shapes = [np.zeros(1, dtype=np.complex)]
    else:
        #split the entry into the shape and one or more TAZ
        if alignment == "left":
            padEntry = create_padding_LL(padLength)
            entries = entries + [padEntry]
        elif alignment == "right":
            padEntry = create_padding_LL(padLength)
            entries = [padEntry] + entries
        else:
            padEntry1 = create_padding_LL(np.floor(padLength/2))
            padEntry2 = create_padding_LL(np.ceil(padLength/2))
            entries = [padEntry1] + entries + [padEntry2]
    return shapes, entries

def flatten_and_separate(seq):
    '''
    Given a (potentially nested) list of instructions, flatten the list into a
    main sequence and a (flattened) list of branches
    '''
    branchSeqs = []
    stack = []
    idx = 0
    # Walk through the sequence, pruning branches and putting them on a stack
    while idx < len(seq):
        node = seq[idx]
        # 3 possibilities: have a plain node, have a nested list, or have a tuple of lists
        if isinstance(node, tuple):
            # treat the first element as the main branch and push the remaining elements on the stack
            seq[idx:idx+1] = node[0] #insert the first element into seq
            stack += list(node[1:])
        elif isinstance(node, list):
            seq[idx:idx+1] = node
        else:
            idx += 1

    while len(stack) > 0:
        node = stack.pop()
        # same 3 possibilities as above
        if isinstance(node, tuple):
            # make the first element the new top entry on the stack
            stack += list(node[1:])
            stack.append(node[0])
        elif isinstance(node, list):
            # add list elements back to stack in reverse order (so first element is on "top")
            stack += node[::-1]
        else:
            branchSeqs.append(node)

    return seq, branchSeqs

# from Stack Overflow: http://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists-in-python/2158532#2158532
def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def resolve_offsets(seqs):
    # create symbol look-up table
    symbols = {}
    for i, seq in enumerate(seqs):
        for j, entry in enumerate(seq):
            if entry.label and entry.label not in symbols:
                symbols[entry.label] = (i, j)
    # re-label targets with offsets
    for seq in seqs:
        for entry in seq:
            if hasattr(entry, 'target') and entry.target and entry.target.offset != 0:
                noOffsetLabel = copy(entry.target)
                noOffsetLabel.offset = 0
                baseidx = symbols[noOffsetLabel]
                targetidx = (baseidx[0], baseidx[1]+entry.target.offset)
                # targets are allowed to point beyond the end of the current sequence
                while targetidx[1] >= len(seqs[targetidx[0]]):
                    targetidx = (targetidx[0]+1, targetidx[1]-len(seqs[targetidx[0]]))
                    assert targetidx[0] < len(seqs), "invalid target"
                entry.target = label(seqs[targetidx[0]][targetidx[1]:])

