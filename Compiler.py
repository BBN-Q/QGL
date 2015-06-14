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
import operator
from warnings import warn
from copy import copy

import config
import PatternUtils
import Channels
from PulsePrimitives import Id
import PulseSequencer
import ControlFlow
import BlockLabel
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

def map_logical_to_physical(wires):
    # construct a mapping of physical channels to lists of logical channels
    # (there will be more than one logical channel if multiple logical
    # channels share a physical channel)
    physicalChannels = {}
    for logicalChan in wires.keys():
        physChan = logicalChan.physChan
        if physChan not in physicalChannels:
            physicalChannels[physChan] = [logicalChan]
        else:
            physicalChannels[physChan].append(logicalChan)

    # loop through the physical channels
    physicalWires = {}
    for physChan, logicalChan in physicalChannels.items():
        if len(logicalChan) > 1:
            physicalWires[physChan] = merge_channels(wires, logicalChan)
        else:
            physicalWires[physChan] = wires[logicalChan[0]]

    return physicalWires

def merge_channels(linkLists, wfLib, channels):
    chan = channels[0]
    newLinkList = [[] for _ in range(len(linkLists[chan]))]
    newWfLib = {}
    for ct, segment in enumerate(newLinkList):
        entryIterators = [iter(linkLists[ch][ct]) for ch in channels]
        while True:
            try:
                entries = [e.next() for e in entryIterators]
                # control flow on any channel should pass thru
                if any(isinstance(e, ControlFlow.ControlInstruction) for e in entries):
                    # for the moment require uniform control flow so that we
                    # can pull from the first channel
                    assert all(e == entries[0] for e in entries), "Non-uniform control flow"
                    segment.append(entries[0])
                    continue
                # at this point we have at least one waveform instruction
                blocklength = pull_uniform_entries(entries, entryIterators, wfLib, channels)
                newentry = copy(entries[0])
                newentry.length = blocklength
                # sum waveforms
                wfnew = reduce(operator.add, [wfLib[channel][e.key] for channel, e in zip(channels, entries)])
                newentry.key = PatternUtils.hash_pulse(wfnew)
                newWfLib[newentry.key] = wfnew
                segment.append(newentry)
            except StopIteration:
                break
    return newLinkList, newWfLib

def pull_uniform_entries(entries, entryIterators, wfLib, channels):
    '''
    Given entries from a set of logical channels (that share a physical
    channel), pull enough entries from each channel so that the total pulse
    length matches. e.g. if your entry data was like this:
        ch1 : | A1 | B1 | C1 |      D1     | ...
        ch2 : |      A2      |      B2     | ...

    and initially entries = [A1, A2], we would construct
        A1* = A1 + B1 + C1
    and update entries such that entries = [A1*, A2].
    The function returns the resulting block length.
    '''
    for ct in range(len(entries)):
        while len(entries[ct]) < max(len(e) for e in entries):
            # concatenate with following entry to make up the length difference
            try:
                nextentry = entryIterators[ct].next()
            except StopIteration:
                raise NameError("Could not find a uniform section of entries")
            entries[ct] = concatenate_entries(entries[ct], nextentry, wfLib[channels[ct]])
    return max(len(e) for e in entries)

def concatenate_entries(entry1, entry2, wfLib):
    newentry = copy(entry1)
    # TA waveforms with the same amplitude can be merged with a just length update
    # otherwise, need to concatenate the pulse shapes
    if not (entry1.isTimeAmp and entry2.isTimeAmp and entry1.key == entry2.key and entry1.phase == entry2.phase and entry1.frameChange == 0):
        # otherwise, need to expand pulses and stack them
        wf = np.hstack((np.resize(wfLib[entry1.key], len(entry1)),
                        np.resize(wfLib[entry2.key], len(entry2))))
        newentry.isTimeAmp = False
        newentry.key = PatternUtils.hash_pulse(wf)
        newentry.frameChange += entry2.frameChange
        wfLib[newentry.key] = wf
    newentry.length = len(entry1) + len(entry2)
    return newentry

def generate_waveforms(physicalWires):
    wfs = {ch : {} for ch in physicalWires.keys()}
    for ch, wire in physicalWires.items():
        for pulse in flatten(wire):
            if not isinstance(pulse, PulseSequencer.Pulse):
                continue
            if pulse.hashshape() not in wfs[ch]:
                wfs[ch][pulse.hashshape()] = pulse.shape
    return wfs

def pulses_to_waveforms(physicalWires):
    wireOuts = {ch : [] for ch in physicalWires.keys()}
    for ch, seqs in physicalWires.items():
        for seq in seqs:
            wireOuts[ch].append([])
            for pulse in seq:
                if not isinstance(pulse, PulseSequencer.Pulse):
                    wireOuts[ch][-1].append(pulse)
                else:
                    wireOuts[ch][-1].append(Waveform(pulse))
    return wireOuts

def channel_delay_map(physicalWires):
    chanDelays = {chan : chan.delay + chan.AWG.delay for chan in physicalWires.keys()}
    return PatternUtils.normalize_delays(chanDelays)

def setup_awg_channels(physicalChannels):
    awgs = set([])
    for chan in physicalChannels:
        awgs.add(chan.AWG)

    data = {awg.label:get_empty_channel_set(awg) for awg in awgs}
    for awgdata in data.values():
        for chan in awgdata.keys():
            awgdata[chan] = {'linkList': [], 'wfLib': {}}
    return data

def bundle_wires(physWires, wfs):
    awgData = setup_awg_channels(physWires.keys())
    for chan in physWires.keys():
        _, awgChan = chan.label.split('-')
        awgChan = 'ch' + awgChan
        awgData[chan.AWG.label][awgChan]['linkList'] = physWires[chan]
        awgData[chan.AWG.label][awgChan]['wfLib'] = wfs[chan]
    return awgData

def compile_to_hardware(seqs, fileName, suffix=''):
    '''
    Compiles 'seqs' to a hardware description and saves it to 'fileName'. Other inputs:
        suffix : string to append to end of fileName (e.g. with fileNames = 'test' and suffix = 'foo' might save to test-APSfoo.h5)
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

    #Compile all the pulses/pulseblocks to sequences of pulses and control flow
    wireSeqs = compile_sequences(seqs, channels)

    if not validate_linklist_channels(wireSeqs.keys()):
        print "Compile to hardware failed"
        return        

    # apply gating constraints
    for chan, seq in wireSeqs.items():
        if isinstance(chan, Channels.LogicalMarkerChannel):
            wireSeqs[chan] = PatternUtils.apply_gating_constraints(chan.physChan, seq)
    # map logical to physical channels
    physWires = map_logical_to_physical(wireSeqs)

    # generate wf library (base shapes), replacing Pulse objects with Waveforms
    wfs = generate_waveforms(physWires)

    physWires = pulses_to_waveforms(physWires)

    # construct channel delay map
    delays = channel_delay_map(physWires)

    # apply delays
    for chan, wire in physWires.items():
        PatternUtils.delay(wire, delays[chan], chan.samplingRate)

    # bundle wires on instruments
    awgData = bundle_wires(physWires, wfs)

    # convert to hardware formats
    fileList = []
    for awgName, data in awgData.items():
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
        if not isinstance(seq[0], ControlFlow.Wait):
            seq.insert(0, ControlFlow.Wait())
    # turn into a loop, by appending GOTO(0) at end of last sequence
    if not isinstance(seqs[-1][-1], ControlFlow.Goto):
        seqs[-1].append(ControlFlow.Goto(BlockLabel.label(seqs[0])))

    # use seqs[0] as prototype in case we were not given a set of channels
    wires = compile_sequence(seqs[0], channels)
    wireSeqs = {chan: [seq] for chan, seq in wires.items()}
    for seq in seqs[1:]:
        wires = compile_sequence(seq, channels)
        for chan in wireSeqs.keys():
            wireSeqs[chan].append(wires[chan])

    #Print a message so for the experiment we know how many sequences there are
    print('Compiled {} sequences.'.format(len(seqs)))
    return wireSeqs

def compile_sequence(seq, channels=None):
    '''
    Takes a list of control flow and pulses, and returns aligned blocks
    separated into individual abstract channels (wires).
    '''

    #Find the set of logical channels used here and initialize them
    if not channels:
        channels = find_unique_channels(seq)

    wires = {chan: [] for chan in channels}

    for block in normalize(flatten(seq), channels):
        # labels and control flow instructions broadcast to all channels
        if isinstance(block, (BlockLabel.BlockLabel, ControlFlow.ControlInstruction)):
            for chan in channels:
                wires[chan] += [copy(block)]
            continue
        # drop length 0 blocks but push frame change onto previous entry
        if block.length == 0:
            for chan in channels:
                if len(wires[chan]) > 0:
                    wires[chan][-1].frameChange += block.pulses[chan].frameChange
                else:
                    warn("Dropping initial frame change")
            continue
        # schedule the block
        for chan in channels:
            # add aligned Pulses (if the block contains a composite pulse, may get back multiple pulses)
            wires[chan] += schedule(chan, block.pulses[chan], block.length, block.alignment)

    return wires

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
        blocklen = block.length
        emptyChannels = channels - set(block.pulses.keys())
        for ch in emptyChannels:
            block.pulses[ch] = Id(ch, length=blocklen)
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

class Waveform(object):
    '''
    IQ LL elements for quadrature mod channels.
    '''
    def __init__(self, pulse=None, label=None):
        self.label = label
        self.frequency = 0
        if pulse is None:
            self.key = None
            self.amp = 0
            self.length = 0
            self.phase = 0
            self.frameChange = 0
            self.isTimeAmp = False
        else:
            # self.key = PatternUtils.hash_pulse(pulse.shape)
            self.key = pulse.hashshape()
            self.amp = pulse.shapeParams['amp']
            self.length = pulse.shapeParams['length']
            self.phase = pulse.phase
            self.frameChange = pulse.frameChange
            self.isTimeAmp = pulse.isTimeAmp

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        labelPart = "{0}: ".format(self.label) if self.label else ""
        if self.isTimeAmp:
            TA = 'HIGH' if self.shapeParams['amp'] != 0 else 'LOW'
            return labelPart + "Waveform-TA(" + TA + ", " + str(self.length) + ")"
        else:
            return labelPart + "Waveform(" + str(self.key)[:6] + ", " + str(self.length) + ")"

    def __len__(self):
        return self.length

    @property
    def isZero(self):
        return self.key == PatternUtils.TAZKey

def schedule(channel, pulse, blockLength, alignment):
    '''
    Converts a Pulse or a CompositePulses into an aligned sequence of Pulses
    by injecting TAPulses before and/or after such that the resulting sequence
    duration is `blockLength`.
        alignment = "left", "right", or "center"
    '''
    # make everything look like a sequence
    if isinstance(pulse, PulseSequencer.CompositePulse):
        pulses = [pulse.pulses]
    else:
        pulses = [pulse]

    padLength = blockLength - pulse.length
    if padLength == 0:
        # no padding element required
        return pulses
    elif alignment == "left":
        return pulses + [Id(channel, padLength)]
    elif alignment == "right":
        return [Id(channel, padLength)] + pulses
    else: # center
        return [Id(channel, padLength/2)] + pulses + [Id(channel, padLength/2)]

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
                entry.target = BlockLabel.label(seqs[targetidx[0]][targetidx[1]:])

def validate_linklist_channels(linklistChannels):
    errors = []
    channels = channelLib.channelDict
    for channel in linklistChannels:
        if channel.label not in channels.keys() and channel.label not in errors:
            print "{0} not found in channel library".format(repr(channel))
            errors.append(channel.label)

    if errors != []:
        return False
    return True
