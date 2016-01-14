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
import itertools
import operator
from warnings import warn
from copy import copy

import config
import PatternUtils
from PatternUtils import flatten
import Channels
from PulsePrimitives import Id
import PulseSequencer
import ControlFlow
import BlockLabel
import instruments


# global parameter libraries
channelLib = {}
instrumentLib = {}

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

def merge_channels(wires, channels):
    chan = channels[0]
    mergedWire = [[] for _ in range(len(wires[chan]))]
    shapeFunLib = {}
    for ct, segment in enumerate(mergedWire):
        entryIterators = [iter(wires[ch][ct]) for ch in channels]
        while True:
            try:
                entries = [e.next() for e in entryIterators]
                # control flow on any channel should pass thru
                if any(isinstance(e, (ControlFlow.ControlInstruction, BlockLabel.BlockLabel)) for e in entries):
                    # for the moment require uniform control flow so that we
                    # can pull from the first channel
                    assert all(e == entries[0] for e in entries), "Non-uniform control flow"
                    segment.append(entries[0])
                    continue
                # at this point we have at least one waveform instruction
                blocklength = pull_uniform_entries(entries, entryIterators)
                newentry = copy(entries[0])
                #TODO properly deal with constant pulses
                newentry.amp = 1.0
                newentry.isTimeAmp = all([e.isTimeAmp for e in entries])
                if all([e.amp == 0 for e in entries]):
                    newentry.amp = 0
                else:
                    assert np.count_nonzero([e.amp * e.channel.frequency for e in entries]) <= 1, "Unable to handle merging more than one non-zero entry with non-zero frequency."

                #If there is a non-zero SSB frequency copy it to the new entry
                nonZeroSSBChan = np.nonzero([e.amp * e.channel.frequency for e in entries])[0]
                if nonZeroSSBChan:
                    newentry.channel = entries[nonZeroSSBChan[0]].channel

                newentry.phase = 0

                pulsesHash = tuple([e.hashshape() for e in entries])
                if pulsesHash not in shapeFunLib:
                    # create closure to sum waveforms
                    def sum_shapes(entries=entries, **kwargs):
                        return reduce(operator.add, [e.amp * np.exp(1j*e.phase) * e.shape for e in entries])
                    shapeFunLib[pulsesHash] = sum_shapes
                newentry.shapeParams = {'shapeFun':shapeFunLib[pulsesHash], 'length':blocklength}
                newentry.label = "*".join([e.label for e in entries])
                segment.append(newentry)


            except StopIteration:
                break
    return mergedWire

def pull_uniform_entries(entries, entryIterators):
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
    numChan = len(entries)
    iterDone = [False]*numChan #keep track of how many entry iterators are used up
    ct = 0
    while True:
        #If we've used up all the entries on all the channels we're done
        if all(iterDone):
            raise StopIteration("Unable to find a uniform set of entries")

        #If all the entry lengths are the same we are finished
        entryLengths = [e.length for e in entries]
        if all(x==entryLengths[0] for x in entryLengths):
            break

        #Otherwise try to concatenate on entries to match lengths
        while entries[ct].length < max(e.length for e in entries):
            # concatenate with following entry to make up the length difference
            try:
                nextentry = entryIterators[ct].next()
            except StopIteration:
                iterDone[ct] = True

            entries[ct] = concatenate_entries(entries[ct], nextentry)

        ct = (ct + 1) % numChan

    return max(e.length for e in entries)

def concatenate_entries(entry1, entry2):
    newentry = copy(entry1)
    # TA waveforms with the same amplitude can be merged with a just length update
    # otherwise, need to concatenate the pulse shapes
    if not (entry1.isTimeAmp and entry2.isTimeAmp and entry1.amp == entry2.amp and entry1.phase == (entry1.frameChange + entry2.phase)):
        # otherwise, need to build a closure to stack them
        def stack_shapes(entry1=entry1, entry2=entry2, **kwargs):
            return np.hstack((entry1.amp * np.exp(1j*entry1.phase) * entry1.shape,
                              entry2.amp * np.exp(1j*(entry1.frameChange + entry2.phase)) * entry2.shape))

        newentry.isTimeAmp = False
        newentry.shapeParams = {'shapeFun' : stack_shapes}
        newentry.label = entry1.label + '+' + entry2.label
    newentry.frameChange += entry2.frameChange
    newentry.length = entry1.length + entry2.length
    newentry.amp = 1.0

    return newentry

def generate_waveforms(physicalWires):
    wfs = {ch : {} for ch in physicalWires.keys()}
    for ch, wire in physicalWires.items():
        for pulse in flatten(wire):
            if not isinstance(pulse, PulseSequencer.Pulse):
                continue
            if pulse.hashshape() not in wfs[ch]:
                if pulse.isTimeAmp:
                    wfs[ch][pulse.hashshape()] = np.ones(1, dtype=np.complex)
                else:
                    wfs[ch][pulse.hashshape()] = pulse.shape
    return wfs

def pulses_to_waveforms(physicalWires):
    wireOuts = {ch : [] for ch in physicalWires.keys()}
    for ch, seqs in physicalWires.iteritems():
        for seq in seqs:
            wireOuts[ch].append([])
            for pulse in seq:
                if not isinstance(pulse, PulseSequencer.Pulse):
                    wireOuts[ch][-1].append(pulse)
                else:
                    wf = Waveform(pulse)
                    wireOuts[ch][-1].append(wf)
    return wireOuts

def channel_delay_map(physicalWires):
    chanDelays = {chan : chan.delay + chan.AWG.delay for chan in physicalWires.keys()}
    return PatternUtils.normalize_delays(chanDelays)

def setup_awg_channels(physicalChannels):
    awgs = set([])
    for chan in physicalChannels:
        awgs.add(chan.AWG)

    data = {awg.label:awg.get_empty_channel_set() for awg in awgs}
    for awgdata in data.values():
        for chan in awgdata.keys():
            awgdata[chan] = {'linkList': [], 'wfLib': {}, 'correctionT': np.identity(2)}
    return data

def bundle_wires(physWires, wfs):
    awgData = setup_awg_channels(physWires.keys())
    for chan in physWires.keys():
        _, awgChan = chan.label.split('-')
        awgChan = 'ch' + awgChan
        awgData[chan.AWG.label][awgChan]['linkList'] = physWires[chan]
        awgData[chan.AWG.label][awgChan]['wfLib'] = wfs[chan]
        if hasattr(chan, 'correctionT'):
            awgData[chan.AWG.label][awgChan]['correctionT'] = chan.correctionT
    return awgData

def collect_specializations(seqs):
    '''
    Collects function definitions for all targets of Call instructions
    '''
    targets = [x.target for x in flatten(seqs) if isinstance(x, ControlFlow.Call)]
    funcDefs = []
    for target in targets:
        funcDefs += ControlFlow.qfunction_specialization(target)
    return funcDefs

def compile_to_hardware(seqs, fileName, suffix=''):
    '''
    Compiles 'seqs' to a hardware description and saves it to 'fileName'. Other inputs:
        suffix : string to append to end of fileName (e.g. with fileNames = 'test' and suffix = 'foo' might save to test-APSfoo.h5)
    '''

    # Add the digitizer trigger to measurements
    PatternUtils.add_digitizer_trigger(seqs)

    # Add gating/blanking pulses
    PatternUtils.add_gate_pulses(seqs)

    # Add the slave trigger
    PatternUtils.add_slave_trigger(seqs, channelLib['slaveTrig'])

    # find channel set at top level to account for individual sequence channel variability
    channels = set([])
    for seq in seqs:
        channels |= find_unique_channels(seq)

    # Compile all the pulses/pulseblocks to sequences of pulses and control flow
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

    # construct channel delay map
    delays = channel_delay_map(physWires)

    # apply delays
    for chan, wire in physWires.items():
        PatternUtils.delay(wire, delays[chan])

    # generate wf library (base shapes)
    wfs = generate_waveforms(physWires)

    # replace Pulse objects with Waveforms
    physWires = pulses_to_waveforms(physWires)

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
        instrumentLib[awgName].write_sequence_file(data, fullFileName)

        fileList.append(fullFileName)

    # Return the filenames we wrote
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
    # inject function definitions prior to sequences
    funcDefs = collect_specializations(seqs)
    if funcDefs:
        # inject GOTO to jump over definitions
        funcDefs.insert(0, ControlFlow.Goto(BlockLabel.label(seqs[0])))
        seqs.insert(0, funcDefs)

    # use seqs[0] as prototype in case we were not given a set of channels
    wires = compile_sequence(seqs[0], channels)
    if not channels:
        channels = set(wires.keys())
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
                    wires[chan][-1] = copy(wires[chan][-1])
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
        if not hasattr(step, 'channel'):
            continue
        if isinstance(step.channel, Channels.Channel):
            channels |= set([step.channel])
        else:
            channels |= set(step.channel)
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
            self.key = pulse.hashshape()
            self.amp = pulse.amp
            self.length = pulse.shapeParams['length']
            self.phase = pulse.phase
            self.frameChange = pulse.frameChange
            self.isTimeAmp = pulse.isTimeAmp
            if hasattr(pulse.channel, 'frequency'):
                self.frequency = pulse.channel.frequency

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        labelPart = "{0}: ".format(self.label) if self.label else ""
        if self.isTimeAmp:
            TA = 'HIGH' if self.amp != 0 else 'LOW'
            return labelPart + "Waveform-TA(" + TA + ", " + str(self.length) + ")"
        else:
            return labelPart + "Waveform(" + str(self.key)[:6] + ", " + str(self.length) + ")"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(frozenset(self.__dict__.iteritems()))

    @property
    def isZero(self):
        return self.amp == 0

def schedule(channel, pulse, blockLength, alignment):
    '''
    Converts a Pulse or a CompositePulses into an aligned sequence of Pulses
    by injecting TAPulses before and/or after such that the resulting sequence
    duration is `blockLength`.
        alignment = "left", "right", or "center"
    '''
    # make everything look like a sequence
    if isinstance(pulse, PulseSequencer.CompositePulse):
        pulses = pulse.pulses
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
