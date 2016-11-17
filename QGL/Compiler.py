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
import logging
import numpy as np
import os
import operator
from warnings import warn
from copy import copy
from functools import reduce
from importlib import import_module
import json

from . import config
from . import PatternUtils
from .PatternUtils import flatten, has_gate
from . import Channels
from . import ChannelLibrary
from . import PulseShapes
from .PulsePrimitives import Id
from .PulseSequencer import Pulse, PulseBlock, CompositePulse
from . import ControlFlow
from . import BlockLabel

logger = logging.getLogger(__name__)

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
        entry_iterators = [iter(wires[ch][ct]) for ch in channels]
        while True:
            try:
                entries = [next(e) for e in entry_iterators]
            except StopIteration:
                break
            # control flow on any channel should pass thru
            if any(isinstance(e, (ControlFlow.ControlInstruction,
                                  BlockLabel.BlockLabel)) for e in entries):
                # for the moment require uniform control flow so that we
                # can pull from the first channel
                assert all(e == entries[0]
                           for e in entries), "Non-uniform control flow"
                segment.append(entries[0])
                continue
            # at this point we have at least one waveform instruction
            entries, block_length = pull_uniform_entries(entries, entry_iterators)

            # look for the simplest case of at most one non-identity
            if len(entries) == 1:
                segment.extend(entries[0])
                continue

            # create a new Pulse object for the merged pulse

            # If there is a non-zero SSB frequency copy it to the new entry
            nonZeroSSBChan = np.nonzero(
                [e.amp * e.frequency for e in entries])[0]
            assert len(nonZeroSSBChan) <= 1, \
                "Unable to handle merging more than one non-zero entry with non-zero frequency."
            if nonZeroSSBChan:
                frequency = entries[nonZeroSSBChan[0]].frequency
            else:
                frequency = 0.0

            if all([e.shapeParams['shapeFun'] == PulseShapes.constant for e in entries]):
                phasor = np.sum([e.amp * np.exp(1j * e.phase) for e in entries])
                amp = np.abs(phasor)
                phase = np.angle(phasor)
                shapeFun = PulseShapes.constant
            else:
                amp = 1.0
                phase = 0.0
                pulsesHash = tuple([(e.hashshape(), e.amp, e.phase) for e in entries])
                if pulsesHash not in shapeFunLib:
                    # create closure to sum waveforms
                    def sum_shapes(entries=entries, **kwargs):
                        return reduce(operator.add,
                            [e.amp * np.exp(1j * e.phase) * e.shape for e in entries])

                    shapeFunLib[pulsesHash] = sum_shapes
                shapeFun = shapeFunLib[pulsesHash]

            shapeParams = {"shapeFun": shapeFun, "length": block_length}

            label = "*".join([e.label for e in entries])
            segment.append(Pulse(label, entries[0].channel, shapeParams, amp, phase))

    return mergedWire


def pull_uniform_entries(entries, entry_iterators):
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
    #keep track of how many entry iterators are used up
    iterDone = [ False ] * numChan
    ct = 0
    lengths = [e.length for e in entries]
    entries_stack = [[e] for e in entries]

    while True:
        #If we've used up all the entries on all the channels we're done
        if all(iterDone):
            raise StopIteration("Unable to find a uniform set of entries")

        if all(x == lengths[0] for x in lengths):
            break

        #Otherwise try to concatenate on entries to match lengths
        while lengths[ct] < max(lengths):
            # concatenate with following entry to make up the length difference
            try:
                next_entry = next(entry_iterators[ct])
            except StopIteration:
                iterDone[ct] = True

            entries_stack[ct].append(next_entry)
            lengths[ct] += next_entry.length

        ct = (ct + 1) % numChan


    #if we have all zeros or a single non zero we return that as a list of entries
    allzero_chan = [all([e.isZero for e in stack]) for stack in entries_stack ]
    if all(allzero_chan):
        entries = [ entries_stack[0] ]
    elif np.sum(allzero_chan) == len(entries) - 1:
        entries = [ entries_stack[allzero_chan.index(False)] ]
    else:
        entries = []
        for stack in entries_stack:
            entries.append( reduce(concatenate_entries, stack) )

    return entries, max(lengths)


def concatenate_entries(entry1, entry2):
    # TA waveforms with the same amplitude can be merged with a just length update
    # otherwise, need to concatenate the pulse shapes
    shapeParams = copy(entry1.shapeParams)
    shapeParams["length"] = entry1.length + entry2.length
    label = entry1.label
    amp = entry1.amp
    phase = entry1.phase
    frameChange = entry1.frameChange + entry2.frameChange
    if not (entry1.isTimeAmp and entry2.isTimeAmp and entry1.amp == entry2.amp
            and entry1.phase == (entry1.frameChange + entry2.phase)):
        # otherwise, need to build a closure to stack them
        def stack_shapes(entry1=entry1, entry2=entry2, **kwargs):
            return np.hstack((
                entry1.amp * np.exp(1j * entry1.phase) * entry1.shape,
                entry2.amp * np.exp(1j * (entry1.frameChange + entry2.phase)) *
                entry2.shape))

        shapeParams['shapeFun'] = stack_shapes
        label = entry1.label + '+' + entry2.label
        amp = 1.0
        phase = 0.0

    return Pulse(label, entry1.channel, shapeParams, amp, phase, frameChange)


def generate_waveforms(physicalWires):
    wfs = {ch: {} for ch in physicalWires.keys()}
    for ch, wire in physicalWires.items():
        for pulse in flatten(wire):
            if not isinstance(pulse, Pulse):
                continue
            if pulse.hashshape() not in wfs[ch]:
                if pulse.isTimeAmp:
                    wfs[ch][pulse.hashshape()] = np.ones(1, dtype=np.complex)
                else:
                    wfs[ch][pulse.hashshape()] = pulse.shape
    return wfs


def pulses_to_waveforms(physicalWires):
    logger.debug("Converting pulses_to_waveforms:")
    wireOuts = {ch: [] for ch in physicalWires.keys()}
    for ch, seqs in physicalWires.items():
        logger.debug('')
        logger.debug("Channel '%s':", ch)
        for seq in seqs:
            wireOuts[ch].append([])
            for pulse in seq:
                if not isinstance(pulse, Pulse):
                    wireOuts[ch][-1].append(pulse)
                    logger.debug(" %s", pulse)
                else:
                    wf = Waveform(pulse)
                    wireOuts[ch][-1].append(wf)
                    logger.debug(" %s", wf)
    return wireOuts


def channel_delay_map(physicalWires):
    chanDelays = {chan: chan.delay for chan in physicalWires.keys()}
    return PatternUtils.normalize_delays(chanDelays)


def setup_awg_channels(physicalChannels):
    translators = {}
    for chan in physicalChannels:
        translators[chan.AWG] = import_module('QGL.drivers.' + chan.translator)

    data = {awg: translator.get_empty_channel_set()
            for awg, translator in translators.items()}
    for name, awg in data.items():
        for chan in awg.keys():
            awg[chan] = {
                'linkList': [],
                'wfLib': {},
                'correctionT': np.identity(2)
            }
        data[name]['translator'] = translators[name]
        data[name]['seqFileExt'] = translators[name].get_seq_file_extension()
    return data


def bundle_wires(physWires, wfs):
    awgData = setup_awg_channels(physWires.keys())
    for chan in physWires.keys():
        _, awgChan = chan.label.split('-')
        awgChan = 'ch' + awgChan
        awgData[chan.AWG][awgChan]['linkList'] = physWires[chan]
        awgData[chan.AWG][awgChan]['wfLib'] = wfs[chan]
        if hasattr(chan, 'correctionT'):
            awgData[chan.AWG][awgChan]['correctionT'] = chan.correctionT
    return awgData


def collect_specializations(seqs):
    '''
    Collects function definitions for all targets of Call instructions
    '''
    targets = [x.target for x in flatten(seqs)
               if isinstance(x, ControlFlow.Call)]
    funcs = []
    done = []
    #Manually keep track of done (instead of `set`) to keep calling order
    for target in targets:
        if target not in done:
            funcs.append(ControlFlow.qfunction_specialization(target))
            done.append(target)
    return funcs


def compile_to_hardware(seqs,
                        fileName,
                        suffix='',
                        axis_descriptor=None,
                        qgl2=False,
                        addQGL2SlaveTrigger=False,
                        edgesToCompile = None,
                        qubitToCompile = None):
    '''
    Compiles 'seqs' to a hardware description and saves it to 'fileName'.
    Other inputs:
        suffix (optional): string to append to end of fileName, e.g. with
            fileNames = 'test' and suffix = 'foo' might save to test-APSfoo.h5
        axis_descriptor (optional): a list of dictionaries describing the effective
            axes of the measurements that the sequence will yield. For instance,
            if `seqs` generates a Ramsey experiment, axis_descriptor would describe
            the time delays between pulses.
        qgl2 (optional): When True, run QGL2 specific debug code and only
            selectively add the slaveTrigger; when False (default), always add
            the slave trigger
        addQGL2SlaveTrigger (optional): When qgl2=True only add the slave trigger
            when this is also True
        edgesToCompile (optional): When not None, only compile edges (and their
            gates) on this list; else compile all edges and their gates (default)
        qubitToCompile (optional): When not None, only compile the given Qubit
            (and its gate), else compile all Qubits and their gates (default)
    '''
    logger.debug("Compiling %d sequence(s)", len(seqs))

    # save input code to file
    if not qgl2:
        save_code(seqs, fileName + suffix)

    # all sequences should start with a WAIT for synchronization
    for seq in seqs:
        if not isinstance(seq[0], ControlFlow.Wait):
            logger.debug("Adding a WAIT - first sequence element was %s", seq[0])
            seq.insert(0, ControlFlow.Wait())

    # Add the digitizer trigger to measurements
    logger.debug("Adding digitizer trigger")
    PatternUtils.add_digitizer_trigger(seqs)

    # Add gating/blanking pulses
    logger.debug("Adding blanking pulses")
    PatternUtils.add_gate_pulses(seqs)

    if not qgl2 or addQGL2SlaveTrigger:
        # Add the slave trigger
        logger.debug("Adding slave trigger")
        PatternUtils.add_slave_trigger(seqs,
                                       ChannelLibrary.channelLib['slaveTrig'])
    else:
        logger.debug("Not adding slave trigger")

    # find channel set at top level to account for individual sequence channel variability
    channels = set()
    for seq in seqs:
        channels |= find_unique_channels(seq)

    # Compile all the pulses/pulseblocks to sequences of pulses and control flow
    wireSeqs = compile_sequences(seqs, channels, edgesToCompile=edgesToCompile, qubitToCompile=qubitToCompile)

    if not validate_linklist_channels(wireSeqs.keys()):
        print("Compile to hardware failed")
        return

    logger.debug('')
    logger.debug("Now after gating constraints:")
    # apply gating constraints
    for chan, seq in wireSeqs.items():
        if isinstance(chan, Channels.LogicalMarkerChannel):
            wireSeqs[chan] = PatternUtils.apply_gating_constraints(
                chan.physChan, seq)
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug('')
            logger.debug("Channel '%s':", chan)
            for elem in seq:
                if isinstance(elem, list):
                    for e2 in elem:
                        logger.debug(" %s", e2)
                    logger.debug('')
                else:
                    logger.debug(" %s", elem)

    # save number of measurements for meta info
    num_measurements = count_measurements(wireSeqs)

    # map logical to physical channels
    physWires = map_logical_to_physical(wireSeqs)

    # construct channel delay map
    delays = channel_delay_map(physWires)

    # apply delays
    for chan, wire in physWires.items():
        PatternUtils.delay(wire, delays[chan])
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug('')
            logger.debug("Delayed wire for chan '%s':", chan)
            for elem in wire:
                if isinstance(elem, list):
                    for e2 in elem:
                        logger.debug(" %s", e2)
                    logger.debug('')
                else:
                    logger.debug(" %s", elem)

    # generate wf library (base shapes)
    wfs = generate_waveforms(physWires)

    # replace Pulse objects with Waveforms
    physWires = pulses_to_waveforms(physWires)

    # bundle wires on instruments
    awgData = bundle_wires(physWires, wfs)

    # convert to hardware formats
    files = {}
    for awgName, data in awgData.items():
        # create the target folder if it does not exist
        targetFolder = os.path.split(os.path.normpath(os.path.join(
            config.AWGDir, fileName)))[0]
        if not os.path.exists(targetFolder):
            os.mkdir(targetFolder)
        fullFileName = os.path.normpath(os.path.join(
            config.AWGDir, fileName + '-' + awgName + suffix + data[
                'seqFileExt']))
        data['translator'].write_sequence_file(data, fullFileName)

        files[awgName] = fullFileName

    # create meta output
    if not axis_descriptor:
        axis_descriptor = [{
            'name': 'segment',
            'unit': None,
            'points': list(range(1, 1 + num_measurements)),
            'partition': 1
        }]
    meta = {
        'instruments': files,
        'num_sequences': len(seqs),
        'num_measurements': num_measurements,
        'axis_descriptor': axis_descriptor,
    }
    metafilepath = os.path.join(config.AWGDir, fileName + '-meta.json')
    with open(metafilepath, 'w') as FID:
        json.dump(meta, FID, indent=2, sort_keys=True)

    # Return the filenames we wrote
    return list(files.values())


def compile_sequences(seqs, channels=set(), edgesToCompile=None, qubitToCompile=None):
    '''
    Main function to convert sequences to miniLL's and waveform libraries.
        edgesToCompile : When not None, only compile edges (and their gates) on
            this list; else compile all edges and their gates (default)
        qubitToCompile : When not None, only compile the given Qubit (and its
            gate), else compile all Qubits and their gates (default)
    '''

    # turn into a loop, by appending GOTO(0) at end of last sequence
    if not isinstance(seqs[-1][-1], ControlFlow.Goto):
        seqs[-1].append(ControlFlow.Goto(BlockLabel.label(seqs[0])))
        logger.debug("Appending a GOTO at end to loop")

    # append function specialization to sequences
    subroutines = collect_specializations(seqs)
    seqs += subroutines

    #expand the channel definitions for anything defined in subroutines
    for func in subroutines:
        channels |= find_unique_channels(subroutines)

    # use seqs[0] as prototype in case we were not given a set of channels
    wires = compile_sequence(seqs[0], channels, edgesToCompile, qubitToCompile)
    if not channels:
        channels = set(wires.keys())
    wireSeqs = {chan: [seq] for chan, seq in wires.items()}
    for seq in seqs[1:]:
        wires = compile_sequence(seq, channels, edgesToCompile, qubitToCompile)
        for chan in wireSeqs.keys():
            wireSeqs[chan].append(wires[chan])
    #Print a message so for the experiment we know how many sequences there are
    print('Compiled {} sequences.'.format(len(seqs) - len(subroutines)))

    # Debugging:
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug('')
        logger.debug("Returning from compile_sequences():")
        for chan in wireSeqs:
            logger.debug('')
            logger.debug("Channel '%s':", chan)
            seq = wireSeqs[chan]
            for elem in seq:
                if isinstance(elem, list):
                    for e2 in elem:
                        logger.debug(" %s", e2)
                else:
                    logger.debug(" %s", elem)

    return wireSeqs

def channels_to_compile(channels, edgesToCompile=None, qubitToCompile=None):
    '''
    Filter input set of channels to exclude those not on given lists
        edgesToCompile : When not None, only compile edges (and their gates) on
            this list; else compile all edges and their gates (default)
        qubitToCompile : When not None, only compile the given Qubit (and its
            gate), else compile all Qubits and their gates (default)
    '''
    # Edited set of channels to return
    newchans = set()

    # Skip the gates for any Edges/Qubits that we are skipping
    gatesToSkip = set()
    for chan in channels:
        if isinstance(chan, Channels.Edge):
            if edgesToCompile is not None and chan not in edgesToCompile:
                if has_gate(chan):
                    gatesToSkip.add(chan.gateChan)
                    logger.debug("%s to be excluded so excluding its gate %s", chan, chan.gateChan)
                else:
                    logger.debug("%s didn't have a gate", chan)
        elif isinstance(chan, Channels.Qubit):
            if qubitToCompile is not None and chan != qubitToCompile:
                if has_gate(chan):
                    gatesToSkip.add(chan.gateChan)
                    logger.debug("%s to be excluded so excluding its gate %s", chan, chan.gateChan)
                else:
                    logger.debug("%s didn't have a gate", chan)

    # Now for each channel skip the relevant channels, building
    # a new set of channels
    for chan in channels:
        if chan in gatesToSkip:
            logger.debug("Compilation skipping gate %s", chan)
            continue
        if isinstance(chan, Channels.Edge):
            if edgesToCompile is None or chan in edgesToCompile:
                # add chan
                logger.debug("Compilation including %s", chan)
                newchans.add(chan)
                continue
            else:
                # skip this edge
                logger.debug("Compilation skipping %s", chan)
                continue
        if isinstance(chan, Channels.Qubit):
            if qubitToCompile is None or chan == qubitToCompile:
                # add
                newchans.add(chan)
                logger.debug("Compilation including %s", chan)
                continue
            else:
                # skip
                logger.debug("Compilation skipping %s", chan)
                continue
        newchans.add(chan)
        logger.debug("Compilation including %s", chan)
    return newchans

def compile_sequence(seq, channels=None, edgesToCompile=None, qubitToCompile=None):
    '''
    Takes a list of control flow and pulses, and returns aligned blocks
    separated into individual abstract channels (wires).
        edgesToCompile : When not None, only compile edges (and their gates) on
            this list; else compile all edges and their gates (default)
        qubitToCompile : When not None, only compile the given Qubit (and its
            gate), else compile all Qubits and their gates (default)
    '''
    logger.debug('')
    logger.debug("In compile_sequence:")
    #Find the set of logical channels used here and initialize them
    if not channels:
        logger.debug("<Had no channels>")
        channels = find_unique_channels(seq)

    # Filter list of channels to only include designated edges and Qubit (and their gates)
    if edgesToCompile is not None or qubitToCompile is not None:
        channels = channels_to_compile(channels, edgesToCompile, qubitToCompile)

    wires = {chan: [] for chan in channels}

    # Debugging: what does the sequence look like?
    if logger.isEnabledFor(logging.DEBUG):
        for elem in seq:
            logger.debug(" %s", elem)
        logger.debug('')
        logger.debug("Channels:")
        for chan in channels:
            logger.debug(" %s", chan)

    logger.debug('')
    logger.debug("Sequence before normalizing:")
    for block in normalize(flatten(seq), channels):
        logger.debug(" %s", block)
        # labels and control flow instructions broadcast to all channels
        if isinstance(block,
                      (BlockLabel.BlockLabel, ControlFlow.ControlInstruction)):
            for chan in channels:
                wires[chan] += [copy(block)]
            continue
        # propagate frame change from nodes to edges
        for chan in channels:
            if block.pulses[chan].frameChange == 0:
                continue
            if chan in ChannelLibrary.channelLib.connectivityG.nodes():
                logger.debug("Doing propagate_node_frame_to_edges()")
                wires = propagate_node_frame_to_edges(
                    wires, chan, block.pulses[chan].frameChange)
        # drop length 0 blocks but push nonzero frame changes onto previous entries
        if block.length == 0:
            for chan in channels:
                if block.pulses[chan].frameChange == 0:
                    continue
                if len(wires[chan]) == 0:
                    warn("Dropping initial frame change")
                    continue
                logger.debug("Modifying pulse on %s: %s", chan, wires[chan][-1])
                # search for last non-TA entry
                for ct in range(1,len(wires[chan])):
                    if hasattr(wires[chan][-ct], 'isTimeAmp') and not wires[chan][-ct].isTimeAmp:
                        updated_frameChange = wires[chan][-ct].frameChange + block.pulses[chan].frameChange
                        wires[chan][-ct] = wires[chan][-ct]._replace(frameChange=updated_frameChange)
                        break
            continue

        # schedule the block
        for chan in channels:
            # add aligned Pulses (if the block contains a composite pulse, may get back multiple pulses)
            wires[chan] += schedule(chan, block.pulses[chan], block.length,
                                    block.alignment)
    if logger.isEnabledFor(logging.DEBUG):
        for chan in wires:
            logger.debug('')
            logger.debug("compile_sequence() return for channel '%s':", chan)
            for elem in wires[chan]:
                logger.debug(" %s", elem)
    return wires


def propagate_node_frame_to_edges(wires, chan, frameChange):
    '''
    Propagate frame change in node to relevant edges (for CR gates)
    '''
    for predecessor in ChannelLibrary.channelLib.connectivityG.predecessors(
            chan):
        edge = ChannelLibrary.channelLib.connectivityG.edge[predecessor][chan][
            'channel']
        if edge in wires:
            # search for last non-TA entry
            for ct in range(1,len(wires[edge])):
                if hasattr(wires[edge][-ct], 'isTimeAmp') and not wires[edge][-ct].isTimeAmp:
                    updated_frameChange = wires[edge][-ct].frameChange + frameChange
                    wires[edge][-ct] = wires[edge][-ct]._replace(frameChange=updated_frameChange)
                    break
    return wires


def find_unique_channels(seq):
    channels = set()
    for step in flatten(seq):
        if not hasattr(step, 'channel'):
            continue
        if isinstance(step.channel, Channels.Channel):
            channels.add(step.channel)
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
    for block in filter(lambda x: isinstance(x, PulseBlock), seq):
        blocklen = block.length
        emptyChannels = channels - set(block.pulses.keys())
        for ch in emptyChannels:
            block.pulses[ch] = Id(ch, blocklen)
    return seq

class Waveform(object):
    """
    Simplified channel independent version of a Pulse with a key into waveform library.
    """


    def __init__(self, pulse=None):
        if pulse is None:
            self.label = ""
            self.key = None
            self.amp = 0
            self.length = 0
            self.phase = 0
            self.frameChange = 0
            self.isTimeAmp = False
            self.frequency = 0
        else:
            self.label = pulse.label
            self.key = pulse.hashshape()
            self.amp = pulse.amp
            self.length = pulse.shapeParams['length']
            self.phase = pulse.phase
            self.frameChange = pulse.frameChange
            self.isTimeAmp = pulse.isTimeAmp
            self.frequency = pulse.frequency

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self.isTimeAmp:
            TA = 'HIGH' if self.amp != 0 else 'LOW'
            return "Waveform-TA(" + TA + ", " + str(self.length) + ")"
        else:
            return "Waveform(" + self.label + ", " + str(
                self.key)[:6] + ", " + str(self.length) + ")"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(frozenset(self.__dict__.items()))

    @property
    def isZero(self):
        return self.amp == 0


def schedule(channel, pulse, blockLength, alignment):
    '''
    Converts a Pulse or a CompositePulses into an aligned sequence of Pulses
    by injecting Ids before and/or after such that the resulting sequence
    duration is `blockLength`.
        alignment = "left", "right", or "center"
    '''
    # make everything look like a sequence
    if isinstance(pulse, CompositePulse):
        pulses = pulse.pulses
    else:
        pulses = [pulse]

    padLength = blockLength - pulse.length
    if padLength == 0:
        # logger.debug("   schedule on chan '%s' made no change", channel)
        pass
    else:
        logger.debug("   schedule on chan '%s' adding Id len %d align %s",
                     channel, padLength, alignment)
    if padLength == 0:
        # no padding element required
        return pulses
    elif alignment == "left":
        return pulses + [Id(channel, padLength)]
    elif alignment == "right":
        return [Id(channel, padLength)] + pulses
    else:  # center
        return [Id(channel, padLength / 2)] + pulses + [Id(channel, padLength /
                                                           2)]


def validate_linklist_channels(linklistChannels):
    errors = []
    channels = ChannelLibrary.channelLib.channelDict
    for channel in linklistChannels:
        if channel.label not in channels.keys(
        ) and channel.label not in errors:
            print("{0} not found in channel library".format(repr(channel)))
            errors.append(channel.label)

    if errors != []:
        return False
    return True


def set_log_level(loggerName='QGL.Compiler',
                  levelDesired=logging.DEBUG,
                  formatStr='%(message)s'):
    '''Set the python logging level for a logger.
    Format messages with just the log message by default.
    Sets QGL.Compiler messages to DEBUG by default.
    '''
    import logging
    import sys
    # Do basicConfig to be safe, but ask for console to be STDOUT
    # so it doesn't look like an error in a jupyter notebook
    logging.basicConfig(level=levelDesired,
                        format=formatStr,
                        stream=sys.stdout)

    # Enable our logger at specified level
    logger = logging.getLogger(loggerName)
    logger.setLevel(levelDesired)

    # Find the first stream handler to STDOUT and set its log level & format
    handlers = logger.handlers
    if not handlers:
        handlers = logging.getLogger().handlers
    updatedH = False
    for handler in handlers:
        if isinstance(handler, logging.StreamHandler) and \
           handler.stream == sys.stdout:
            handler.setLevel(levelDesired)
            handler.setFormatter(logging.Formatter(formatStr))
            return

    # No existing handler, so add a local one
    handler = logging.StreamHandler(
        stream=sys.stdout)  # Note we request STDOUT
    handler.setLevel(levelDesired)
    handler.setFormatter(logging.Formatter(formatStr))
    logger.addHandler(handler)
    # Without this when running in notebook, console gets log stmts at default format
    logger.propagate = 0


def save_code(seqs, filename):
    from IPython.lib.pretty import pretty
    import io  #needed for writing unicode to file in Python 2.7
    # create the target folder if it does not exist
    targetFolder = os.path.split(os.path.normpath(os.path.join(config.AWGDir,
                                                               filename)))[0]
    if not os.path.exists(targetFolder):
        os.mkdir(targetFolder)
    fullname = os.path.normpath(os.path.join(config.AWGDir, filename +
                                             '-code.py'))
    with io.open(fullname, "w", encoding="utf-8") as FID:
        FID.write(u'seqs =\n')
        FID.write(pretty(seqs))

def count_measurements(wireSeqs):
    # count number of measurements per sequence as the max over the the number
    # of measurements per wire

    # pick an arbitrary key to determine sequence length
    seq_len = len(wireSeqs[list(wireSeqs)[0]])
    seq_measurements = [0 for _ in range(seq_len)]
    for ct in range(seq_len):
        for wire, seqs in wireSeqs.items():
            seq_measurements[ct] = max(
                seq_measurements[ct],
                sum(PatternUtils.contains_measurement(e) for e in seqs[ct]))
    return sum(seq_measurements)
