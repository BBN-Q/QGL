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
from . import ChannelLibraries
from . import PulseShapes
from .PulsePrimitives import Id, clear_pulse_cache
from .PulseSequencer import Pulse, PulseBlock, CompositePulse
from . import ControlFlow
from . import BlockLabel
from . import TdmInstructions # only for APS2-TDM
import gc

logger = logging.getLogger(__name__)

def map_logical_to_physical(wires):
    """construct a mapping of physical channels to lists of logical channels
    (there will be more than one logical channel if multiple logical
    channels share a physical channel)"""
    physicalChannels = {}
    for logicalChan in wires.keys():
        phys_chan = logicalChan.phys_chan
        if not phys_chan:
            raise ValueError("LogicalChannel {} does not have a PhysicalChannel".format(logicalChan))
        if phys_chan not in physicalChannels:
            physicalChannels[phys_chan] = [logicalChan]
        else:
            physicalChannels[phys_chan].append(logicalChan)

    # loop through the physical channels
    physicalWires = {}
    for phys_chan, logicalChan in physicalChannels.items():
        if len(logicalChan) > 1:
            physicalWires[phys_chan] = merge_channels(wires, logicalChan)
        else:
            physicalWires[phys_chan] = wires[logicalChan[0]]

    return physicalWires


def merge_channels(wires, channels):
    chan = channels[0]
    mergedWire = [[] for _ in range(len(wires[chan]))]
    shape_funLib = {}
    for ct, segment in enumerate(mergedWire):
        entry_iterators = [iter(wires[ch][ct]) for ch in channels]
        while True:
            try:
                entries = [next(e) for e in entry_iterators]
            except StopIteration:
                break
            # control flow on any channel should pass thru
            if any(isinstance(e, (ControlFlow.ControlInstruction,
                                  BlockLabel.BlockLabel, TdmInstructions.WriteAddrInstruction,
                                  TdmInstructions.CustomInstruction,
                                  TdmInstructions.LoadCmpVramInstruction)) for e in entries):
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
            if nonZeroSSBChan.size > 0:
                frequency = entries[nonZeroSSBChan[0]].frequency
            else:
                frequency = 0.0

            if all([e.shapeParams['shape_fun'] == PulseShapes.constant for e in entries]):
                phasor = np.sum([e.amp * np.exp(1j * e.phase) for e in entries])
                amp = np.abs(phasor)
                phase = np.angle(phasor)
                shape_fun = PulseShapes.constant
            else:
                amp = 1.0
                phase = 0.0
                pulsesHash = tuple([(e.hashshape(), e.amp, e.phase) for e in entries])
                if pulsesHash not in shape_funLib:
                    # create closure to sum waveforms
                    def sum_shapes(entries=entries, **kwargs):
                        return reduce(operator.add,
                            [e.amp * np.exp(1j * e.phase) * e.shape for e in entries])

                    shape_funLib[pulsesHash] = sum_shapes
                shape_fun = shape_funLib[pulsesHash]

            shapeParams = {"shape_fun": shape_fun, "length": block_length}

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
    lengths = np.array([e.length for e in entries])
    entries_stack = [[e] for e in entries]

    while True:
        #If we've used up all the entries on all the channels we're done
        if all(iterDone):
            raise StopIteration("Unable to find a uniform set of entries")

        #if all(np.isclose(x, lengths[0], atol=1e-10) for x in lengths):
        if np.all(np.less(np.abs(lengths - lengths[0]), 1e-10)):
            break

        #Otherwise try to concatenate on entries to match lengths
        while not np.isclose(lengths[ct], max(lengths), atol=1e-10):
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

        shapeParams['shape_fun'] = stack_shapes
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
        translators[chan.instrument] = import_module('QGL.drivers.' + chan.translator)

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
        _, awgChan = chan.label.rsplit('-', 1)
        if awgChan[0] != 'm':
            awgChan = 'ch' + awgChan
        awgData[chan.instrument][awgChan]['linkList'] = physWires[chan]
        awgData[chan.instrument][awgChan]['wfLib'] = wfs[chan]
        if hasattr(chan, 'correctionT'):
            awgData[chan.instrument][awgChan]['correctionT'] = chan.correctionT
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
                        library_version=None,
                        suffix='',
                        axis_descriptor=None,
                        add_slave_trigger=True,
                        extra_meta=None,
                        tdm_seq=False):
                        meas_qs=None,
                        meas_decoupled_qs=None,
                        CR_chs=None,
                        CR_decoupled_chs=None):
    '''
    Compiles 'seqs' to a hardware description and saves it to 'fileName'.
    Other inputs:
        library_version (optional): string or ChannelLibrary instance to pack in the
            metafile. This will be the version of the library loaded during program
            execution. Default None uses the current working version.
        suffix (optional): string to append to end of fileName, e.g. with
            fileNames = 'test' and suffix = 'foo' might save to test-APSfoo.h5
        axis_descriptor (optional): a list of dictionaries describing the effective
            axes of the measurements that the sequence will yield. For instance,
            if `seqs` generates a Ramsey experiment, axis_descriptor would describe
            the time delays between pulses.
        add_slave_trigger (optional): add the slave trigger(s)
        tdm_seq (optional): compile for TDM
    '''
    ChannelLibraries.channelLib.update_channelDict()
    clear_pulse_cache()

    logger.debug("Compiling %d sequence(s)", len(seqs))

    # save input code to file
    save_code(seqs, fileName + suffix)

    # add decoupling pulses
    PatternUtils.decouple_seqs(seqs, meas_qs, meas_decoupled_qs, CR_chs, CR_decoupled_chs)
    save_code(seqs, fileName + suffix)

    # all sequences should start with a WAIT for synchronization
    for seq in seqs:
        if not isinstance(seq[0], ControlFlow.Wait):
            logger.debug("Adding a WAIT - first sequence element was %s", seq[0])
            seq.insert(0, ControlFlow.Wait())

    # Add the digitizer trigger to measurements
    logger.info("Adding digitizer trigger")
    PatternUtils.add_digitizer_trigger(seqs)

    # Add gating/blanking pulses
    logger.info("Adding blanking pulses")
    for seq in seqs:
        PatternUtils.add_parametric_pulses(seq)
        PatternUtils.add_gate_pulses(seq)

    if add_slave_trigger and 'slave_trig' in ChannelLibraries.channelLib:
        # Add the slave trigger
        logger.debug("Adding slave trigger")
        PatternUtils.add_slave_trigger(seqs,
                                       ChannelLibraries.channelLib['slave_trig'])
    else:
        logger.info("Not adding slave trigger")

    # find channel set at top level to account for individual sequence channel variability
    logger.info("Finding unique channels.")
    channels = set()
    for seq in seqs:
        channels |= find_unique_channels(seq)

    # Compile all the pulses/pulseblocks to sequences of pulses and control flow
    logger.info("Compiling sequences.")
    wireSeqs = compile_sequences(seqs, channels)

    if not validate_linklist_channels(wireSeqs.keys()):
        print("Compile to hardware failed")
        return

    logger.debug('')
    logger.debug("Now after gating constraints:")
    # apply gating constraints
    logger.info("Applying gating constraints")
    for chan, seq in wireSeqs.items():
        if isinstance(chan, Channels.LogicalMarkerChannel):
            wireSeqs[chan] = PatternUtils.apply_gating_constraints(
                chan.phys_chan, seq)
    debug_print(wireSeqs, 'Gated sequence')

    # save number of measurements for meta info
    logger.info("Counting measurements.")
    num_measurements = count_measurements(wireSeqs)
    wire_measurements = count_measurements_per_wire(wireSeqs)

    # map logical to physical channels, physWires is a list of
    # PhysicalQuadratureChannels and PhysicalMarkerChannels
    # for the APS, the naming convention is:
    # ASPName-12, or APSName-12m1
    logger.info("Mapping logical to physical channels.")
    physWires = map_logical_to_physical(wireSeqs)

    # Pave the way for composite instruments, not useful yet...
    files = {}
    label_to_inst   = {}
    label_to_chan   = {}
    old_wire_names  = {}
    old_wire_instrs = {}
    for wire, pulses in physWires.items():
        pattern_module = import_module('QGL.drivers.' + wire.translator)
        if pattern_module.SEQFILE_PER_CHANNEL:
            inst_name = wire.transmitter.label
            chan_name = wire.label
            has_non_id_pulses = any([len([p for p in ps if isinstance(p,Pulse) and p.label!="Id"]) > 0 for ps in pulses])
            label_to_inst[wire.label] = inst_name
            if has_non_id_pulses:
                label_to_chan[wire.label] = chan_name
            # Change the name/inst for uniqueness, but we must restore this later!
            old_wire_names[wire] = wire.label
            old_wire_instrs[wire] = wire.instrument
            wire.instrument = wire.label
            wire.label = chan_name
            #files[inst_name] = {}

    # construct channel delay map
    logger.info("Constructing delay map.")
    delays = channel_delay_map(physWires)

    # apply delays
    logger.info("Applying delays.")
    for chan, wire in physWires.items():
        PatternUtils.delay(wire, delays[chan])
    debug_print(physWires, 'Delayed wire')

    # generate wf library (base shapes)
    logger.info("Generating waveform library.")
    wfs = generate_waveforms(physWires)

    # replace Pulse objects with Waveforms
    logger.info("Replacing pulses with waveforms")
    physWires = pulses_to_waveforms(physWires)

    # bundle wires on instruments, or channels depending
    # on whether we have one sequence per channel
    logger.info("Bundling wires.")
    awgData = bundle_wires(physWires, wfs)
    del wireSeqs
    gc.collect()

    # convert to hardware formats
    # files = {}
    awg_metas = {}
    for awgName in list(awgData.keys()):
        data = awgData[awgName]
        # create the target folder if it does not exist
        targetFolder = os.path.split(os.path.normpath(os.path.join(
            config.AWGDir, fileName)))[0]
        if not os.path.exists(targetFolder):
            os.mkdir(targetFolder)
        fullFileName = os.path.normpath(os.path.join(
            config.AWGDir, fileName + '-' + awgName + suffix + data[
                'seqFileExt']))
        logger.info("Writing sequence file for: {}".format(awgName))
        new_meta = data['translator'].write_sequence_file(data, fullFileName)
        if new_meta:
            awg_metas[awgName] = new_meta
            ChannelLibraries.channelLib[awgName].extra_meta = new_meta

        # Allow for per channel and per AWG seq files
        if awgName in label_to_inst:
            if awgName in label_to_chan:
                files[label_to_chan[awgName]] = fullFileName
        else:
            files[awgName] = fullFileName

        del data
        del awgData[awgName]
        gc.collect()

    # generate TDM sequences FIXME: what's the best way to identify the need for a TDM seq.? Support for single TDM
    if tdm_seq and 'APS2Pattern' in [wire.translator for wire in physWires]:
            aps2tdm_module = import_module('QGL.drivers.APS2Pattern') # this is redundant with above
            tdm_instr = aps2tdm_module.tdm_instructions(seqs)
            files['TDM'] = os.path.normpath(os.path.join(
                config.AWGDir, fileName + '-' + 'TDM' + suffix + '.aps2'))
            aps2tdm_module.write_tdm_seq(tdm_instr, files['TDM'])

    if extra_meta:
        extra_meta.update(awg_metas)
    else:
        extra_meta = awg_metas
    # create meta output
    db_info = {
        'db_provider': ChannelLibraries.channelLib.db_provider,
        'db_resource_name': ChannelLibraries.channelLib.db_resource_name,
        'library_name': 'working',
        'library_id': ChannelLibraries.channelLib.channelDatabase.id
    }
    if not axis_descriptor:
        axis_descriptor = [{
            'name': 'segment',
            'unit': None,
            'points': list(range(1, 1 + num_measurements)),
            'partition': 1
        }]
    receiver_measurements = {}
    for wire, n in wire_measurements.items():
        if wire.receiver_chan and n>0:
            receiver_measurements[wire.receiver_chan.label] = n
    meta = {
        'database_info': db_info,
        'instruments': files,
        'num_sequences': len(seqs),
        'num_measurements': num_measurements,
        'axis_descriptor': axis_descriptor,
        'qubits': [c.label for c in channels if isinstance(c, Channels.Qubit)],
        'measurements': [c.label for c in channels if isinstance(c, Channels.Measurement)],
        'receivers': receiver_measurements,
        'edges': [c.label for c in channels if isinstance(c, Channels.Edge)]
    }
    if extra_meta:
        meta.update({'extra_meta': extra_meta})
    metafilepath = os.path.join(config.AWGDir, fileName + '-meta.json')
    with open(metafilepath, 'w') as FID:
        json.dump(meta, FID, indent=2, sort_keys=True)

    # Restore the wire info
    for wire in old_wire_names.keys():
        wire.label = old_wire_names[wire]
    for wire in old_wire_instrs.keys():
        wire.instrument = old_wire_instrs[wire]

    # Return the filenames we wrote
    return metafilepath


def compile_sequences(seqs, channels=set()):
    '''
    Main function to convert sequences to miniLL's and waveform libraries.
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
    wires = compile_sequence(seqs[0], channels)
    if not channels:
        channels = set(wires.keys())
    wireSeqs = {chan: [seq] for chan, seq in wires.items()}
    for seq in seqs[1:]:
        wires = compile_sequence(seq, channels)
        for chan in wireSeqs.keys():
            wireSeqs[chan].append(wires[chan])
    #Print a message so for the experiment we know how many sequences there are
    print('Compiled {} sequences.'.format(len(seqs) - len(subroutines)))

    # Debugging:
    debug_print(wireSeqs, 'Return from compile_sequences()')

    return wireSeqs

def compile_sequence(seq, channels=None):
    '''
    Takes a list of control flow and pulses, and returns aligned blocks
    separated into individual abstract channels (wires).
    '''
    logger.debug('')
    logger.debug("In compile_sequence:")
    #Find the set of logical channels used here and initialize them
    if not channels:
        logger.debug("<Had no channels>")
        channels = find_unique_channels(seq)

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
    # make everything look like a sequence of pulses
    def flatten_to_pulses(obj):
        if isinstance(obj, Pulse):
            yield obj
        else:
            for pulse in obj.pulses:
                yield from flatten_to_pulses(pulse)
    for block in normalize(flatten(seq), channels):
        logger.debug(" %s", block)
        # labels broadcast to all channels
        if isinstance(block, BlockLabel.BlockLabel):
            for chan in channels:
                wires[chan] += [copy(block)]
            continue
        # control flow broadcasts to all channels if channel attribute is None
        if (isinstance(block, ControlFlow.ControlInstruction) or
                isinstance(block, TdmInstructions.WriteAddrInstruction) or
                isinstance(block, TdmInstructions.CustomInstruction) or
                isinstance(block, TdmInstructions.LoadCmpVramInstruction)):
            # Need to deal with attaching measurements and edges to control
            # instruction. Until we have a proper solution for that, we will
            # always broadcast control instructions to all channels
            # block_channels = block.channels if block.channels else channels
            block_channels = channels
            for chan in block_channels:
                wires[chan] += [copy(block)]
            continue
        # propagate frame change from nodes to edges
        for chan in channels:
            if block.pulses[chan].frameChange == 0:
                continue
            if chan in ChannelLibraries.channelLib.connectivityG.nodes():
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
                for ct in range(1,len(wires[chan])+1):
                    if hasattr(wires[chan][-ct], 'isTimeAmp') and not wires[chan][-ct].isTimeAmp:
                        updated_frameChange = wires[chan][-ct].frameChange + block.pulses[chan].frameChange
                        wires[chan][-ct] = wires[chan][-ct]._replace(frameChange=updated_frameChange)
                        break
                    # if no non-TA entry, add frame change at first available opportunity, from the start moving forward
                    if ct == len(wires[chan]):
                        for ct2 in range(len(wires[chan])):
                            if hasattr(wires[chan][ct2], 'frameChange'):
                                if hasattr(wires[chan][ct2], 'isTimeAmp') and  wires[chan][ct2].isTimeAmp:
                                    updated_frameChange = wires[chan][ct2].frameChange + block.pulses[chan].frameChange
                                    wires[chan][ct2] = wires[chan][ct2]._replace(frameChange=updated_frameChange)
                                break # if the first pulse is not an Id left, it should have already been updated
            continue
        # add the pulses per channel
        for chan in channels:
            wires[chan] += list(flatten_to_pulses(block.pulses[chan]))

    debug_print(wires, 'compile_sequence() return')

    return wires


def propagate_node_frame_to_edges(wires, chan, frameChange):
    '''
    Propagate frame change in node to relevant edges (for CR gates)
    '''
    for predecessor in ChannelLibraries.channelLib.connectivityG.predecessors(
            chan):
        edge = ChannelLibraries.channelLib.connectivityG.edges[predecessor, chan]['channel']
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
        if not hasattr(step, 'channel') or step.channel is None:
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
    seq = [p.promote(PulseBlock) for p in seq]

    if not channels:
        channels = find_unique_channels(seq)

    # inject Id's for PulseBlocks not containing every channel
    for block in filter(lambda x: isinstance(x, PulseBlock), seq):
        blocklen = block.length
        emptyChannels = channels - set(block.pulses.keys())
        for ch in emptyChannels:
            block.pulses[ch] = Id(ch, blocklen)
    return seq

class MemoizedObject(type):
    """Metaclass for memoizing objects. Designed to save memory when processing
    very long sequences. Overrides __call__ to prevent creating a new instance.
    Idea from https://stackoverflow.com/questions/47785795/memoized-objects-still-have-their-init-invoked
    """
    def __init__(self, name, bases, namespace):
        super().__init__(name, bases, namespace)
        self.cache = {}
    def __call__(self, pulse=None):
        hash = pulse.hashshape()
        if pulse is not None and hash not in self.cache:
            self.cache[hash] = super().__call__(pulse=pulse)
        return self.cache[hash]

class Waveform(object):
    """
    Simplified channel independent version of a Pulse with a key into waveform library.
    """

    #Use slots to create attributes to save on memory.
    __slots__ = ["label", "key", "amp", "length", "phase", "frameChange",
                    "isTimeAmp", "frequency", "logicalChan", "maddr", "startTime"]

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
            self.logicalChan = ""
            self.maddr = (-1, 0)
        else:
            self.label = pulse.label
            self.key = pulse.hashshape()
            self.amp = pulse.amp
            self.length = pulse.shapeParams['length']
            self.phase = pulse.phase
            self.frameChange = pulse.frameChange
            self.isTimeAmp = pulse.isTimeAmp
            self.frequency = pulse.frequency
            self.logicalChan = pulse.channel
            self.maddr = pulse.maddr

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
            #No __dict__ property so we check all properties.
            return all((getattr(self, attr, None) == getattr(other, attr, None) for attr in self.__slots__))
        return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(frozenset((attr, getattr(self, attr, None)) for attr in self.__slots__))

    @property
    def isZero(self):
        return self.amp == 0

def validate_linklist_channels(linklistChannels):
    errors = []
    channels = ChannelLibraries.channelLib.channelDict
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

def debug_print(seqs, label):
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug('')
        for chan, seq in seqs.items():
            logger.debug("%s for chan '%s':", label, chan)
            for elem in seq:
                if isinstance(elem, list):
                    for e2 in elem:
                        logger.debug(" %s", e2)
                    logger.debug('')
                else:
                    logger.debug(" %s", elem)

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
        seq_measurements[ct] = \
            reduce(max, count_measurements_per_wire_idx(wireSeqs, ct).values())
    return sum(seq_measurements)

def count_measurements_per_wire(wireSeqs):
    # pick an arbitrary key to determine sequence length
    seq_len = len(wireSeqs[list(wireSeqs)[0]])
    meas_wires = list(filter(lambda x: isinstance(x, Channels.Measurement), wireSeqs))
    measurements = {wire: 0 for wire in meas_wires}
    for ct in range(seq_len):
        seq_measurements = count_measurements_per_wire_idx(wireSeqs, ct)
        for wire in meas_wires:
            measurements[wire] += seq_measurements[wire]
    return measurements

def count_measurements_per_wire_idx(wireSeqs, idx):
    measurements = {
        wire: sum(PatternUtils.contains_measurement(e) for e in seqs[idx])
              for wire, seqs in wireSeqs.items()
    }
    return measurements
