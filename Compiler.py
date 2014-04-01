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
import hashlib
import os
from warnings import warn

import config
import PatternUtils
import Channels
from PulsePrimitives import Id
import PulseSequencer
from instruments.AWGs import get_empty_channel_set, APS, Tek5014

SEQUENCE_PADDING = 480 #2800 #480
from APSPattern import write_APS_file
from TekPattern import write_Tek_file

# global parameter libraries
channelLib = {}
instrumentLib = {}

def hash_pulse(shape):
    return hashlib.sha1(shape.tostring()).hexdigest()
    
TAZKey = hash_pulse(np.zeros(1, dtype=np.complex))
markerHighKey = hash_pulse(np.ones(1, dtype=np.bool))

markerWFLib = {TAZKey:np.zeros(1, dtype=np.bool), markerHighKey:np.ones(1, dtype=np.bool) }

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
        # dig in an grab the associated gate channel, too
        if not isinstance(chan, Channels.LogicalMarkerChannel):
            awgs.add(chan.physChan.gateChan.AWG)
    return {awg.label:get_empty_channel_set(awg) for awg in awgs}


def map_logical_to_physical(linkLists, wfLib):
    physicalChannels = {chan: channelLib[get_channel_label(chan)].physChan.label for chan in linkLists.keys()}
    awgData = setup_awg_channels(linkLists.keys())

    for chan in linkLists.keys():
        awgName, awgChan = physicalChannels[chan].split('-')
        awgData[awgName]['ch'+awgChan] = {'linkList': linkLists[chan], 'wfLib': wfLib[chan]}

    return awgData


def compile_to_hardware(seqs, fileName=None, suffix='', alignMode="right", nbrRepeats=1):
    #Add the digitizer trigger to each sequence
    #TODO: Make this more sophisticated.
    PatternUtils.add_digitizer_trigger(seqs, channelLib['digitizerTrig'])
    
    # normalize sequences
    channels = set([])
    for seq in seqs:
        channels |= find_unique_channels(seq)
    seqs = [normalize(seq, channels) for seq in seqs]

    #Compile all the pulses/pulseblocks to linklists and waveform libraries
    linkLists, wfLib = compile_sequences(seqs)

    # align channels
    # this horrible line finds the longest miniLL across all channels
    longestLL = max([sum([entry.totLength for entry in miniLL]) for LL in linkLists.values() for miniLL in LL])
    
    for chan, LL in linkLists.items():
        PatternUtils.align(LL, alignMode, longestLL+SEQUENCE_PADDING)

    #Add the slave trigger
    #TODO: only add to slave devices
    linkLists[channelLib['slaveTrig']], wfLib[channelLib['slaveTrig']] = PatternUtils.slave_trigger(len(seqs))

    # map logical to physical channels
    awgData = map_logical_to_physical(linkLists, wfLib)

    # for each physical channel need to:
    # 1) delay
    # 2) apply SSB if necessary
    # 3) mixer correct
    for awgName, awg in awgData.items():
        for chanName, chanData in awg.items():
            if chanData:
                # construct IQkey using existing convention
                IQkey = awgName + '-' + chanName[2:]
                chanObj = channelLib[IQkey]
                #We handle marker and quadrature channels differently
                #For now that is all we handle
                if isinstance(chanObj, Channels.PhysicalQuadratureChannel):
                    #Apply mixer corrections and channel delay 
                    PatternUtils.delay(chanData['linkList'], chanObj.delay + chanObj.AWG.delay, chanObj.samplingRate)

                    #At this point we finally have the timing of all the pulses so we can apply SSB
                    if hasattr(chanObj, 'SSBFreq') and abs(chanObj.SSBFreq) > 0:
                        PatternUtils.apply_SSB(chanData['linkList'], chanData['wfLib'], chanObj.SSBFreq, chanObj.samplingRate)

                    PatternUtils.correctMixer(chanData['wfLib'], chanObj.correctionT)

                    # add gate pulses on the marker channel
                    # note that the marker may be on an entirely different AWG
                    markerAwgName, markerKey = chanObj.gateChan.label.split('-')
                    markerKey = 'ch' + markerKey
                    markerAwg = awgData[markerAwgName]
                    genObj = chanObj.generator
                    #TODO: check if this actually catches overwriting markers
                    if markerAwg[markerKey]:
                        warn('Reuse of marker gating channel: {0}'.format(markerKey))
                    markerAwg[markerKey] = {'linkList':None, 'wfLib':markerWFLib}
                    markerAwg[markerKey]['linkList'] = PatternUtils.create_gate_seqs(
                        chanData['linkList'], genObj.gateBuffer, genObj.gateMinWidth, chanObj.samplingRate)
                    markerDelay = genObj.gateDelay + chanObj.gateChan.delay + (chanObj.AWG.delay - chanObj.gateChan.AWG.delay)
                    PatternUtils.delay(markerAwg[markerKey]['linkList'], markerDelay, chanObj.gateChan.samplingRate )


                elif isinstance(chanObj, Channels.PhysicalMarkerChannel):
                    PatternUtils.delay(chanData['linkList'], chanObj.delay+chanObj.AWG.delay, chanObj.samplingRate)

                else:
                    raise NameError('Unable to handle channel type.')

                #Remove unused waveforms
                compress_wfLib(chanData['linkList'], chanData['wfLib'])

    #Loop back through to fill empty channels and write to file
    fileList = []
    for awgName, awg in awgData.items():
        #If all the channels are empty then do not bother writing the file
        if not all([chan is None for chan in awg.values()]):
            for chan in awg.keys():
                if not awg[chan]:
                    #"Seems hackish but check for marker
                    if chan[-2] == 'm':
                        awg[chan] = {'linkList': [[create_padding_LL(SEQUENCE_PADDING//2), create_padding_LL(SEQUENCE_PADDING//2)]],
                                 'wfLib': markerWFLib}
                    else:
                        awg[chan] = {'linkList': [[create_padding_LL(SEQUENCE_PADDING//2), create_padding_LL(SEQUENCE_PADDING//2)]],
                                 'wfLib': {TAZKey:np.zeros(1, dtype=np.complex)}}

            # convert to hardware formats
            # create the target folder if it does not exist
            targetFolder = os.path.split(os.path.normpath(os.path.join(config.AWGDir, fileName)))[0]
            if not os.path.exists(targetFolder):
                os.mkdir(targetFolder)
            fullFileName = os.path.normpath(os.path.join(config.AWGDir, fileName + '-' + awgName + suffix + instrumentLib[awgName].seqFileExt))
            if isinstance(instrumentLib[awgName], APS):
                write_APS_file(awg, fullFileName, nbrRepeats)
            elif isinstance(instrumentLib[awgName], Tek5014):
                assert nbrRepeats == 1, 'nbrRepeats > 1 not implemented for the Tek'
                write_Tek_file(awg, fullFileName, fileName)
            else:
                raise NameError('Unknown AWG type')
            fileList.append(fullFileName)

    #Return the filenames we wrote
    return fileList


def compile_sequences(seqs):
    '''
    Main function to convert sequences to miniLL's and waveform libraries.
    '''
    if isinstance(seqs[0], list):
        # nested sequences
        wfLib = {}
        # use seqs[0] as prototype for finding channels (assume every miniLL operates on the same set of channels)
        miniLL, wfLib = compile_sequence(seqs[0], wfLib)
        linkLists = {chan: [LL] for chan, LL in miniLL.items()}
        for seq in seqs[1:]:
            miniLL, wfLib = compile_sequence(seq, wfLib)
            for chan in linkLists.keys():
                linkLists[chan].append(miniLL[chan])
    else:
        miniLL, wfLib = compile_sequence(seqs)
        linkLists = {chan: [LL] for chan, LL in miniLL.items()}

    #Print a message so for the experiment we know how many sequences there are
    print('Compiled {} sequences.'.format(len(seqs)))
    return linkLists, wfLib

def compile_sequence(seq, wfLib={} ):
    '''
    Converts a single sequence into a miniLL and waveform library.
    Returns a single-entry list of a miniLL and the updated wfLib
    '''

    #Find the set of logical channels used here and initialize them
    channels = find_unique_channels(seq)

    logicalLLs = {}        
    for chan in channels:
        logicalLLs[chan] = []
        if chan not in wfLib:
            if isinstance(chan, Channels.LogicalMarkerChannel):
                wfLib[chan] = markerWFLib
            else:
                wfLib[chan] = {TAZKey:  np.zeros(1, dtype=np.complex)}
    carriedPhase = {ch: 0 for ch in channels}
    for block in seq:
        #Align the block 
        blockLength = block.maxPts
        # drop length 0 blocks but push frame change onto next non-zero entry
        if blockLength == 0:
            carriedPhase = {ch: carriedPhase[ch]+block.pulses[ch].frameChange for ch in channels}
            continue
        for chan in channels:
            # add aligned LL entry(ies) (if the block contains a composite pulse, may get back multiple waveforms and LL entries)
            wfs, LLentries = align(block.pulses[chan], blockLength, block.alignment)
            for wf in wfs:
                if isinstance(chan, Channels.LogicalMarkerChannel):
                    wf = wf.astype(np.bool)
                if hash_pulse(wf) not in wfLib:
                    wfLib[chan][hash_pulse(wf)] = wf
            # Frame changes are then propagated through
            logicalLLs[chan] += propagate_frame(LLentries, carriedPhase[chan])
        carriedPhase = {ch: 0 for ch in channels}

    # loop through again to find phases, frame changes, and SSB modulation for quadrature channels
    for chan, miniLL in logicalLLs.items():
        if isinstance(chan, Channels.Qubit) or isinstance(chan, Channels.Measurement):
            curFrame = 0
            for entry in miniLL:
                # frame update
                shape = np.copy(wfLib[chan][entry.key])

                # See if we can turn into a TA pair
                # fragile: if you buffer a square pulse it will not be constant valued
                if np.all(shape == shape[0]):
                    entry.isTimeAmp = True
                    shape = shape[:1]
                    # convert near zeros to TAZKey
                    if abs(shape[0]) < 1e-6:
                        entry.key = TAZKey

                #Rotate for phase and frame change (don't rotate zeros...)
                if entry.key != TAZKey:
                    shape *= np.exp(1j*(entry.phase+curFrame))
                    shapeHash = hash_pulse(shape)
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
    usedKeys = set([TAZKey, markerHighKey])
    for miniLL in seqs:
        for entry in miniLL:
            usedKeys.add(entry.key)

    unusedKeys = set(wfLib.keys()) - usedKeys
    for key in unusedKeys:
        del wfLib[key]

def find_unique_channels(seq):
    channels = set([])
    for step in seq:
        if isinstance(step, PulseSequencer.Pulse):
            if isinstance(step.qubits, Channels.Channel):
                channels |= set([step.qubits])
            else:
                channels |= set(step.qubits)
        else:
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
    for block in seq:
        emptyChannels = channels - set(block.pulses.keys())
        for ch in emptyChannels:
            block.pulses[ch] = Id(ch, length=0)
    return seq



class LLElement(object):
    '''
    IQ LL elements for quadrature mod channels.
    '''
    def __init__(self, pulse=None):
        self.repeat = 1

        if pulse is None:
            self.key = None
            self.length = 0
            self.phase = 0
            self.frameChange = 0
            self.isTimeAmp = False
        else:
            self.key = hash_pulse(pulse.shape)
            self.length = pulse.length
            self.phase = pulse.phase
            self.frameChange = pulse.frameChange
            self.isTimeAmp = pulse.isTimeAmp
    
    @property
    def hasMarker(self):
        return self.markerDelay1 or self.markerDelay2

    @property
    def isZero(self):
        return self.key == TAZKey

    @property
    def totLength(self):
        return self.length*self.repeat

def create_padding_LL(length=SEQUENCE_PADDING, high=False):
    tmpLL = LLElement()
    tmpLL.isTimeAmp = True
    tmpLL.key = markerHighKey if high else TAZKey
    tmpLL.length = length
    return tmpLL

def align(pulse, blockLength, alignment, cutoff=12):
    # check for composite pulses
    if hasattr(pulse, 'pulses'):
        entries = [LLElement(p) for p in pulse.pulses]
        shapes = [p.shape for p in pulse.pulses]
    else:
        entries = [LLElement(pulse)]
        shapes = [pulse.shape]
    padLength = blockLength - pulse.length
    if padLength == 0:
        # no padding element required
        return shapes, entries
    if (padLength < cutoff) and (alignment == "left" or alignment == "right"):
        # pad the first/last shape on one side
        if alignment == "left":
            shapes[-1] = np.hstack((shapes[-1], np.zeros(padLength)))
            entries[-1].key = hash_pulse(shapes[-1])
        else: #right alignment
            shapes[0] = np.hstack((np.zeros(padLength), shapes[0]))
            entries[0].key = hash_pulse(shapes[0])
    elif (padLength < 2*cutoff and alignment == "center"):
        # pad the both sides of the shape(s)
        if len(shapes) == 1:
            shapes[0] = np.hstack(( np.zeros(np.floor(padLength/2)), shapes[0], np.zeros(np.ceil(padLength/2)) ))
            entries[0].key = hash_pulse(shapes[0])
        else:
            shapes[0] = np.hstack(( np.zeros(np.floor(padLength/2)), shapes[0]))
            shapes[-1] = np.hstack(( np.zeroes(np.ceil(padLength/2)), shapes[-1]))
            entries[0].key = hash_pulse(shapes[0])
            entries[-1].key = hash_pulse(shapes[-1])
    elif padLength == blockLength:
        #Here we have a zero-length sequence which just needs to be expanded
        entries[0].key = TAZKey
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
