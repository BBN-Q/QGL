'''
functions for compiling lists of pulses/pulseBlocks down to the hardware level.
'''

import numpy as np
import json
import AWG
import PatternUtils
import config
from Channels import ChannelDict
import Channels
from PulsePrimitives import Id

from warnings import warn

from APSPattern import  write_APS_file

SEQUENCE_PADDING = 244

def get_channel_name(chanKey):
    ''' Takes in a channel key and returns a channel name '''
    if type(chanKey) != tuple:
        return chanKey.name
    else:
        return "".join([chan.name for chan in chanKey])

def setup_awg_channels(logicalChannels):
    awgs = set([])
    for chan in logicalChannels:
        awgs.add(ChannelDict[get_channel_name(chan)].physicalChannel.AWG)
    return {awg.name: getattr(AWG, awg.model)().channels() for awg in awgs}

def map_logical_to_physical(linkLists, wfLib):
    physicalChannels = {chan: ChannelDict[get_channel_name(chan)].physicalChannel.name for chan in linkLists.keys()}
    awgData = setup_awg_channels(linkLists.keys())
    
    for chan in linkLists.keys():
        awgName, awgChan = physicalChannels[chan].split('-')
        awgData[awgName]['ch'+awgChan] = {'linkList': linkLists[chan], 'wfLib': wfLib[chan]}
    
    return awgData

def compile_to_hardware(seqs, fileName=None, suffix='', alignMode="right"):
    #Add the digitizer trigger to each sequence
    #TODO: Make this more sophisticated.
    PatternUtils.add_digitizer_trigger(seqs, ChannelDict['digitizerTrig'])

    # normalize sequences
    seqs = [normalize(seq) for seq in seqs]

    #Compile all the pulses/pulseblocks to linklists and waveform libraries
    linkLists, wfLib = compile_sequences(seqs)

    # align channels
    # this horrible line finds the longest miniLL across all channels
    longestLL = max([sum([entry.length*entry.repeat for entry in miniLL]) for LL in linkLists.values() for miniLL in LL])
    for chan, LL in linkLists.items():
        PatternUtils.align(LL, alignMode, longestLL+SEQUENCE_PADDING)

    #Add the slave trigger
    #TODO: only add to slave devices
    linkLists[ChannelDict['slaveTrig']], wfLib[ChannelDict['slaveTrig']] = PatternUtils.slave_trigger(len(seqs))

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
                chanObj = ChannelDict[IQkey]
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
                    genObj = chanObj.generator
                    markerKey = 'ch' + genObj.gateChannel.name.split('-')[1]
                    #TODO: check if this actually catches overwriting markers
                    if awg[markerKey]:
                        warn('Reuse of marker gating channel: {0}'.format(markerKey))
                    awg[markerKey] = {'linkList':None, 'wfLib':None}
                    awg[markerKey]['linkList'] = PatternUtils.create_gate_seqs(
                        chanData['linkList'], genObj.gateBuffer, genObj.gateMinWidth, chanObj.samplingRate)
                    PatternUtils.delay(awg[markerKey]['linkList'], genObj.gateDelay+chanObj.delay, genObj.gateChannel.samplingRate )

                elif isinstance(chanObj, Channels.PhysicalMarkerChannel):
                    PatternUtils.delay(chanData['linkList'], chanObj.delay+chanObj.AWG.delay, chanObj.samplingRate)

                else:
                    raise NameError('Unable to handle channel type.')

    #Loop back through to fill empty channels and write to file
    fileList = []
    for awgName, awg in awgData.items():
        #If all the channels are empty then do not bother writing the file
        if not all([chan is None for chan in awg.values()]):
            for chan in awg.keys():
                if not awg[chan]:
                    awg[chan] = {'linkList': [[create_padding_LL(SEQUENCE_PADDING//2), create_padding_LL(SEQUENCE_PADDING//2)]],
                                 'wfLib': {TAZKey:np.zeros(1, dtype=np.complex)}}

            # convert to hardware formats
            if ChannelDict[awgName].model == 'BBNAPS':
                tmpFileName = config.AWGDir + fileName + '-' + awgName + suffix + '.h5'
                write_APS_file(awg, tmpFileName )
                fileList.append(tmpFileName)

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

    #Compress the waveform library
    for chan in linkLists.keys():
        compress_wfLib(linkLists[chan], wfLib[chan])

    #Print a message so for the experiment we know how many sequences there are
    print('Compiled {} sequences.'.format(len(seqs)))
    return linkLists, wfLib

def compile_sequence(seq, wfLib = {} ):
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
            wfLib[chan] = {TAZKey:  np.zeros(1, dtype=np.complex)}
    carriedPhase = {ch: 0 for ch in channels}
    for block in seq:
        #Align the block 
        blockLength = block.maxPts
        # drop length 0 blocks but push frame change onto previous entry
        if blockLength == 0:
            # Push frame change forward through the next pulse
            carriedPhase = {ch: carriedPhase[ch]+block.pulses[ch].frameChange for ch in channels}
            # continue to drop the entry
            continue
        for chan in channels:
            # add aligned LL entry
            wf, LLentry = align(block.pulses[chan], blockLength, block.alignment)
            if hash_pulse(wf) not in wfLib:
                wfLib[chan][hash_pulse(wf)] = wf
            LLentry[0].phase -= carriedPhase[chan]
            LLentry[0].frameChange += carriedPhase[chan]
            logicalLLs[chan] += LLentry
        carriedPhase = {ch: 0 for ch in channels}

    # loop through again to find phases, frame changes, and SSB modulation
    for chan, miniLL in logicalLLs.items():
        curFrame = 0
        for entry in miniLL:
            # frame update
            shape = np.copy(wfLib[chan][entry.key])

            # See if we can turn into a TA pair
            # fragile: if you buffer a square pulse it will not be constant valued
            if np.all(shape == shape[0]):
                entry.isTimeAmp = True
                shape = shape[:1]

            #Rotate for phase and frame change 
            shape *= np.exp(1j*(entry.phase+curFrame))
            shapeHash = hash_pulse(shape)
            if shapeHash not in wfLib[chan]:
                wfLib[chan][shapeHash] = shape
            entry.key = shapeHash
            curFrame += entry.frameChange 

    return logicalLLs, wfLib

def compress_wfLib(seqs, wfLib):
    '''
    Helper function to remove unused waveforms from the library.
    '''
    usedKeys = set()
    for miniLL in seqs:
        for entry in miniLL:
            usedKeys.add(entry.key)

    unusedKeys = set(wfLib.keys()) - usedKeys
    for key in unusedKeys:
        del wfLib[key]

def find_unique_channels(seq):
    channels = set([])
    for step in seq:
        channels |= set(step.pulses.keys())
    return channels

def normalize(seq):
    '''
    For mixed lists of Pulses and PulseBlocks, converts to list of PulseBlocks
    with uniform channels on each PulseBlock. We inject Id's where necessary.
    '''
    # promote to PulseBlocks
    seq = [p.promote() for p in seq]

    channels = set(find_unique_channels(seq))

    # inject Id's for PulseBlocks not containing every channel
    for block in seq:
        emptyChannels = channels - set(block.pulses.keys())
        for ch in emptyChannels:
            block.pulses[ch] = Id(ch, length=0)
    return seq

def hash_pulse(shape):
    # if we need more speed, this version is about 10x faster in my tests on arrays of length 2000
    #return hashlib.sha1(shape.view(np.uint8)).hexdigest()
    return hash(tuple(shape))

TAZKey = hash_pulse(np.zeros(1, dtype=np.complex))
markerHighKey = hash_pulse([True])

class LLElement(object):
    '''
    IQ LL elements for quadrature mod channels.
    '''
    def __init__(self, pulse=None):
        self.repeat = 1
        self.isTimeAmp = False

        if pulse is None:
            self.key = None
            self.length = 0
            self.phase = 0
            self.frameChange = 0
        else:
            self.key = hash_pulse(pulse.shape)
            self.length = len(pulse.shape)
            self.phase = pulse.phase
            self.frameChange = pulse.frameChange
    
    @property
    def hasMarker(self):
        return self.markerDelay1 or self.markerDelay2

    @property
    def isZero(self):
        return self.key == TAZKey

def create_padding_LL(length, high=False):
    tmpLL = LLElement()
    tmpLL.isTimeAmp = True
    tmpLL.key = markerHighKey if high else TAZKey
    tmpLL.length = length
    return tmpLL

def align(pulse, blockLength, alignment, cutoff=12):
    entry = LLElement(pulse)
    entry.length = pulse.shape.size
    entry.key = hash_pulse(pulse.shape)
    entry.phase = pulse.phase
    entry.frameChange = pulse.frameChange
    padLength = blockLength - pulse.shape.size
    shape = pulse.shape
    if padLength == 0:
        # can do everything with a single LLentry
        return shape, [entry]
    if (padLength < cutoff) and (alignment == "left" or alignment == "right"):
        # pad the shape on one side
        if alignment == "left":
            shape = np.hstack((shape, np.zeros(padLength)))
        else: #right alignment
            shape = np.hstack((np.zeros(padLength), shape))
        entry.key = hash_pulse(shape)
        return shape, [entry]
    elif (padLength < 2*cutoff and alignment == "center"):
        # pad the shape on each side
        shape = np.hstack(( np.zeros(np.floor(padLength/2)), shape, np.zeros(np.ceil(padLength/2)) ))
        entry.key = hash_pulse(shape)
        return shape, [entry]
    elif padLength == blockLength:
        #Here we have a zero-length sequence which just needs to be expanded
        entry.key = TAZKey
        entry.length = blockLength
        return np.zeros(1, dtype=np.complex), [entry]
    else:
        #split the entry into the shape and one or more TAZ
        if alignment == "left":
            padEntry = create_padding_LL(padLength)
            return shape, [entry, padEntry]
        elif alignment == "right":
            padEntry = create_padding_LL(padLength)
            return shape, [padEntry, entry]
        else:
            padEntry1 = create_padding_LL(np.floor(padLength/2))
            padEntry2 = create_padding_LL(np.ceil(padLength/2))
            return shape, [padEntry1, entry, padEntry2]
