'''
Module for writing hdf5 APS files from LL's and patterns
'''

import h5py
import numpy as np
from warnings import warn
from itertools import chain

#Some constants
ADDRESS_UNIT = 4 #everything is done in units of 4 timesteps
MIN_ENTRY_LENGTH = 16
MIN_LL_ENTRY_COUNT = 2 #minimum length of mini link list
MAX_WAVEFORM_PTS = 2**15 #maximum size of waveform memory
MAX_WAVEFORM_VALUE = 2**13-1 #maximum waveform value i.e. 14bit DAC
MAX_LL_ENTRIES = 8192 #maximum number of LL entries in a bank
MAX_REPEAT_COUNT = 2^10-1;

#APS bit masks
START_MINILL_BIT = 15;
END_MINILL_BIT = 14;
WAIT_TRIG_BIT = 13;
TA_PAIR_BIT = 12;

TAZKey = hash(tuple(np.zeros(1, dtype=np.complex)))

def preprocess_APS(miniLL, wfLib):
    '''
    Helper function to deal with LL elements less than minimum LL entry count  
    by trying to concatenate them into neighbouring entries
    '''
    newMiniLL = []
    entryct = 0
    while entryct < len(miniLL):
        curEntry = miniLL[entryct]
        if curEntry.length > MIN_ENTRY_LENGTH:
            newMiniLL.append(curEntry)
            entryct += 1
        else:
            nextEntry = miniLL[entryct+1]
            previousEntry = miniLL[entryct-1] if entryct > 0 else None
            
            #For short TA pairs we see if we can add them to the next waveform
            if curEntry.key == TAZKey and not nextEntry.key == TAZKey:
                print("Got here")
                #Concatenate the waveforms                
                paddedWF = np.hstack((wfLib[curEntry.key]*np.ones(curEntry.length), wfLib[nextEntry.key]))
                #Hash the result to generate a new unique key and add
                newKey = hash(tuple(paddedWF))
                wfLib[newKey] = paddedWF
                nextEntry.key = newKey
                nextEntry.length = wfLib[newKey].size
                newMiniLL.append(nextEntry)
                entryct += 2
            
            #For short pulses we see if we can steal some padding from the previous or next entry
            elif previousEntry.key == TAZKey and previousEntry.length > 2*MIN_ENTRY_LENGTH:
                padLength = MIN_ENTRY_LENGTH - curEntry.length
                newMiniLL[-1].length -= padLength
                #Concatenate the waveforms                
                paddedWF = np.hstack((np.zeros(padLength, dtype=np.complex), wfLib[curEntry.key]))
                #Hash the result to generate a new unique key and add
                newKey = hash(tuple(paddedWF))
                wfLib[newKey] = paddedWF
                curEntry.key = newKey
                curEntry.length = wfLib[newKey].size
                newMiniLL.append(curEntry)
                entryct += 1

            elif nextEntry.key == TAZKey and nextEntry.length > 2*MIN_ENTRY_LENGTH:
                padLength = MIN_ENTRY_LENGTH - curEntry.length
                newMiniLL[-1].length -= padLength
                #Concatenate the waveforms                
                paddedWF = np.hstack((wfLib[curEntry.key], np.zeros(padLength, dtype=np.complex)))
                #Hash the result to generate a new unique key and add
                newKey = hash(tuple(paddedWF))
                wfLib[newKey] = paddedWF
                curEntry.key = newKey
                curEntry.length = wfLib[newKey].size
                newMiniLL.append(curEntry)
                entryct += 1

            else:
                warn("Unable to handle too short LL element, dropping.")
                entryct += 1

    #Update the miniLL 
    return newMiniLL

def create_wf_vector(wfLib):
    '''
    Helper function to create the wf vector and offsets into it.
    '''
    wfVec = np.zeros(MAX_WAVEFORM_PTS, dtype=np.int16)
    offsets = {}
    idx = 0
    for key, wf in wfLib.items():
        #Clip the wf
        wf[wf>1] = 1.0
        wf[wf<-1] = -1.0
        #TA pairs need to be repeated ADDRESS_UNIT times
        if wf.size == 1:
            wf = wf.repeat(ADDRESS_UNIT)
        #Ensure the wf is an integer number of ADDRESS_UNIT's 
        trim = wf.size%ADDRESS_UNIT
        if trim:
            wf = wf[:-trim]
        wfVec[idx:idx+wf.size] = np.uint16(MAX_WAVEFORM_VALUE*wf)
        offsets[key] = idx
        idx += wf.size
                    
    #Trim the waveform 
    wfVec = waveformLib[0:idx] 

    return wfVec, offsets

def create_LL_data(LLs, offsets):
    '''
    Helper function to create LL data vectors from a list of miniLL's and an offset dictionary
    keyed on the wf keys.
    '''

    #Preallocate the bank data and do some checking for miniLL lengths
    seqLengths = np.array([len(miniLL) for miniLL in LLs])
    assert np.all(seqLengths >= 3), 'Oops! mini LL''s needs to have at least three elements.'
    assert np.all(seqLengths < MAX_BANK_SIZE), 'Oops! mini LL''s cannot have length greater than {0}, you have {1} entries'.format(MAX_BANK_SIZE, len(miniLL))
    numEntries = sum(seqLengths)
    LLData = {label: np.zeros(numEntries, dtype=np.uint16) for label in ['addr','count', 'trigger1', 'trigger2', 'repeat']}

    #Loop over all entries
    TAPairEntries = []
    for ct, entry in enumerate(chain.from_iterable(LLs)):
        LLData['addr'][ct] = offsets[entry.key]
        LLData['count'][ct] = entry.length//ADDRESS_UNIT-1
        LLData['trigger1'][ct], LLData['trigger2'][ct] = calc_trigger(entry)
        LLData['repeat'][ct] = entry.repeat-1
        if entry.isTAPair:
            TAPairEntries.append(ct)


    #Add in the miniLL start/stop and TA pair flags on the upper bits of the repeat entries
    startPts = np.hstack((0, seqLengths))
    endPts = np.hstack((seqLengths-1, numEntries-1)
    LLData['repeat'][startPts] += 2**START_MINILL_BIT + 2**WAIT_TRIG_BIT
    LLData['repeat'][endPts] += 2**END_MINILL_BIT
    LLData['repeat'][TAPairEntries] += 2**TA_PAIR_BIT
    
    return LLData

def write_APS_file(LLs12, wfLib12, LLs34, wfLib34, chanData34, fileName, miniLLRepeat):
    '''
    Main function to pack channel LLs into an APS h5 file.

    '''
    #Open the HDF5 file
    with h5py.File(fileName, 'w') as FID:  
    
        #List of which channels we have data for
        #TODO: actually handle incomplete channel data
        channelDataFor = [1,2,3,4]
        FID['/'].attrs['Version'] = 2.0
        FID['/'].attrs['channelDataFor'] = np.uint16(channelDataFor)
        FID['/'].attrs['miniLLRepeat'] = np.uint16(miniLLRepeat)
   
        #Create the waveform vectors
        wfVecs = []
        offsets = []
        wfVec, offsets = create_wf_vector({key:wf.real for key,wf in wfLib12})
        wfVecs.append(wfVec)
        offsets.append(offsets)
        wfVec, offsets = create_wf_vector({key:wf.imag for key,wf in wfLib12})
        wfVecs.append(wfVec)
        offsets.append(offsets)
        wfVec, offsets = create_wf_vector({key:wf.real for key,wf in wfLib34})
        wfVecs.append(wfVec)
        offsets.append(offsets)
        wfVec, offsets = create_wf_vector({key:wf.imag for key,wf in wfLib34})
        wfVecs.append(wfVec)
        offsets.append(offsets)

        LLData = [LLs12, LLs34]
        #Create the groups and datasets
        for chanct in range(4):
            chanStr = '/chan_{0}'.format(chanct+1)
            chanGroup = FID.create_group(chanStr)
            chanGroup.attrs['isIQMode'] = np.uint8(1)
            #Write the waveformLib to file
            FID.create_dataset('{0}/waveformLib'.format(chanStr), data=waveformLib)

            #For A channels (1 & 3) we write link list data
            if np.mod(chanct,2) == 0:
                chanGroup.attrs['isLinkListData'] = np.uint8(1)
                groupStr = chanStr+'/linkListData'
                LLGroup = FID.create_group(groupStr)
                LLDataVecs = create_LL_data(LLData[chanct//2])
                LLGroup.attrs['length'] = np.uint16(LLDataVecs['addr'].size)
                for key,dataVec in LLDataVecs.items():
                    FID.create_dataset(groupStr+'/' + key, data=dataVec)


def read_APS_file(fileName):
    '''
    Helper function to read back in data from a H5 file and reconstruct the sequence
    '''
    AWGData = {}
    #APS bit masks
    START_MINILL_MASK = 2**START_MINILL_BIT;
    TA_PAIR_MASK = 2**TA_PAIR_BIT;
    REPEAT_MASK = 2**10-1
            
    
    with h5py.File(fileName, 'r') as FID:
        chanct = 0
        for chanct, chanStr in enumerate(chanStrs2):
            #If we're in IQ mode then the Q channel gets its linkListData from the I channel
            if FID[chanStr].attrs['isIQMode'][0]:
                tmpChan = 2*(chanct//2)
                curLLData = FID[chanStrs2[tmpChan]]['linkListData']
            else:
                curLLData = FID[chanStr]['linkListData']
            #Pull out the LL data
            tmpAddr = curLLData['addr'].value
            tmpCount = curLLData['count'].value
            tmpRepeat = curLLData['repeat'].value
            tmpTrigger1 = curLLData['trigger1'].value
            tmpTrigger2 = curLLData['trigger2'].value
            numEntries = curLLData.attrs['length'][0]
   
            #Pull out and scale the waveform data
            wfLib =(1.0/MAX_WAVEFORM_VALUE)*FID[chanStr]['waveformLib'].value.flatten()
            
            #Initialize the lists of sequences
            AWGData[chanStrs[chanct]] = []
            AWGData[mrkStrs[chanct]] = []

            #Loop over LL entries
            for entryct in range(numEntries):
                #If we are starting a new entry push back an empty array
                if START_MINILL_MASK & tmpRepeat[entryct]:
                    AWGData[chanStrs[chanct]].append(np.array([], dtype=np.float64))
                    AWGData[mrkStrs[chanct]].append(np.array([], dtype=np.bool))
                #If it is a TA pair or regular pulse
                curRepeat = (tmpRepeat[entryct] & REPEAT_MASK)+1
                if TA_PAIR_MASK & curLLData['repeat'][entryct][0]:
                    AWGData[chanStrs[chanct]][-1] = np.hstack((AWGData[chanStrs[chanct]][-1], 
                                                    np.tile(wfLib[tmpAddr[entryct]*ADDRESS_UNIT:tmpAddr[entryct]*ADDRESS_UNIT+4], curRepeat*(tmpCount[entryct]+1))))
                else:
                    AWGData[chanStrs[chanct]][-1] = np.hstack((AWGData[chanStrs[chanct]][-1], 
                                                    np.tile(wfLib[tmpAddr[entryct]*ADDRESS_UNIT:tmpAddr[entryct]*ADDRESS_UNIT+4*(tmpCount[entryct]+1)], curRepeat)))
                #Add the trigger pulse
                tmpPulse = np.zeros(ADDRESS_UNIT*curRepeat*(tmpCount[entryct]+1), dtype=np.bool)
                if chanct//2 == 0:
                    if tmpTrigger1[entryct] > 0:
                        tmpPulse[4*tmpTrigger1[entryct]] = True
                else:
                    if tmpTrigger2[entryct] > 0:
                        tmpPulse[4*tmpTrigger2[entryct]] = True
                AWGData[mrkStrs[chanct]][-1] = np.hstack((AWGData[mrkStrs[chanct]][-1], tmpPulse)) 
                
    return AWGData


if __name__ == '__main__':

    pass
