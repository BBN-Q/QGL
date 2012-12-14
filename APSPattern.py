'''
Module for writing hdf5 APS files from LL's and patterns
'''

import h5py
import numpy as np
from warnings import warn

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

    #Update the miniLL in place
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

    for miniLLct, miniLL in enumerate(AWGData[chanStr]['LLs']):
        LLlength = len(miniLL)
        #The minimum miniLL length is two 
        assert LLlength >= 3, 'Oops! mini LL''s needs to have at least three elements.'
        assert LLlength < MAX_BANK_SIZE, 'Oops! mini LL''s cannot have length greater than {0}, you have {1} entries'.format(MAX_BANK_SIZE, len(miniLL))
        #If we need to allocate a new bank
        if entryct + len(miniLL) > MAX_BANK_SIZE:
            #Fix the final entry as we no longer have to indicate the next enty is the start of a miniLL
            tmpBank['offset'][entryct-1] -= ELL_FIRST_ENTRY
            #Write the current bank to file
            write_bank_to_file(FID, tmpBank, '/{0}/linkListData/bank{1}'.format(chanStrs2[chanct], bankct), entryct)
            #Allocate a new bank
            tmpBank = create_empty_bank()
            bankct += 1
            #Reset the entry count
            entryct = 0
        
        #Otherwise enter each LL entry into the bank arrays
        for ct, LLentry in enumerate(miniLL):
            tmpBank['offset'][entryct] = calc_offset(LLentry, offsets, entryct==0 or (ct==LLlength-1 and miniLLct<numMiniLLs-1) , ct==LLlength-2)
            tmpBank['count'][entryct] = LLentry.length//ADDRESS_UNIT-1
            tmpBank['trigger'][entryct] = calc_trigger(LLentry)
            tmpBank['repeat'][entryct] = LLentry.repeat-1
            entryct += 1


def write_APS_file(LL12, wfLib12, LL34, wfLib34, chanData34, fileName):
    '''
    Main function to pack channel LLs into an APS h5 file.

    '''
    #Open the HDF5 file
    with h5py.File(fileName, 'w') as FID:  
    
        #List of which channels we have data for
        # channelDataFor = []

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

        #Create the groups and datasets
        for chanct in range(4):
            chanStr = '/ch{0}'.format(chanct+1)
            chanGroup = FID.create_group(chanStr)

            #Write the waveformLib to file
            FID.create_dataset('{0}/waveformLib'.format(chanStr), data=waveformLib)

            #For A channels (1 & 3) we write link list data
            if np.mod(chanct,2) == 0:
                chanGroup.attrs['isLinkListData'] = np.int16(1)



                
        # #Loop over the channels
        # for chanct, chanStr in enumerate(chanStrs):
        #     if chanStr in AWGData:
        #         channelDataFor.append(chanct+1)
    
        #         if AWGData[chanStr]['WFLibrary']:
                    
                
        #         #Create the LL data group
        #         LLGroup = FID.create_group('/'+chanStrs2[chanct] + '/linkListData')
                
        #         #Create the necessary number of banks as we step through the mini LL
        #         entryct = 0
        #         tmpBank = create_empty_bank()
        #         bankct = 1
        #         numMiniLLs = len(AWGData[chanStr]['LLs'])
                        
        #         #Write the final bank
        #         write_bank_to_file(FID, tmpBank, '/{0}/linkListData/bank{1}'.format(chanStrs2[chanct], bankct), entryct)
                
        #         LLGroup.attrs['numBanks'] = np.uint16(bankct)
        #         LLGroup.attrs['repeatCount'] = np.uint16(0)
                    
                        
            
        # FID['/'].attrs['Version'] = 2.0
        # FID['/'].attrs['channelDataFor'] = np.int16(channelDataFor)
    
    
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
