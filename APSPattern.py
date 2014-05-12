'''
Module for writing hdf5 APS files from LL's and patterns

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

import h5py
import os
import numpy as np
from warnings import warn
from itertools import chain, izip_longest
from copy import copy, deepcopy
import hashlib

#Some constants
ADDRESS_UNIT = 4 #everything is done in units of 4 timesteps
MIN_ENTRY_LENGTH = 12
MIN_LL_ENTRY_COUNT = 2 #minimum length of mini link list
MAX_WAVEFORM_PTS = 2**15 #maximum size of waveform memory
MAX_WAVEFORM_VALUE = 2**13-1 #maximum waveform value i.e. 14bit DAC
MAX_LL_ENTRIES = 8192 #maximum number of LL entries in a bank
MAX_REPEAT_COUNT = 2**10-1;
MAX_TRIGGER_COUNT = 2**16-1

#APS bit masks
START_MINILL_BIT = 15;
END_MINILL_BIT = 14;
WAIT_TRIG_BIT = 13;
TA_PAIR_BIT = 12;

def hash_pulse(shape):
    return hashlib.sha1(shape.tostring()).hexdigest()

TAZKey = hash_pulse(np.zeros(1, dtype=np.complex))

def preprocess_APS(miniLL, wfLib):
	'''
	Helper function to deal with LL elements less than minimum LL entry count  
	by trying to concatenate them into neighbouring entries
	'''
	newMiniLL = []
	entryct = 0
	while entryct < len(miniLL):
		curEntry = miniLL[entryct]
		if curEntry.length >= MIN_ENTRY_LENGTH:
			newMiniLL.append(curEntry)
			entryct += 1
			continue

		if entryct == len(miniLL) - 1:
			# we've run out of entries to append to. drop it?
			warn("Unable to handle too short LL element, dropping.")
			break
		nextEntry = miniLL[entryct+1]
		previousEntry = miniLL[entryct-1] if entryct > 0 else None
		
		#For short TA pairs we see if we can add them to the next waveform
		if curEntry.key == TAZKey and not nextEntry.key == TAZKey:
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
			if curEntry.isTimeAmp:
				paddedWF = np.hstack((np.zeros(padLength, dtype=np.complex), wfLib[curEntry.key]*np.ones(curEntry.length)))
				if curEntry.key != TAZKey:
					curEntry.isTimeAmp = False
			else:
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
			nextEntry.length -= padLength
			#Concatenate the waveforms
			if curEntry.isTimeAmp:
				paddedWF = np.hstack((wfLib[curEntry.key]*np.ones(curEntry.length), np.zeros(padLength, dtype=np.complex)))
				if curEntry.key != TAZKey:
					curEntry.isTimeAmp = False
			else:
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
		assert idx + wf.size < MAX_WAVEFORM_PTS, 'Oops! You have exceeded the waveform memory of the APS'
		wfVec[idx:idx+wf.size] = np.uint16(np.round(MAX_WAVEFORM_VALUE*wf))
		offsets[key] = idx
		idx += wf.size 
					
	#Trim the waveform 
	wfVec = wfVec[0:idx] 

	return wfVec, offsets

def calc_marker_delay(entry):
	#The firmware cannot handle 0 delay markers so push out one clock cycle
	if entry.markerDelay1 is not None:
		if entry.markerDelay1 < ADDRESS_UNIT:
			entry.markerDelay1 = ADDRESS_UNIT
		markerDelay1 = entry.markerDelay1//ADDRESS_UNIT
	else:
		markerDelay1 = 0

	if entry.markerDelay2 is not None:
		if entry.markerDelay2 < ADDRESS_UNIT:
			entry.markerDelay2 = ADDRESS_UNIT
		markerDelay2 = entry.markerDelay2//ADDRESS_UNIT
	else:
		markerDelay2 = 0

	return markerDelay1, markerDelay2


def create_LL_data(LLs, offsets, AWGName=''):
	'''
	Helper function to create LL data vectors from a list of miniLL's and an offset dictionary
	keyed on the wf keys.
	'''

	#Preallocate the bank data and do some checking for miniLL lengths
	seqLengths = np.array([len(miniLL) for miniLL in LLs])
	assert np.all(seqLengths >= MIN_LL_ENTRY_COUNT), 'Oops! mini LL''s needs to have at least two elements.'
	assert np.all(seqLengths < MAX_LL_ENTRIES), 'Oops! mini LL''s cannot have length greater than {0}, you have {1} entries'.format(MAX_BANK_SIZE, len(miniLL))
	numEntries = sum(seqLengths)
	LLData = {label: np.zeros(numEntries, dtype=np.uint16) for label in ['addr','count', 'trigger1', 'trigger2', 'repeat']}

	#Loop over all entries
	TAPairEntries = []
	for ct, entry in enumerate(chain.from_iterable(LLs)):
		LLData['addr'][ct] = offsets[entry.key]//ADDRESS_UNIT
		LLData['count'][ct] = entry.length//ADDRESS_UNIT-1
		LLData['trigger1'][ct], LLData['trigger2'][ct] = calc_marker_delay(entry)
		LLData['repeat'][ct] = entry.repeat-1
		if entry.isTimeAmp:
			TAPairEntries.append(ct)

	#Add in the miniLL start/stop and TA pair flags on the upper bits of the repeat entries
	startPts = np.hstack((0, np.cumsum(seqLengths[:-1])))
	endPts = np.cumsum(seqLengths)-1
	LLData['repeat'][startPts] += 2**START_MINILL_BIT + 2**WAIT_TRIG_BIT
	LLData['repeat'][endPts] += 2**END_MINILL_BIT
	LLData['repeat'][TAPairEntries] += 2**TA_PAIR_BIT
	
	#Check streaming requirements
	numEntries = np.uint16(LLData['addr'].size)
	if numEntries > MAX_LL_ENTRIES:
		print('Streaming will be necessary for {}'.format(AWGName))
		#Get the length of the longest LL
		llLengths = np.sort([len(miniLL) for miniLL in LLs])[-2:]
		if sum(llLengths) > MAX_LL_ENTRIES:
			print('Oops!  It seems the longest two sequences do not fit in memory at the same time. Make sure you know what you are doing.')
		timePerEntry = .050/4096; # measured 46ms average for 4096 entries, use 50 ms as a conservative estimate
		maxRepInterval = timePerEntry*llLengths[1]
		print('Maximum suggested sequence rate is {:.3f}ms, or for 100us rep. rate this would be {} miniLL repeats'.format(1e3*maxRepInterval, int(maxRepInterval/100e-6)))

	return LLData, numEntries

def merge_APS_markerData(IQLL, markerLL, markerNum):
	'''
	Helper function to merge two marker channels into an IQ channel.
	'''

	markerAttr = 'markerDelay' + str(markerNum)

	# expand link lists to the same length (copying first element of shorter one)
	for miniLL_IQ, miniLL_m in izip_longest(IQLL, markerLL):
		if not miniLL_IQ:
			IQLL.append(deepcopy(IQLL[0]))
		if not miniLL_m:
			markerLL.append(deepcopy(markerLL[0]))

	#Step through the all the miniLL's together
	for miniLL_IQ, miniLL_m in zip(IQLL, markerLL):
		#Find the cummulative length for each entry of IQ channel
		timePts = np.cumsum([0] + [entry.totLength for entry in miniLL_IQ])

		#Find the switching points of the marker channels
		switchPts = []
		curIndex = 0
		for curEntry, nextEntry in zip(miniLL_m[:-1], miniLL_m[1:]):
			curIndex += curEntry.totLength
			if curEntry.key != nextEntry.key:
				switchPts.append(curIndex)
				
		#Assume switch pts seperated by 1 point are single trigger blips
		blipPts = (np.diff(switchPts) == 1).nonzero()[0]
		for pt in blipPts[::-1]:
			del switchPts[pt+1]
		#Ensure the IQ LL is long enough to support the blips
		if switchPts:
			if max(switchPts) >= timePts[-1]:
				assert miniLL_IQ[-1].isTimeAmp
				miniLL_IQ[-1].length += max(switchPts) - timePts[-1] + 8 

		#Now map onto linklist elements
		curIQIdx = 0
		trigQueue = []
		for switchPt in switchPts:
			#If the trigger count is too long we need to move to the next IQ entry
			while ((switchPt - timePts[curIQIdx]) > ADDRESS_UNIT * MAX_TRIGGER_COUNT) or (len(trigQueue) > 1):
				# update the trigger queue, dropping triggers that have played
				trigQueue = [t - miniLL_IQ[curIQIdx].length for t in trigQueue]
				trigQueue = [t for t in trigQueue if t >= 0]
				curIQIdx += 1
			#Push on the trigger count
			#If are switch point is before the start of the LL entry then we are in trouble...
			if switchPt - timePts[curIQIdx] <= 0:
				#See if the previous entry was a TA pair and whether we can split it
				needToShift = switchPt - timePts[curIQIdx-1]
				if miniLL_IQ[curIQIdx-1].isTimeAmp and miniLL_IQ[curIQIdx-1].length > needToShift + MIN_ENTRY_LENGTH:
					miniLL_IQ.insert(curIQIdx, deepcopy(miniLL_IQ[curIQIdx-1]))
					miniLL_IQ[curIQIdx-1].length = needToShift-ADDRESS_UNIT
					miniLL_IQ[curIQIdx].length -= needToShift-ADDRESS_UNIT
					miniLL_IQ[curIQIdx].markerDelay1 = None
					miniLL_IQ[curIQIdx].markerDelay2 = None
					setattr(miniLL_IQ[curIQIdx], markerAttr, ADDRESS_UNIT)
					#Recalculate the timePts
					timePts = np.cumsum([0] + [entry.totLength for entry in miniLL_IQ])
				else:
					setattr(miniLL_IQ[curIQIdx], markerAttr, 0)
					print("Had to push marker blip out to start of next entry.")
			else:
				setattr(miniLL_IQ[curIQIdx], markerAttr, switchPt - timePts[curIQIdx])
				trigQueue.insert(0, switchPt - timePts[curIQIdx])
			# update the trigger queue
			trigQueue = [t - miniLL_IQ[curIQIdx].length for t in trigQueue]
			trigQueue = [t for t in trigQueue if t >= 0]
			curIQIdx += 1

	#Replace any remaining empty entries with None
	for miniLL_IQ in IQLL:
		for entry in miniLL_IQ:
			if not hasattr(entry, markerAttr):
				setattr(entry, markerAttr, None)

def write_APS_file(awgData, fileName, miniLLRepeat=1):
	'''
	Main function to pack channel LLs into an APS h5 file.
	'''

	#Preprocess the LL data to handle APS restrictions
	LLs12 = [preprocess_APS(miniLL, awgData['ch12']['wfLib']) for miniLL in awgData['ch12']['linkList']]
	LLs34 = [preprocess_APS(miniLL, awgData['ch34']['wfLib']) for miniLL in awgData['ch34']['linkList']]

	#Merge the the marker data into the IQ linklists
	merge_APS_markerData(LLs12, awgData['ch1m1']['linkList'], 1)
	merge_APS_markerData(LLs12, awgData['ch2m1']['linkList'], 2)
	merge_APS_markerData(LLs34, awgData['ch3m1']['linkList'], 1)
	merge_APS_markerData(LLs34, awgData['ch4m1']['linkList'], 2)
	
	#Open the HDF5 file
	if os.path.isfile(fileName):
		os.remove(fileName)
	with h5py.File(fileName, 'w') as FID:  
	
		#List of which channels we have data for
		#TODO: actually handle incomplete channel data
		channelDataFor = [1,2,3,4]
		FID['/'].attrs['Version'] = 2.0
		FID['/'].attrs['channelDataFor'] = np.uint16(channelDataFor)
		FID['/'].attrs['miniLLRepeat'] = np.uint16(miniLLRepeat - 1)
   
		#Create the waveform vectors
		wfInfo = []
		for wfLib in (awgData['ch12']['wfLib'], awgData['ch34']['wfLib']):
			wfInfo.append(create_wf_vector({key:wf.real for key,wf in wfLib.items()}))
			wfInfo.append(create_wf_vector({key:wf.imag for key,wf in wfLib.items()}))

		LLData = [LLs12, LLs34]
		#Create the groups and datasets
		for chanct in range(4):
			chanStr = '/chan_{0}'.format(chanct+1)
			chanGroup = FID.create_group(chanStr)
			chanGroup.attrs['isIQMode'] = np.uint8(1)
			#Write the waveformLib to file
			FID.create_dataset('{0}/waveformLib'.format(chanStr), data=wfInfo[chanct][0])

			#For A channels (1 & 3) we write link list data
			if np.mod(chanct,2) == 0:
				chanGroup.attrs['isLinkListData'] = np.uint8(1)
				groupStr = chanStr+'/linkListData'
				LLGroup = FID.create_group(groupStr)
				LLDataVecs, numEntries = create_LL_data(LLData[chanct//2], wfInfo[chanct][1], os.path.basename(fileName))
				LLGroup.attrs['length'] = numEntries
				for key,dataVec in LLDataVecs.items():
					FID.create_dataset(groupStr+'/' + key, data=dataVec)
			else:
				chanGroup.attrs['isLinkListData'] = np.uint8(0)


def read_APS_file(fileName):
	'''
	Helper function to read back in data from a H5 file and reconstruct the sequence
	'''
	AWGData = {}
	#APS bit masks
	START_MINILL_MASK = 2**START_MINILL_BIT
	END_MINILL_MASK = 2**END_MINILL_BIT
	TA_PAIR_MASK = 2**TA_PAIR_BIT
	REPEAT_MASK = 2**10-1
	
	chanStrs = ['ch1','ch2', 'ch3', 'ch4']
	chanStrs2 = ['chan_1', 'chan_2', 'chan_3', 'chan_4']
	mrkStrs = ['ch1m1', 'ch2m1', 'ch3m1', 'ch4m1']

	with h5py.File(fileName, 'r') as FID:
		for chanct, chanStr in enumerate(chanStrs2):
			#If we're in IQ mode then the Q channel gets its linkListData from the I channel
			if FID[chanStr].attrs['isIQMode']:
				tmpChan = 2*(chanct//2)
				curLLData = FID[chanStrs2[tmpChan]]['linkListData']
			else:
				curLLData = FID[chanStr]['linkListData']
			#Pull out the LL data
			#Matlab puts our column vectors so need to flatten too
			tmpAddr = curLLData['addr'].value.flatten()
			tmpCount = curLLData['count'].value.flatten()
			tmpRepeat = curLLData['repeat'].value.flatten()
			tmpTrigger1 = curLLData['trigger1'].value.flatten()
			tmpTrigger2 = curLLData['trigger2'].value.flatten()
			numEntries = curLLData.attrs['length']
   
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
					triggerDelays = []
					
				#Record the trigger delays
				if np.mod(chanct,2) == 0:
					if tmpTrigger1[entryct] > 0:
						triggerDelays.append(AWGData[chanStrs[chanct]][-1].size + ADDRESS_UNIT*tmpTrigger1[entryct])
				else:
					if tmpTrigger2[entryct] > 0:
						triggerDelays.append(AWGData[chanStrs[chanct]][-1].size + ADDRESS_UNIT*tmpTrigger2[entryct])

				#If it is a TA pair or regular pulse
				curRepeat = (tmpRepeat[entryct] & REPEAT_MASK)+1
				if TA_PAIR_MASK & tmpRepeat[entryct]:
					AWGData[chanStrs[chanct]][-1] = np.hstack((AWGData[chanStrs[chanct]][-1], 
													np.tile(wfLib[tmpAddr[entryct]*ADDRESS_UNIT:tmpAddr[entryct]*ADDRESS_UNIT+4], curRepeat*(tmpCount[entryct]+1))))
				else:
					AWGData[chanStrs[chanct]][-1] = np.hstack((AWGData[chanStrs[chanct]][-1], 
													np.tile(wfLib[tmpAddr[entryct]*ADDRESS_UNIT:tmpAddr[entryct]*ADDRESS_UNIT+ADDRESS_UNIT*(tmpCount[entryct]+1)], curRepeat)))
				#Add the trigger pulse
				if END_MINILL_MASK & tmpRepeat[entryct]:
					triggerSeq = np.zeros(AWGData[chanStrs[chanct]][-1].size, dtype=np.bool)
					triggerSeq[triggerDelays] = True
					AWGData[mrkStrs[chanct]].append(triggerSeq)
	return AWGData


if __name__ == '__main__':

	pass
