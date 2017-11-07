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
from itertools import chain
from future.moves.itertools import zip_longest
from QGL import Compiler, ControlFlow, BlockLabel, PatternUtils
from QGL.PatternUtils import hash_pulse, flatten
from copy import copy, deepcopy

#Some constants
SAMPLING_RATE = 1.2e9
ADDRESS_UNIT = 4  #everything is done in units of 4 timesteps
MIN_ENTRY_LENGTH = 12
MIN_LL_ENTRY_COUNT = 2  #minimum length of mini link list
MAX_WAVEFORM_PTS = 2**15  #maximum size of waveform memory
MAX_WAVEFORM_VALUE = 2**13 - 1  #maximum waveform value i.e. 14bit DAC
MAX_LL_ENTRIES = 8192  #maximum number of LL entries in a bank
MAX_REPEAT_COUNT = 2**10 - 1
MAX_TRIGGER_COUNT = 2**16 - 1

#APS bit masks
START_MINILL_BIT = 15
END_MINILL_BIT = 14
WAIT_TRIG_BIT = 13
TA_PAIR_BIT = 12

# Do we want a pulse file per instrument or per channel
SEQFILE_PER_CHANNEL = False

def get_empty_channel_set():
    return {'ch12': {},
            'ch34': {},
            'ch1m1': {},
            'ch2m1': {},
            'ch3m1': {},
            'ch4m1': {}}


def get_seq_file_extension():
    return '.h5'


def is_compatible_file(filename):
    with h5py.File(filename, 'r') as FID:
        target = FID['/'].attrs['target hardware']
        if isinstance(target, str):
            target = target.encode('utf-8')
        if target == b'APS1':
            return True
    return False

class APSWaveform(object):
    """
    More specific APS version of a waveform
    """
    def __init__(self, waveform):
        self.label        = waveform.label
        self.key          = waveform.key
        self.amp          = waveform.amp
        self.length       = PatternUtils.convert_length_to_samples(waveform.length, SAMPLING_RATE, ADDRESS_UNIT)
        self.phase        = waveform.phase
        self.frameChange  = waveform.frameChange
        self.isTimeAmp    = waveform.isTimeAmp
        self.frequency    = waveform.frequency
        self.repeat       = 1
        self.markerDelay1 = None
        self.markerDelay2 = None

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self.isTimeAmp:
            TA = 'HIGH' if self.amp != 0 else 'LOW'
            return "APSWaveform-TA(" + TA + ", " + str(self.length) + ")"
        else:
            return "APSWaveform(" + self.label + ", " + str(
                self.key)[:6] + ", " + str(self.length) + ")"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self == other

    @property
    def isZero(self):
        return self.amp == 0


def preprocess(seqs, shapeLib, T):
    for seq in seqs:
        for ct,e in enumerate(seq):
            if isinstance(e, Compiler.Waveform):
                seq[ct] = APSWaveform(e)
    seqs, miniLLrepeat = unroll_loops(seqs)
    for seq in seqs:
        PatternUtils.propagate_frame_changes(seq, wf_type=APSWaveform)
    PatternUtils.quantize_phase(seqs, 1.0 / 2**13, wf_type=APSWaveform)
    compress_sequences(seqs)
    wfLib = build_waveforms(seqs, shapeLib)
    PatternUtils.correct_mixers(wfLib, T)
    for ct in range(len(seqs)):
        seqs[ct] = apply_min_pulse_constraints(seqs[ct], wfLib)
    return seqs, miniLLrepeat, wfLib


def compress_sequences(seqs):
    '''
	Drop zero-length pulses and combine adjacent TA pairs into single entries
	'''
    for seq in seqs:
        ct = 1
        while ct < len(seq):
            prevEntry = seq[ct - 1]
            curEntry = seq[ct]
            if isinstance(curEntry, APSWaveform) and curEntry.length == 0:
                del seq[ct]
            elif isinstance(prevEntry, APSWaveform) and isinstance(curEntry, APSWaveform) and \
               prevEntry.isTimeAmp and curEntry.isTimeAmp and \
               prevEntry.amp == curEntry.amp and \
               prevEntry.phase == curEntry.phase:
                prevEntry.length += curEntry.length
                prevEntry.frameChange += curEntry.frameChange
                del seq[ct]
            ct += 1


def build_waveforms(seqs, shapeLib):
    # apply amplitude, phase, and modulation and add the resulting waveforms to the library
    wfLib = {wf_sig(padding_entry(0)): TAZShape}
    for wf in flatten(seqs):
        if isinstance(wf, APSWaveform) and wf_sig(wf) not in wfLib:
            shape = np.exp(1j * wf.phase) * wf.amp * shapeLib[wf.key]
            if wf.frequency != 0 and wf.amp != 0:
                shape *= np.exp(
                    -1j * 2 * np.pi * wf.frequency * np.arange(wf.length) /
                    SAMPLING_RATE)  #minus from negative frequency qubits
            wfLib[wf_sig(wf)] = shape
    return wfLib


def wf_sig(wf):
    '''
	Compute a signature of a Compiler.Waveform that identifies the relevant properties for
	two Waveforms to be considered "equal" in the waveform library. For example, we ignore
	length of TA waveforms.
	'''
    # 2nd condition necessary until we support RT SSB
    if wf.isZero or (wf.isTimeAmp and wf.frequency == 0 ):
        return (wf.amp, wf.phase)
    else:
        return (wf.key, wf.amp, round(wf.phase * 2**13), wf.length,
                wf.frequency)


TAZShape = np.zeros(1, dtype=np.complex)
TAZKey = hash_pulse(TAZShape)


def padding_entry(length):
    entry = Compiler.Waveform()
    entry.length = length / SAMPLING_RATE
    entry.key = TAZKey
    entry.isTimeAmp = True
    return APSWaveform(entry)


def apply_min_pulse_constraints(miniLL, wfLib):
    '''
	Helper function to deal with LL elements less than minimum LL entry count
	by trying to concatenate them into neighbouring entries
	'''

    newMiniLL = []
    entryct = 0
    while entryct < len(miniLL):
        curEntry = miniLL[entryct]
        if not isinstance(curEntry, APSWaveform) or \
                curEntry.length >= MIN_ENTRY_LENGTH:
            newMiniLL.append(curEntry)
            entryct += 1
            continue

        if entryct == len(miniLL) - 1:
            # we've run out of entries to append to. drop it?
            warn("Unable to handle too short LL element, dropping.")
            break
        nextEntry = miniLL[entryct + 1]
        previousEntry = miniLL[entryct - 1] if entryct > 0 else None

        # For short TA pairs we see if we can add them to the next waveform
        if curEntry.isZero and not nextEntry.isZero:
            # Concatenate the waveforms
            paddedWF = np.hstack((np.zeros(curEntry.length,
                                           dtype=np.complex),
                                  wfLib[wf_sig(nextEntry)]))
            # Generate a new key
            nextEntry.key = hash_pulse(paddedWF)
            nextEntry.length = paddedWF.size
            wfLib[wf_sig(nextEntry)] = paddedWF
            newMiniLL.append(nextEntry)
            entryct += 2

        # For short pulses we see if we can steal some padding from the previous or next entry
        elif isinstance(previousEntry, APSWaveform) and \
                previousEntry.isZero and \
                previousEntry.length > 2 * MIN_ENTRY_LENGTH:
            padLength = MIN_ENTRY_LENGTH - curEntry.length
            newMiniLL[-1].length -= padLength
            # Concatenate the waveforms
            if curEntry.isZero:
                curEntry.length += padLength
                entryct += 1
                curEntry.isTimeAmp = True
                continue
            elif curEntry.isTimeAmp:  # non-zero
                paddedWF = np.hstack(
                    (np.zeros(padLength, dtype=np.complex),
                     wfLib[wf_sig(curEntry)] * np.ones(curEntry.length)))
                curEntry.isTimeAmp = False
            else:
                paddedWF = np.hstack((np.zeros(padLength,
                                               dtype=np.complex),
                                      wfLib[wf_sig(curEntry)]))
            # Generate a new key
            curEntry.key = hash_pulse(paddedWF)
            curEntry.length = paddedWF.size
            wfLib[wf_sig(curEntry)] = paddedWF
            newMiniLL.append(curEntry)
            entryct += 1

        elif isinstance(nextEntry, APSWaveform) and \
                nextEntry.isZero and \
                nextEntry.length > 2 * MIN_ENTRY_LENGTH:
            padLength = MIN_ENTRY_LENGTH - curEntry.length
            nextEntry.length -= padLength
            # Concatenate the waveforms
            if curEntry.isZero:
                curEntry.length += padLength
                entryct += 1
                curEntry.isTimeAmp = True
                continue
            elif curEntry.isTimeAmp:  #non-zero
                paddedWF = np.hstack(
                    (wfLib[curEntry.key] * np.ones(curEntry.length),
                     np.zeros(padLength, dtype=np.complex)))
                curEntry.isTimeAmp = False
            else:
                paddedWF = np.hstack((wfLib[curEntry.key],
                                      np.zeros(padLength,
                                               dtype=np.complex)))
            # Generate a new key
            curEntry.key = hash_pulse(paddedWF)
            curEntry.length = paddedWF.size
            wfLib[wf_sig(curEntry)] = paddedWF
            newMiniLL.append(curEntry)
            entryct += 1

        else:
            warn("Unable to handle too short LL element, dropping.")
            entryct += 1

    # Update the miniLL
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
        wf[wf > 1] = 1.0
        wf[wf < -1] = -1.0
        #TA pairs need to be repeated ADDRESS_UNIT times
        if wf.size == 1:
            wf = wf.repeat(ADDRESS_UNIT)
        #Ensure the wf is an integer number of ADDRESS_UNIT's
        trim = wf.size % ADDRESS_UNIT
        if trim:
            wf = wf[:-trim]
        assert idx + wf.size < MAX_WAVEFORM_PTS, 'Oops! You have exceeded the waveform memory of the APS'
        wfVec[idx:idx + wf.size] = np.uint16(np.round(MAX_WAVEFORM_VALUE * wf))
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
        markerDelay1 = entry.markerDelay1 // ADDRESS_UNIT
    else:
        markerDelay1 = 0

    if entry.markerDelay2 is not None:
        if entry.markerDelay2 < ADDRESS_UNIT:
            entry.markerDelay2 = ADDRESS_UNIT
        markerDelay2 = entry.markerDelay2 // ADDRESS_UNIT
    else:
        markerDelay2 = 0

    return markerDelay1, markerDelay2


class Instruction(object):
    def __init__(self, addr=0, count=0, trig1=0, trig2=0, repeat=0):
        self.addr = int(addr)
        self.count = int(count)
        self.trig1 = int(trig1)
        self.trig2 = int(trig2)
        self.repeat = int(repeat)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return ("Instruction(" + str(self.addr) + ", " + str(self.count) + ", "
                + str(self.trig1) + ", " + str(self.trig2) + ", " +
                str(self.repeat) + ")")

    @property
    def start(self):
        return self.repeat & (1 << START_MINILL_BIT)

    @start.setter
    def start(self, value):
        self.repeat |= (value & 0x1) << START_MINILL_BIT

    @property
    def end(self):
        return self.repeat & (1 << END_MINILL_BIT)

    @end.setter
    def end(self, value):
        self.repeat |= (value & 0x1) << END_MINILL_BIT

    @property
    def wait(self):
        return self.repeat & (1 << WAIT_TRIG_BIT)

    @wait.setter
    def wait(self, value):
        self.repeat |= (value & 0x1) << WAIT_TRIG_BIT

    @property
    def TAPair(self):
        return self.repeat & (1 << TA_PAIR_BIT)

    @TAPair.setter
    def TAPair(self, value):
        self.repeat |= (value & 0x1) << TA_PAIR_BIT

    def flatten(self):
        return (self.addr << 16 * 4) | (self.count << 16 * 3) | (
            self.trig1 << 16 * 2) | (self.trig2 << 16 * 1) | self.repeat


def create_LL_data(LLs, offsets, AWGName=''):
    '''
	Helper function to create LL data vectors from a list of miniLL's and an offset dictionary
	keyed on the wf keys.
	'''

    # Do some checking on miniLL lengths
    seqLengths = np.array([len(miniLL) for miniLL in LLs])
    assert np.all(
        seqLengths <
        MAX_LL_ENTRIES), 'Oops! mini LL' 's cannot have length greater than {0}, you have {1} entries'.format(
            MAX_BANK_SIZE, len(miniLL))

    for miniLL in LLs:
        # add one because we need at least one control instruction (WAIT) plus MIN_LL_ENTRY_COUNT waveforms
        while len(miniLL) < MIN_LL_ENTRY_COUNT + 1:
            miniLL.append(padding_entry(MIN_ENTRY_LENGTH))

    instructions = []
    waitFlag = False
    for miniLL in LLs:
        miniStart = True
        for entry in miniLL:
            if isinstance(entry, BlockLabel.BlockLabel):
                # don't emit any instructions for labels
                continue
            elif isinstance(entry, ControlFlow.ControlInstruction):
                if isinstance(entry, ControlFlow.Wait):
                    waitFlag = True
                    continue
                elif isinstance(
                        entry, ControlFlow.Goto) and entry.target == LLs[0][0]:
                    # can safely skip a goto with a target of the first instruction
                    continue
                else:
                    warn("skipping instruction {0}".format(entry))
            else:  # waveform instructions
                t1, t2 = calc_marker_delay(entry)
                instr = Instruction(addr=offsets[wf_sig(entry)] //
                                    ADDRESS_UNIT,
                                    count=entry.length // ADDRESS_UNIT - 1,
                                    trig1=t1,
                                    trig2=t2,
                                    repeat=entry.repeat - 1)
                # set flags
                instr.TAPair = entry.isTimeAmp or entry.isZero
                instr.wait = waitFlag
                instr.start = miniStart
                waitFlag = False
                miniStart = False
                instructions.append(instr)
        instructions[-1].end = True

    # convert to LLData structure
    numEntries = len(instructions)
    LLData = {label: np.zeros(numEntries, dtype=np.uint16)
              for label in ['addr', 'count', 'trigger1', 'trigger2', 'repeat']}
    for ct in range(numEntries):
        LLData['addr'][ct] = instructions[ct].addr
        LLData['count'][ct] = instructions[ct].count
        LLData['trigger1'][ct] = instructions[ct].trig1
        LLData['trigger2'][ct] = instructions[ct].trig2
        LLData['repeat'][ct] = instructions[ct].repeat
    #Check streaming requirements
    if numEntries > MAX_LL_ENTRIES:
        print('Streaming will be necessary for {}'.format(AWGName))
        #Get the length of the longest LL
        llLengths = np.sort([len(miniLL) for miniLL in LLs])[-2:]
        if sum(llLengths) > MAX_LL_ENTRIES:
            print(
                'Oops!  It seems the longest two sequences do not fit in memory at the same time. Make sure you know what you are doing.')
        timePerEntry = .050 / 4096
        # measured 46ms average for 4096 entries, use 50 ms as a conservative estimate
        maxRepInterval = timePerEntry * llLengths[1]
        print(
            'Maximum suggested sequence rate is {:.3f}ms, or for 100us rep. rate this would be {} miniLL repeats'.format(
                1e3 * maxRepInterval, int(maxRepInterval / 100e-6)))

    return LLData, numEntries


def merge_APS_markerData(IQLL, markerLL, markerNum):
    '''
	Helper function to merge two marker channels into an IQ channel.
	'''
    if len(markerLL) == 0:
        return
    assert len(IQLL) <= len(markerLL), "Sequence length mismatch"
    if len(IQLL) < len(markerLL):
        for ct in range(len(markerLL) - len(IQLL)):
            IQLL.append([])

    for seq in markerLL:
        PatternUtils.convert_lengths_to_samples(seq, SAMPLING_RATE,
                                                ADDRESS_UNIT, Compiler.Waveform)

    markerAttr = 'markerDelay' + str(markerNum)

    #Step through the all the miniLL's together
    for miniLL_IQ, miniLL_m in zip_longest(IQLL, markerLL):
        #Find the switching points of the marker channels
        switchPts = []
        prevAmplitude = 0
        t = 0
        for entry in miniLL_m:
            if hasattr(entry, 'amp') and prevAmplitude != entry.amp:
                switchPts.append(t)
                prevAmplitude = entry.amp
            t += entry.length

        if len(switchPts) == 0:
            # need at least a WAIT on an empty IQ LL in order to match segment sequencing
            if len(miniLL_IQ) == 0:
                miniLL_IQ.append(ControlFlow.qwait())
            continue

        # Push on an extra switch point if we have an odd number of switches (to maintain state)
        if len(switchPts) % 2 == 1:
            switchPts.append(t)

        #Assume switch pts seperated by 0 or 1 point are single trigger blips
        blipPts = (np.diff(switchPts) <= 1).nonzero()[0]
        for pt in blipPts[::-1]:
            del switchPts[pt + 1]

        # if the IQ sequence is empty, make an ideally length-matched sequence
        if len(miniLL_IQ) == 0:
            miniLL_IQ.append(ControlFlow.qwait())
            miniLL_IQ.append(padding_entry(max(switchPts[0],
                                               MIN_ENTRY_LENGTH)))
            for length in np.diff(switchPts):
                miniLL_IQ.append(padding_entry(max(length, MIN_ENTRY_LENGTH)))

        #Find the cummulative length for each entry of IQ channel
        timePts = np.cumsum([0] + [entry.length for entry in miniLL_IQ])

        #Ensure the IQ LL is long enough to support the blips
        if max(switchPts) >= timePts[-1]:
            dt = max(switchPts) - timePts[-1]
            if hasattr(miniLL_IQ[-1], 'isTimeAmp') and miniLL_IQ[-1].isTimeAmp:
                miniLL_IQ[-1].length += dt + 4
            else:
                # inject before any control flow statements at the end of the sequence
                idx = len(miniLL_IQ)
                while idx > 0 and isinstance(miniLL_IQ[idx - 1],
                                             ControlFlow.ControlInstruction):
                    idx -= 1
                miniLL_IQ.insert(idx,
                                 padding_entry(max(dt + 4, MIN_ENTRY_LENGTH)))

        #Now map onto linklist elements
        curIQIdx = 0
        trigQueue = []
        for ct, switchPt in enumerate(switchPts):
            # skip if:
            #   1) control-flow instruction or label (i.e. not a waveform)
            #   2) the trigger count is too long
            #   3) the previous trigger pulse entends into the current entry
            while (not isinstance(miniLL_IQ[curIQIdx], APSWaveform) or
                   (switchPt - timePts[curIQIdx]) >
                   (ADDRESS_UNIT * MAX_TRIGGER_COUNT) or len(trigQueue) > 1):
                # update the trigger queue, dropping triggers that have played
                trigQueue = [t - miniLL_IQ[curIQIdx].length for t in trigQueue]
                trigQueue = [t for t in trigQueue if t >= 0]
                curIQIdx += 1
                # add padding pulses if needed
                if curIQIdx >= len(miniLL_IQ):
                    if len(trigQueue) > 0:
                        pad = max(MIN_ENTRY_LENGTH, min(trigQueue, 0))
                    else:
                        pad = MIN_ENTRY_LENGTH
                    miniLL_IQ.append(padding_entry(pad))
            #Push on the trigger count

            #If our switch point is before the start of the LL entry then we are in trouble...
            if switchPt - timePts[curIQIdx] < 0:
                #See if the previous entry was a TA pair and whether we can split it
                needToShift = switchPt - timePts[curIQIdx - 1]
                assert needToShift > MIN_ENTRY_LENGTH + ADDRESS_UNIT, "Sequential marker blips too close together."
                if isinstance(miniLL_IQ[curIQIdx-1], APSWaveform) and \
                 miniLL_IQ[curIQIdx-1].isTimeAmp and \
                 miniLL_IQ[curIQIdx-1].length > (needToShift + MIN_ENTRY_LENGTH):

                    miniLL_IQ.insert(curIQIdx,
                                     deepcopy(miniLL_IQ[curIQIdx - 1]))
                    miniLL_IQ[curIQIdx - 1].length = needToShift - ADDRESS_UNIT
                    miniLL_IQ[curIQIdx].length -= needToShift - ADDRESS_UNIT
                    miniLL_IQ[curIQIdx].markerDelay1 = None
                    miniLL_IQ[curIQIdx].markerDelay2 = None
                    setattr(miniLL_IQ[curIQIdx], markerAttr, ADDRESS_UNIT)
                    #Recalculate the timePts
                    timePts = np.cumsum([0] + [entry.length
                                               for entry in miniLL_IQ])
                else:
                    setattr(miniLL_IQ[curIQIdx], markerAttr, 0)
                    print(
                        "Had to push marker blip out to start of next entry.")

            else:
                setattr(miniLL_IQ[curIQIdx], markerAttr,
                        switchPt - timePts[curIQIdx])
                trigQueue.insert(0, switchPt - timePts[curIQIdx])
            # update the trigger queue
            trigQueue = [t - miniLL_IQ[curIQIdx].length for t in trigQueue]
            trigQueue = [t for t in trigQueue if t >= 0]
            curIQIdx += 1

            # add padding pulses if needed
            if ct + 1 < len(switchPts) and curIQIdx >= len(miniLL_IQ):
                if len(trigQueue) > 0:
                    pad = max(MIN_ENTRY_LENGTH, min(trigQueue, 0))
                else:
                    pad = MIN_ENTRY_LENGTH
                miniLL_IQ.append(padding_entry(pad))


def unroll_loops(LLs):
    '''
	Unrolls repeated sequences in place, unless the sequence can be unrolled with a miniLL repeat
	attribute. Returns the (potentially) modified sequence and the miniLL repeat value.
	'''
    # if all sequences start and end with LOAD and REPEAT, respectively, and all load values
    # are the same, we can just drop these instructions and return a miniLLrepeat value
    if not LLs or not LLs[0] or not LLs[0][0]:
        return LLs, 0
    elif isinstance(LLs[0][0], ControlFlow.LoadRepeat):
        repeats = LLs[0][0].value
    else:
        repeats = -1

    simpleUnroll = True
    for seq in LLs:
        if not isinstance(seq[0], ControlFlow.LoadRepeat) or \
         not isinstance(seq[-1], ControlFlow.Repeat) or \
         seq[0].value != repeats or \
         seq[-1].target != seq[1]:
            simpleUnroll = False

    if simpleUnroll:
        return LLs, repeats

    # otherwise, we need to manually unroll any repeated section
    instructions = []
    for seq in LLs:
        symbols = {}
        ct = 0
        while ct < len(seq):
            entry = seq[ct]
            # fill symbol table
            if isinstance(entry, BlockLabel.BlockLabel) and \
             entry not in symbols:
                symbols[entry] = ct
            # look for the end of a repeated block
            if isinstance(entry, ControlFlow.Repeat):
                repeatedBlock = seq[symbols[entry.target] + 1:ct]
                numRepeats = seq[symbols[entry.target] - 1].value
                # unroll the block (dropping the LOAD and REPEAT)
                if len(repeatedBlock) == 1:
                    repeatedBlock[0].repeat = numRepeats
                    seq[symbols[entry.target] - 1:ct + 1] = repeatedBlock
                    # dropped two instructions and a label
                    ct -= 3
                else:
                    seq[symbols[entry.target] - 1:ct +
                        1] = repeatedBlock * numRepeats
                    # advance the count (minus 3 for dropped instructions and label)
                    ct += (numRepeats - 1) * len(repeatedBlock) - 3
            ct += 1
        # add unrolled sequence to instruction list
        instructions.append(seq)
    return instructions, 0


def write_sequence_file(awgData, fileName, miniLLRepeat=1):
    '''
	Main function to pack channel LLs into an APS h5 file.
	'''
    #Preprocess the sequence data to handle APS restrictions
    LLs12, repeat12, wfLib12 = preprocess(awgData['ch12']['linkList'],
                                          awgData['ch12']['wfLib'],
                                          awgData['ch12']['correctionT'])
    LLs34, repeat34, wfLib34 = preprocess(awgData['ch34']['linkList'],
                                          awgData['ch34']['wfLib'],
                                          awgData['ch34']['correctionT'])
    assert repeat12 == repeat34, 'Failed to unroll sequence'
    if repeat12 != 0:
        miniLLRepeat *= repeat12

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
        channelDataFor = [1, 2] if LLs12 else []
        channelDataFor += [3, 4] if LLs34 else []
        FID['/'].attrs['Version'] = 2.2
        FID['/'].attrs['target hardware'] = 'APS1'
        FID['/'].attrs['channelDataFor'] = np.uint16(channelDataFor)
        FID['/'].attrs['miniLLRepeat'] = np.uint16(miniLLRepeat - 1)

        #Create the waveform vectors
        wfInfo = []
        for wfLib in (wfLib12, wfLib34):
            wfInfo.append(create_wf_vector({key: wf.real
                                            for key, wf in wfLib.items()}))
            wfInfo.append(create_wf_vector({key: wf.imag
                                            for key, wf in wfLib.items()}))

        LLData = [LLs12, LLs34]
        repeats = [0, 0]
        #Create the groups and datasets
        for chanct in range(4):
            chanStr = '/chan_{0}'.format(chanct + 1)
            chanGroup = FID.create_group(chanStr)
            chanGroup.attrs['isIQMode'] = np.uint8(1)
            #Write the waveformLib to file
            FID.create_dataset('{0}/waveformLib'.format(chanStr),
                               data=wfInfo[chanct][0])

            #For A channels (1 & 3) we write link list data if we actually have any
            if (np.mod(chanct, 2) == 0) and LLData[chanct // 2]:
                groupStr = chanStr + '/linkListData'
                LLGroup = FID.create_group(groupStr)
                LLDataVecs, numEntries = create_LL_data(
                    LLData[chanct // 2], wfInfo[chanct][1],
                    os.path.basename(fileName))
                LLGroup.attrs['length'] = numEntries
                for key, dataVec in LLDataVecs.items():
                    FID.create_dataset(groupStr + '/' + key, data=dataVec)
            else:
                chanGroup.attrs['isLinkListData'] = np.uint8(0)


def read_sequence_file(fileName):
    """ Helper function to read back in data from a H5 file and reconstruct the
	sequence as a list of (time, amplitude) pairs.
	"""
    AWGData = {}
    #APS bit masks
    START_MINILL_MASK = 2**START_MINILL_BIT
    END_MINILL_MASK = 2**END_MINILL_BIT
    TA_PAIR_MASK = 2**TA_PAIR_BIT
    REPEAT_MASK = 2**10 - 1

    chanStrs = ['ch1', 'ch2', 'ch3', 'ch4']
    chanStrs2 = ['chan_1', 'chan_2', 'chan_3', 'chan_4']
    mrkStrs = ['ch1m1', 'ch2m1', 'ch3m1', 'ch4m1']

    with h5py.File(fileName, 'r') as FID:
        for chanct, chanStr in enumerate(chanStrs2):
            #If we're in IQ mode then the Q channel gets its linkListData from the I channel
            if FID[chanStr].attrs['isIQMode']:
                tmpChan = 2 * (chanct // 2)
                curLLData = FID[chanStrs2[tmpChan]][
                    'linkListData'] if "linkListData" in FID[chanStrs2[
                        tmpChan]] else []
            else:
                curLLData = FID[chanStr][
                    'linkListData'] if "linkListData" in FID[chanStrs2[
                        tmpChan]] else []

            if curLLData:
                #Pull out the LL data in sample units
                #Cast type to avoid uint16 overflow
                addr = (curLLData['addr'].value.astype(np.uint)) * ADDRESS_UNIT
                count = ((curLLData['count'].value + 1).astype(np.uint)
                         ) * ADDRESS_UNIT
                repeat = curLLData['repeat'].value
                trigger1 = (
                    curLLData['trigger1'].value.astype(np.uint)) * ADDRESS_UNIT
                trigger2 = (
                    curLLData['trigger2'].value.astype(np.uint)) * ADDRESS_UNIT

                #Pull out and scale the waveform data
                wf_lib = (
                    1.0 /
                    MAX_WAVEFORM_VALUE) * FID[chanStr]['waveformLib'].value

                #Initialize the lists of time-amplitude pairs
                AWGData[chanStrs[chanct]] = []
                AWGData[mrkStrs[chanct]] = []

                cum_time = 0

                #Loop over LL entries
                for ct in range(curLLData.attrs['length']):
                    #If we are starting a new sequence push back an empty array
                    if START_MINILL_MASK & repeat[ct]:
                        AWGData[chanStrs[chanct]].append([])
                        trigger_delays = [0]
                        cum_time = 0

                    #Record the trigger delays
                    if np.mod(chanct, 2) == 0:
                        if trigger1[ct] > 0:
                            trigger_delays.append(cum_time + trigger1[ct])
                    else:
                        if trigger2[ct] > 0:
                            trigger_delays.append(cum_time + trigger2[ct])

                    #waveforms
                    wf_repeat = (repeat[ct] & REPEAT_MASK) + 1
                    if TA_PAIR_MASK & repeat[ct]:
                        AWGData[chanStrs[chanct]][-1].append(
                            (wf_repeat * count[ct], wf_lib[addr[ct]]))
                    else:
                        for repct in range(wf_repeat):
                            for sample in wf_lib[addr[ct]:addr[ct] + count[
                                    ct]]:
                                AWGData[chanStrs[chanct]][-1].append(
                                    (1, sample))

                    cum_time += count[ct]

                    #Create the trigger sequence
                    if END_MINILL_MASK & repeat[ct]:
                        AWGData[mrkStrs[chanct]].append([])
                        for delay in np.diff(trigger_delays):
                            AWGData[mrkStrs[chanct]][-1].append((delay - 1, 0))
                            AWGData[mrkStrs[chanct]][-1].append((1, 1))

    return AWGData


if __name__ == '__main__':

    pass
