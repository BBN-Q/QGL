'''
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
import Compiler
from warnings import warn
from APSPattern import MIN_ENTRY_LENGTH
from PulseSequencer import Pulse
from math import pi

def delay(linkList, delay, samplingRate):
    '''
    Delays a mini link list by the given delay. Postives delays
    shift right, negative delays shift left.
    '''
    # we assume that link lists are constructed with a padding element at the beginning that can be adjusted
    assert linkList[0][0].isZero, 'Delay error: link list does not start with a padding element'
    sampShift = int(round(delay * samplingRate))
    for miniLL in linkList:
        miniLL[0].length += sampShift
        assert miniLL[0].length > 0, 'Delay error: Negative length padding element after delay.'

def apply_SSB(linkList, wfLib, SSBFreq, samplingRate):
    #Negative because of negative frequency qubits
    phaseStep = -2*pi*SSBFreq/samplingRate
        
    #Bits of phase precision
    #Choose usual DAC vertical precision arbirarily
    phasePrecision = 2**14
    def round_phase(phase, precision):
        """
        Helper function to round a phase to a certain binary precision.
        """
        #Convert radians to portion of circle and then to integer precision round to precision
        intPhase = round(phasePrecision*np.mod(phase/2.0/pi,1))
        return int(intPhase), 2*pi*(intPhase/phasePrecision)

    #Keep a dictionary of pulses and phases
    pulseDict = {}
    for miniLL in linkList:
        curFrame = 0.0
        for entry in miniLL:
            #If it's a zero then just adjust the frame and move on 
            if entry.key == Compiler.TAZKey: 
                curFrame += phaseStep*entry.length
                continue
            # expand time-amplitude pulses in-place
            if entry.isTimeAmp:
                entry.isTimeAmp = False
                shape = wfLib[entry.key][0] * np.ones(entry.length, dtype=np.complex)
            else:
                shape = np.copy(wfLib[entry.key])

            intPhase, truncPhase = round_phase(curFrame, 14)
            pulseTuple = (entry.key, intPhase, entry.length)
            if pulseTuple in pulseDict:
                entry.key = pulseDict[pulseTuple]
            else:
                phaseRamp = phaseStep*np.arange(0.5, shape.size)
                shape *= np.exp(1j*(truncPhase + phaseRamp))
                shapeHash = Compiler.hash_pulse(shape)
                if shapeHash not in wfLib:
                    wfLib[shapeHash] = shape
                pulseDict[pulseTuple] = shapeHash
                entry.key = shapeHash
            curFrame += phaseStep*entry.length

def align(linkList, mode, length):
    for miniLL in linkList:
        miniLL_length = sum([entry.totLength for entry in miniLL])
        paddingEntry = Compiler.create_padding_LL(length - miniLL_length)
        if mode == 'left':
            miniLL.append(paddingEntry)
        elif mode == 'right':
            miniLL.insert(0, paddingEntry)
        else:
            raise NameError("Unknown aligment mode")

def correctMixer(wfLib, T):
    for k, v in wfLib.items():
        # To get the broadcast to work in numpy, need to do the multiplication one row at a time
        iqWF = np.vstack((np.real(v), np.imag(v)))
        wfLib[k] = T[0,:].dot(iqWF) + 1j*T[1,:].dot(iqWF)

def split_multiple_triggers():
	'''
	Split entries with multiple triggers into two entries.
	'''
	pass

def create_gate_seqs(linkList, gateBuffer=0, gateMinWidth=0, samplingRate=1.2e9):
    '''
    Helper function that takes a set of analog channel LL and creates a LL with appropriate 
    blanking on a marker channel. 
    '''

    # convert times into samples
    gateBuffer = int(round(gateBuffer * samplingRate))
    gateMinWidth = int(round(gateMinWidth * samplingRate))
    
    #Initialize list of sequences to return
    gateSeqs = []

    # Time from end of previous LL entry that trigger needs to go
    # high to gate pulse
    startDelay = gateBuffer
    for miniLL in linkList:
        #Initialize a zero-length padding sequence
        gateSeq = [Compiler.create_padding_LL(0)]
        # we need to pad the miniLL with an extra entry if the last entry is not a zero
        if not miniLL[-1].isZero:
            miniLL.append(Compiler.create_padding_LL(MIN_ENTRY_LENGTH))

        #Step through sequence changing state as necessary
        blankHigh = False
        for entry in miniLL:
            #If we are low and the current entry is high then we need to add an element
            if not blankHigh and not entry.isZero:
                gateSeq.append(Compiler.create_padding_LL(entry.totLength, high=True))
                blankHigh = True
            #If we are high and the next entry is low then we need to add an element
            elif blankHigh and entry.isZero:
                gateSeq.append(Compiler.create_padding_LL(entry.totLength, high=False))
                blankHigh = False
            #Otherwise we just continue along in the same state
            else:
                gateSeq[-1].length += entry.totLength


        #Go back through and add the gate buffer to the start of each marker high period.
        #Assume that we start low and alternate low-high-low from the construction above
        #Step through every high pulse and look at the previous one
        for entryct in range(1, len(gateSeq),2):
            #If the previous low pulse is less than the gate buffer then we'll drop it
            #and add its length and the length of the current high entry to the previous high entry 
            if gateSeq[entryct-1].length < gateBuffer:
                #Look for the last valid previous high entry
                goodIdx = entryct-2
                while not gateSeq[goodIdx]:
                    goodIdx -= 2
                gateSeq[goodIdx].length += \
                        gateSeq[entryct-1].totLength + gateSeq[entryct].totLength
                #Mark the two dropped entries as removed by setting them to none
                gateSeq[entryct-1] = None
                gateSeq[entryct] = None
            #Otherwise we subtract the gate buffer from the previous length 
            else:
                gateSeq[entryct-1].length -= gateBuffer
                gateSeq[entryct].length += gateBuffer
            entryct += 2
        #Remove dropped entries
        gateSeq = filter(lambda x : x is not None, gateSeq)

        #Loop through again and make sure that all the low points between pulses are sufficiently long
        #Given the above construction we should have the low-high-low form
        for entryct in range(2, len(gateSeq)-1, 2):
            if gateSeq[entryct].length < gateMinWidth:
                #Consolidate this and the next entry onto the previous one
                #Look for the last valid previous high entry
                goodIdx = entryct-1
                while not gateSeq[goodIdx]:
                    goodIdx -= 2
                gateSeq[goodIdx].length += \
                        gateSeq[entryct].totLength + gateSeq[entryct+1].totLength
                #Mark the two dropped entries as removed by setting them to none
                gateSeq[entryct] = None
                gateSeq[entryct+1] = None
            entryct = 2
        #Remove dropped entries
        gateSeq = filter(lambda x : x is not None, gateSeq)
    
        #Add it on
        gateSeqs.append(gateSeq)

    return gateSeqs

def add_marker_pulse(LL, startPt, length):
    '''
    Helper function to add a marker pulse to a LL from a given startPt and with a given length
    '''
    pass

def add_digitizer_trigger(seqs, trigChan):
    '''
    Add the digitizer trigger.  For now hardcoded but should be loaded from config file.
    '''
    #Assume that last pulse is the measurment pulse for now and tensor on the digitizer trigger pulse
    for seq in seqs:
        #Hack around copied list elements referring to same sequence element
        if isinstance(seq[-1], Pulse) or trigChan not in seq[-1].pulses.keys():
            seq[-1] *= Pulse("digTrig", trigChan, trigChan.pulseParams['shapeFun'](**trigChan.pulseParams), 0.0, 0.0)

def slave_trigger(numSeqs):
    """
    Create slave trigger link lists.
    """
    return [[Compiler.create_padding_LL(1), Compiler.create_padding_LL(1, True), Compiler.create_padding_LL(1)] for _ in range(numSeqs)], Compiler.markerWFLib


    
        



