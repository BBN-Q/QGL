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
from warnings import warn
from PulseSequencer import Pulse, TAPulse
from PulsePrimitives import BLANK
import ControlFlow, Compiler
from math import pi
import hashlib

def hash_pulse(shape):
    return hashlib.sha1(shape.tostring()).hexdigest()

TAZKey = hash_pulse(np.zeros(1, dtype=np.complex))
markerHighKey = hash_pulse(np.ones(1, dtype=np.bool))

def delay(linkList, delay, samplingRate):
    '''
    Delays a mini link list by the given amount.
    '''
    sampShift = int(round(delay * samplingRate))
    if sampShift <= 0: # no need to inject zero delays
        return
    for miniLL in linkList:
        # loop through and look for WAIT instructions
        # use while loop because len(miniLL) will change as we inject delays
        ct = 0
        while ct < len(miniLL):
            if miniLL[ct] == ControlFlow.Wait() or miniLL[ct] == ControlFlow.Sync():
                miniLL.insert(ct+1, Compiler.create_padding_LL(sampShift))
            ct += 1

def normalize_delays(delays):
    '''
    Normalizes a dictionary of channel delays. Postives delays shift right, negative delays shift left.
    Since we cannot delay by a negative amount in hardware, shift all delays until they are positive.
    Takes in a dict of channel:delay pairs and returns a normalized copy of the same.
    '''
    min_delay = min(delays.values())
    out = dict(delays) # copy before modifying
    if min_delay < 0:
        for chan in delays.keys():
            out[chan] += -min_delay
    return out

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
            if not hasattr(entry, 'key') or entry.key == TAZKey:
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
                shapeHash = hash_pulse(shape)
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

def add_gate_pulses(seqs):
    '''
    add gating pulses to Qubit pulses
    '''
    for seq in seqs:
        for ct in range(len(seq)):
            if hasattr(seq[ct], 'pulses'):
                for chan, pulse in seq[ct].pulses.items():
                    if has_gate(chan) and not pulse.isZero and not (chan.gateChan in seq[ct].pulses.keys()):
                        seq[ct] *= BLANK(chan, pulse.length)
            elif hasattr(seq[ct], 'qubits'):
                chan = seq[ct].qubits
                if has_gate(chan) and not seq[ct].isZero:
                    seq[ct] *= BLANK(chan, seq[ct].length)

def has_gate(channel):
    return hasattr(channel, 'gateChan') and channel.gateChan

def apply_gating_constraints(chan, linkList):
    # get channel parameters in samples
    gateBuffer = int(round(chan.gateBuffer * chan.samplingRate))
    gateMinWidth = int(round(chan.gateMinWidth * chan.samplingRate))

    #Initialize list of sequences to return
    gateSeqs = []

    for miniLL in linkList:
        gateSeq = []
        # first pass consolidates entries
        previousEntry = None
        for entry in miniLL:
            if isinstance(entry, ControlFlow.ControlInstruction):
                if previousEntry:
                    gateSeq.append(previousEntry)
                    previousEntry = None
                gateSeq.append(entry)
                continue

            if previousEntry is None:
                previousEntry = entry
                continue

            # matching entry types can be globbed together
            if previousEntry.isZero == entry.isZero:
                previousEntry.length += entry.length
            else:
                gateSeq.append(previousEntry)
                previousEntry = entry

        # push on the last entry if necessary
        if previousEntry:
            gateSeq.append(previousEntry)

        # second pass expands non-zeros by gateBuffer
        for ct in range(len(gateSeq)):
            if isNonZeroWaveform(gateSeq[ct]):
                gateSeq[ct].length += gateBuffer
                # contract the next pulse by the same amount
                if ct + 1 < len(gateSeq) - 1 and not isinstance(gateSeq[ct+1], ControlFlow.ControlInstruction):
                    gateSeq[ct+1].length -= gateBuffer #TODO: what if this becomes negative?

        # third pass ensures gateMinWidth
        ct = 0
        while ct+2 < len(gateSeq):
            # look for pulse, delay, pulse pattern and ensure delay is long enough
            if [isNonZeroWaveform(x) for x in gateSeq[ct:ct+3]] == [True, False, True] and \
                gateSeq[ct+1].length < gateMinWidth and \
                [isinstance(x, ControlFlow.ControlInstruction) for x in gateSeq[ct:ct+3]] == [False, False, False]:
                gateSeq[ct].length += gateSeq[ct+1].length + gateSeq[ct+2].length
                del gateSeq[ct+1:ct+3]
            else:
                ct += 1

        gateSeqs.append(gateSeq)

    return gateSeqs

def isNonZeroWaveform(entry):
    return not isinstance(entry, ControlFlow.ControlInstruction) and not entry.isZero

def add_digitizer_trigger(seqs, trigChan):
    '''
    Add the digitizer trigger to a logical LL (pulse blocks).
    '''
    # Attach a trigger to any pulse block containing a measurement
    pulseLength = trigChan.pulseParams['length'] * trigChan.physChan.samplingRate
    for seq in seqs:
        for ct in range(len(seq)):
            if contains_measurement(seq[ct]) and not (hasattr(seq[ct], 'pulses') and trigChan in seq[ct].pulses.keys()):
                seq[ct] *= TAPulse("TRIG", trigChan, pulseLength, 1.0, 0.0, 0.0)

def contains_measurement(entry):
    '''
    Determines if a LL entry contains a measurement
    '''
    if entry.label == "MEAS":
        return True
    elif hasattr(entry, 'pulses'):
        for p in entry.pulses.values():
            if p.label == "MEAS":
                return True
    return False

def add_slave_trigger(seqs, slaveChan):
    '''
    Add the slave trigger to each sequence.
    '''
    pulseLength = slaveChan.pulseParams['length'] * slaveChan.physChan.samplingRate
    for seq in seqs:
        # skip if the sequence already starts with a slave trig
        if hasattr(seq[0], 'qubits') and seq[0].qubits == slaveChan:
            continue
        seq.insert(0, TAPulse("TRIG", slaveChan, pulseLength, 1.0, 0.0, 0.0))
