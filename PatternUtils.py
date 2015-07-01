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
import ControlFlow, BlockLabel, Compiler
from math import pi
import hashlib, collections

def hash_pulse(shape):
    return hashlib.sha1(shape.tostring()).hexdigest()

TAZKey = hash_pulse(np.zeros(1, dtype=np.complex))

def delay(sequences, delay):
    '''
    Delays a sequence by the given amount.
    '''
    if delay <= 0: # no need to inject zero delays
        return
    for seq in sequences:
        # loop through and look for WAIT instructions
        # use while loop because len(seq) will change as we inject delays
        ct = 0
        while ct < len(seq)-1:
            if seq[ct] == ControlFlow.Wait() or seq[ct] == ControlFlow.Sync():
                seq.insert(ct+1, TAPulse("Id", seq[ct+1].qubits, delay, 0))
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

def correct_mixers(wfLib, T):
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
    if not hasattr(chan,'gateBuffer'):
        raise AttributeError("{0} does not have gateBuffer".format(chan.label))

    if not hasattr(chan,'gateMinWidth'):
        raise AttributeError("{0} does not have gateMinWidth".format(chan.label))
      
    # get channel parameters
    gateBuffer = chan.gateBuffer
    gateMinWidth = chan.gateMinWidth

    #Initialize list of sequences to return
    gateSeqs = []

    for miniLL in linkList:
        gateSeq = []
        # first pass consolidates entries
        previousEntry = None
        for entry in miniLL:
            if isinstance(entry, (ControlFlow.ControlInstruction, BlockLabel.BlockLabel)):
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
                if ct + 1 < len(gateSeq) - 1 and isinstance(gateSeq[ct+1], Pulse):
                    gateSeq[ct+1].length -= gateBuffer #TODO: what if this becomes negative?

        # third pass ensures gateMinWidth
        ct = 0
        while ct+2 < len(gateSeq):
            # look for pulse, delay, pulse pattern and ensure delay is long enough
            if [isNonZeroWaveform(x) for x in gateSeq[ct:ct+3]] == [True, False, True] and \
                gateSeq[ct+1].length < gateMinWidth and \
                [isinstance(x, Pulse) for x in gateSeq[ct:ct+3]] == [True, True, True]:
                gateSeq[ct].length += gateSeq[ct+1].length + gateSeq[ct+2].length
                del gateSeq[ct+1:ct+3]
            else:
                ct += 1

        gateSeqs.append(gateSeq)

    return gateSeqs

def isNonZeroWaveform(entry):
    return isinstance(entry, Pulse) and not entry.isZero

def add_digitizer_trigger(seqs, trigChan):
    '''
    Add the digitizer trigger to a logical LL (pulse blocks).
    '''
    # Attach a trigger to any pulse block containing a measurement
    for seq in seqs:
        for ct in range(len(seq)):
            if contains_measurement(seq[ct]) and not (hasattr(seq[ct], 'pulses') and trigChan in seq[ct].pulses.keys()):
                seq[ct] *= TAPulse("TRIG", trigChan, trigChan.pulseParams['length'], 1.0, 0.0, 0.0)

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
    for seq in seqs:
        # skip if the sequence already starts with a slave trig
        if hasattr(seq[0], 'qubits') and seq[0].qubits == slaveChan:
            continue
        seq.insert(0, TAPulse("TRIG", slaveChan, slaveChan.pulseParams['length'], 1.0, 0.0, 0.0))

def propagate_frame_changes(seq):
    '''
    Propagates all frame changes through sequence
    '''
    frame = 0
    for entry in seq:
        if not isinstance(entry, Compiler.Waveform):
            continue
        entry.phase = np.mod(frame + entry.phase, 2*pi)
        frame += entry.frameChange + (2*np.pi * entry.frequency * entry.length)
    return seq

def quantize_phase(seqs, precision):
    '''
    Quantizes waveform phases with given precision (in radians).
    '''
    for entry in flatten(seqs):
        if not isinstance(entry, Compiler.Waveform):
            continue
        phase = np.mod(entry.phase, 2*np.pi)
        entry.phase = precision * round(phase / precision)
    return seqs

def convert_lengths_to_samples(instructions, samplingRate):
    for entry in flatten(instructions):
        if isinstance(entry, Compiler.Waveform):
            entry.length = int(round(entry.length * samplingRate))
    return instructions

# from Stack Overflow: http://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists-in-python/2158532#2158532
def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el
