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
from math import pi
import hashlib, collections
import pickle
from copy import copy

from .PulseSequencer import Pulse, TAPulse, PulseBlock, CompositePulse
from .PulsePrimitives import BLANK
from . import ControlFlow
from . import BlockLabel
import QGL.drivers

def hash_pulse(shape):
    return hashlib.sha1(shape.tostring()).hexdigest()


TAZKey = hash_pulse(np.zeros(1, dtype=np.complex))


def delay(sequences, delay):
    '''
    Delays a sequence by the given amount.
    '''
    if delay <= 0:  # no need to inject zero delays
        return
    for seq in sequences:
        # loop through and look for WAIT instructions
        # use while loop because len(seq) will change as we inject delays
        ct = 0
        while ct < len(seq) - 1:
            if seq[ct] == ControlFlow.Wait() or seq[ct] == ControlFlow.Sync():
                seq.insert(ct + 1, TAPulse("Id", seq[ct + 1].channel, delay,
                                           0))
            ct += 1


def normalize_delays(delays):
    '''
    Normalizes a dictionary of channel delays. Postive delays shift right, negative delays shift left.
    Since we cannot delay by a negative amount in hardware, shift all delays until they are positive.
    Takes in a dict of channel:delay pairs and returns a normalized copy of the same.
    '''
    out = dict(delays)  # copy before modifying
    if not out or len(out) == 0:
        # Typically an error (empty sequence)
        import logging
        logging.error("normalize_delays() had no delays?")
        return out
    min_delay = min(delays.values())
    if min_delay < 0:
        for chan in delays.keys():
            out[chan] += -min_delay
    return out


def correct_mixers(wfLib, T):
    for k, v in wfLib.items():
        # To get the broadcast to work in numpy, need to do the multiplication one row at a time
        iqWF = np.vstack((np.real(v), np.imag(v)))
        wfLib[k] = T[0, :].dot(iqWF) + 1j * T[1, :].dot(iqWF)


def add_gate_pulses(seqs):
    '''
    add gating pulses to Qubit pulses
    '''
    for seq in seqs:
        for ct in range(len(seq)):
            if isinstance(seq[ct], PulseBlock):
                pb = None
                for chan, pulse in seq[ct].pulses.items():
                    if has_gate(chan) and not pulse.isZero and not (
                            chan.gateChan in seq[ct].pulses.keys()):
                        if pb:
                            pb *= BLANK(chan, pulse.length)
                        else:
                            pb = BLANK(chan, pulse.length)
                if pb:
                    seq[ct] *= pb
            elif hasattr(seq[ct], 'channel'):
                chan = seq[ct].channel
                if has_gate(chan) and not seq[ct].isZero:
                    seq[ct] *= BLANK(chan, seq[ct].length)


def has_gate(channel):
    return hasattr(channel, 'gateChan') and channel.gateChan

def update_pulse_length(pulse, new_length):
    """Return new Pulse with modified length"""
    assert new_length >= 0
    #copy shape parameter dictionary to avoid updating other copies
    new_params = copy(pulse.shapeParams)
    new_params["length"] = new_length
    return pulse._replace(shapeParams=new_params)

def apply_gating_constraints(chan, linkList):
    # get channel parameters in samples
    if not hasattr(chan, 'gateBuffer'):
        raise AttributeError("{0} does not have gateBuffer".format(chan.label))

    if not hasattr(chan, 'gateMinWidth'):
        raise AttributeError("{0} does not have gateMinWidth".format(
            chan.label))

    # get channel parameters
    gateBuffer = chan.gateBuffer
    gateMinWidth = chan.gateMinWidth

    #Initialize list of sequences to return
    gateSeqs = []

    for miniLL in linkList:
        gateSeq = []
        # first pass consolidates entries
        previousEntry = None
        for ct,entry in enumerate(miniLL):
            if isinstance(entry, (ControlFlow.ControlInstruction,
                                  BlockLabel.BlockLabel)):
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
                previousEntry = update_pulse_length(previousEntry, previousEntry.length + entry.length)
            else:
                gateSeq.append(previousEntry)
                previousEntry = entry

        # push on the last entry if necessary
        if previousEntry:
            gateSeq.append(previousEntry)

        # second pass expands non-zeros by gateBuffer
        for ct in range(len(gateSeq)):
            if isNonZeroWaveform(gateSeq[ct]):
                gateSeq[ct] = update_pulse_length(gateSeq[ct], gateSeq[ct].length + gateBuffer)

                # contract the next pulse by the same amount
                if ct + 1 < len(gateSeq) - 1 and isinstance(gateSeq[ct + 1], Pulse):
                    gateSeq[ct+1] = update_pulse_length(gateSeq[ct+1], gateSeq[ct+1].length - gateBuffer)

        # third pass ensures gateMinWidth
        ct = 0
        while ct + 2 < len(gateSeq):
            # look for pulse, delay, pulse pattern and ensure delay is long enough
            if [isNonZeroWaveform(x) for x in gateSeq[ct:ct+3]] == [True, False, True] and \
                gateSeq[ct+1].length < gateMinWidth and \
                [isinstance(x, Pulse) for x in gateSeq[ct:ct+3]] == [True, True, True]:
                gateSeq[ct] = update_pulse_length(gateSeq[ct], gateSeq[ct + 1].length + gateSeq[ct + 2].length)
                del gateSeq[ct + 1:ct + 3]
            else:
                ct += 1
        gateSeqs.append(gateSeq)

    return gateSeqs


def isNonZeroWaveform(entry):
    return isinstance(entry, Pulse) and not entry.isZero


def add_digitizer_trigger(seqs):
    '''
    Add a digitizer trigger to a logical LL (pulse blocks).
    '''
    # Attach a trigger to any pulse block containing a measurement. Each trigger is specific to each measurement
    for seq in seqs:
        for ct in range(len(seq)):
            if not contains_measurement(seq[ct]):
                continue
            #find corresponding digitizer trigger
            chanlist = list(flatten([seq[ct].channel]))
            for chan in chanlist:
                if hasattr(chan, 'trigChan'):
                    trigChan = chan.trigChan
                    if not (hasattr(seq[ct], 'pulses') and
                            trigChan in seq[ct].pulses.keys()):
                        seq[ct] *= TAPulse("TRIG", trigChan,
                                           trigChan.pulseParams['length'], 1.0,
                                           0.0, 0.0)


def contains_measurement(entry):
    """
    Determines if a LL entry contains a measurement
    """
    if entry.label == "MEAS":
        return True
    elif isinstance(entry, PulseBlock):
        for p in entry.pulses.values():
            if p.label == "MEAS":
                return True
    return False


def add_slave_trigger(seqs, slaveChan):
    '''
    Add the slave trigger to each sequence.
    '''
    for seq in seqs:
        # Attach a TRIG immediately after a WAIT.
        ct = 0
        while ct < len(seq) - 1:
            if isinstance(seq[ct], ControlFlow.Wait):
                try:
                    seq[ct + 1] *= TAPulse("TRIG", slaveChan,
                                           slaveChan.pulseParams['length'],
                                           1.0, 0.0, 0.0)
                except TypeError:
                    seq.insert(ct + 1, TAPulse("TRIG", slaveChan,
                                               slaveChan.pulseParams['length'],
                                               1.0, 0.0, 0.0))
                ct += 2  # can skip over what we just modified
            else:
                ct += 1


def propagate_frame_changes(seq, wf_type):
    '''
    Propagates all frame changes through sequence
    '''
    frame = 0
    for entry in seq:
        if not isinstance(entry, wf_type):
            continue
        entry.phase = np.mod(frame + entry.phase, 2 * pi)
        frame += entry.frameChange + (-2 * np.pi * entry.frequency *
                                      entry.length
                                      )  #minus from negative frequency qubits
    return seq


def quantize_phase(seqs, precision, wf_type):
    '''
    Quantizes waveform phases with given precision (in radians).
    '''
    for entry in flatten(seqs):
        if not isinstance(entry, wf_type):
            continue
        phase = np.mod(entry.phase, 2 * np.pi)
        entry.phase = precision * round(phase / precision)
    return seqs


def convert_lengths_to_samples(instructions, samplingRate, quantization=1, wf_type=None):
    for entry in flatten(instructions):
        if isinstance(entry, wf_type):
            entry.length = int(round(entry.length * samplingRate))
            # TODO: warn when truncating?
            entry.length -= entry.length % quantization
    return instructions

def convert_length_to_samples(wf_length, sampling_rate, quantization=1):
    num_samples = int(round(wf_length * sampling_rate))
    num_samples -= num_samples % quantization
    return num_samples

# from Stack Overflow: http://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists-in-python/2158532#2158532
def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, (str, Pulse, CompositePulse)) :
            for sub in flatten(el):
                yield sub
        else:
            yield el


def update_wf_library(pulses, path):
    """
    Update the waveform library in-place.

    Parameters
    ------------
    pulses : iterable of pulse object to update
        e.g. [X90(q1), X(q1), Y90(q1), Y(q1), X90(q2), X(q2), Y90(q2), Y(q2), ZX90_CR(q1, q2)]
    path : path to base name of files to update e.g. /path/to/GST/GST will update files such as
        /path/to/GST/GST-APSII1.h5 and /path/to/GST/GST-APSII2.h5
    """
    #Look through the pulses and figure out what pulses are associated with which APS
    awg_pulses = collections.defaultdict(dict)
    translators = {}

    def flatten_pulses():
        for p in flatten(pulses):
            if isinstance(p, CompositePulse):
                for sub_p in p.pulses:
                    yield sub_p
            else:
                yield p

    pulse_list = list(flatten_pulses())
    for ct, pulse in enumerate(pulse_list):
        awg = pulse.channel.physChan.AWG
        if awg not in translators:
            translators[awg] = getattr(QGL.drivers,
                                       pulse.channel.physChan.translator)
        if pulse.label not in awg_pulses[awg]:
            awg_pulses[awg][pulse.label] = pulse_list[ct]

    for awg, ps in awg_pulses.items():
        #load the offset dictionary for this AWG
        try:
            with open(path + "-" + awg + ".offsets", "rb") as FID:
                offsets = pickle.load(FID)
        except IOError:
            print("Offset file not found for {}, skipping pulses {}".format(
                awg, [str(p) for p in ps.values()]))
            continue
        print("Updating pulses for {}".format(awg))
        translators[awg].update_wf_library(path + "-" + awg + ".h5", ps,
                                           offsets)
