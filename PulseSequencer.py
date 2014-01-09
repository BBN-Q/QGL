'''
Quantum Gate Language Module

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

from copy import copy
import json
import numpy as np
import matplotlib.pyplot as plt

class Pulse(object):
    '''
    A single channel pulse object
        label - name of the pulse
        qubits - array of qubit/channel objects the pulse acts upon
        shape - numpy array pulse shape
        frameChange - accumulated phase from the pulse
    '''
    def __init__(self, label, qubits, shape, phase, frameChange):
        self.label = label
        self.qubits = qubits
        self.shape = shape.astype(np.complex)
        self.phase = phase
        self.frameChange = frameChange
        self.isTimeAmp = False
        self.TALength = 0

    def __repr__(self):
        return str(self)

    def __str__(self):
        if isinstance(self.qubits, tuple):
            return '{0}({1})'.format(self.label, ','.join([qubit.label for qubit in self.qubits]))
        else:
            return '{0}({1})'.format(self.label, self.qubits.label)

    # adding pulses concatenates the pulse shapes
    def __add__(self, other):
        if self.qubits != other.qubits:
            raise NameError("Can only concatenate pulses acting on the same channel")
        return CompositePulse([self, other])

    # unary negation inverts the pulse shape and frame change
    def __neg__(self):
        return Pulse(self.label, self.qubits, -self.shape, self.phase, -self.frameChange)

    def __mul__(self, other):
        return self.promote()*other.promote()

    def promote(self):
        # promote a Pulse to a PulseBlock
        pb =  PulseBlock()
        pb.pulses = {self.qubits: self}
        return pb

    @property
    def length(self):
        if self.isTimeAmp:
            return self.TALength
        else:
            return len(self.shape)

def TAPulse(label, qubits, length, amp, phase, frameChange):
    '''
    Creates a time/amplitude pulse with the given pulse length and amplitude
    '''
    p = Pulse(label, qubits, np.array([amp], np.complex), phase, frameChange)
    p.isTimeAmp = True
    p.TALength = length
    return p

class CompositePulse(object):
    '''
    A sequential series of pulses that reside within one time bin of a pulse block
    '''
    def __init__(self, pulses, label=""):
        self.pulses = pulses
        self.label = label

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.label != "":
            return self.label
        else:
            return "+".join([str(p) for p in self.pulses])

    def __add__(self, other):
        if self.qubits != other.qubits:
            raise NameError("Can only concatenate pulses acting on the same channel")
        if hasattr(other, 'pulses'):
            return CompositePulse(self.pulses + other.pulses)
        else:
            return CompositePulse(self.pulses + [other])

    def __mul__(self, other):
        return self.promote()*other.promote()

    def promote(self):
        # promote a CompositePulse to a PulseBlock
        pb =  PulseBlock()
        pb.pulses = {self.qubits: self}
        return pb

    @property
    def qubits(self):
        # Assume that the first pulse in the composite contains the qubit information
        return self.pulses[0].qubits

    @property
    def length(self):
        return sum([p.length for p in self.pulses])

    @property
    def frameChange(self):
        return sum([p.frameChange for p in self.pulses])


class PulseBlock(object):
    '''
    The basic building block for pulse sequences. This is what we can concatenate together to make sequences.
    We overload the * operator so that we can combine pulse blocks on different channels.
    We overload the + operator to concatenate pulses on the same channel.
    '''
    def __init__(self):
        #Set some default values
        #How multiple channels are aligned.
        self.alignment = 'left'
        self.pulses = {}

    def __repr__(self):
        return u"\u2297 ".join([str(pulse) for pulse in self.pulses.values()]).encode('utf-8')

    def __str__(self):
        return "Pulses " + ";".join([str(pulse) for pulse in self.pulses.values()]) + " alignment: {0}".format(self.alignment).encode('utf-8')

    #Overload the multiplication operator to combine pulse blocks
    def __mul__(self, rhs):
        # make sure RHS is a PulseBlock
        rhs = rhs.promote()
        # we need to go one layer deep in the copy so that we can manipulate self.pulses and self.channels w/o affecting the original object
        # should bundle this behavior into a __copy__ method
        result = copy(self)
        result.pulses = copy(self.pulses)
        
        for (k, v) in rhs.pulses.items():
            if k in result.pulses.keys():
                raise NameError("Attempted to multiply pulses acting on the same space")
            else:
                result.pulses[k] = v
        return result

    #PulseBlocks don't need to be promoted, so just return self
    def promote(self):
        return self

    #A list of the channels used in this block
    @property
    def channelNames(self):
        return [channel.name for channel in self.pulses.keys()]

    #The maximum number of points needed for any channel on this block
    @property
    def maxPts(self):
        return max([p.length for p in self.pulses.values()])

def align(pulseBlock, mode="center"):
    # make sure we have a PulseBlock
    pulseBlock = pulseBlock.promote()
    pulseBlock.alignment = mode
    return pulseBlock

AWGFreq = 1e9

def show(seq):
    from Compiler import compile_sequence, normalize #import here to avoid circular imports 
    #compile
    linkList, wfLib = compile_sequence(normalize(seq))

    # build a concatenated waveform for each channel
    channels = linkList.keys()
    concatShapes = {q: np.array([0], dtype=np.complex128) for q in channels}
    for q in channels:
        for entry in linkList[q]:
            if entry.isTimeAmp:
                concatShapes[q] = np.append(concatShapes[q], wfLib[q][entry.key][0]*np.ones(entry.length*entry.repeat))
            else:
                concatShapes[q] = np.append(concatShapes[q], np.tile(wfLib[q][entry.key], (1, entry.repeat)) )
    # add an extra zero to make things look more normal
    for q in channels:
        concatShapes[q] = np.append(concatShapes[q], 0)
    
    # plot
    for (ct,chan) in enumerate(channels):
        plt.subplot(len(channels),1,ct+1)
        waveformToPlot = concatShapes[chan]
        xpts = np.linspace(0,len(waveformToPlot)/AWGFreq/1e-6,len(waveformToPlot))
        p = plt.plot(xpts, np.real(waveformToPlot), 'r')
        p = plt.plot(xpts, np.imag(waveformToPlot), 'b')
        plt.ylim((-1.05,1.05))
        plt.title(repr(chan))
    plt.show(p)


if __name__ == '__main__':
    pass
