'''

Quantum Gate Language Module

'''

from copy import copy
import json
import numpy as np
import matplotlib.pyplot as plt
import Compiler

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
        self.shape = shape.astype(np.complex) # for now, do this since we are gettern objects from PatternGen rather than lists
        self.phase = phase
        self.frameChange = frameChange

    def __repr__(self):
        return str(self)

    def __str__(self):
        if isinstance(self.qubits, tuple):
            return '{0}({1})'.format(self.label, ','.join([qubit.name for qubit in self.qubits]))
        else:
            return '{0}({1})'.format(self.label, self.qubits.name)

    # adding pulses concatenates the pulse shapes
    def __add__(self, other):
        newLabel = self.label+"+"+other.label
        if self.qubits != other.qubits:
            raise NameError("Can only concatenate pulses acting on the same channel")
        return Pulse(newLabel, self.qubits, np.append(self.shape, other.shape), self.frameChange + other.frameChange)

    # unary negation inverts the pulse shape
    # TODO: does the frame change need to be updated??
    def __neg__(self):
        return Pulse(self.label, self.qubits, -self.shape, self.frameChange)

    def __mul__(self, other):
        return self.promote()*other.promote()

    def promote(self):
        # promote a Pulse to a PulseBlock
        pb =  PulseBlock()
        pb.pulses = {self.qubits: self}
        return pb

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
        # we to go one layer deep in the copy so that we can manipulate self.pulses and self.channels w/o affecting the original object
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
        return max([len(p.shape) for p in self.pulses.values()])

def align(pulseBlock, mode="center"):
    # make sure we have a PulseBlock
    pulseBlock = pulseBlock.promote()
    pulseBlock.alignment = mode
    return pulseBlock

AWGFreq = 1e9

def show(seq):
    #compile
    linkList, wfLib = Compiler.compile_sequence(seq)

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
