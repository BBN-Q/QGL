'''

Quantum Gate Language Module

'''

from copy import copy

import numpy as np
import matplotlib.pyplot as plt

# from Channels import Qubit
# from PulsePrimitives import *
# import PulseSequencePlotter

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
    # normalize sequence to PulseBlocks
    seq = [p.promote() for p in seq]
    # find set of used channels
    qubits = set([])
    for step in seq:
        qubits |= set(step.pulses.keys())

    # initialize empty arrays for each qubit
    concatShapes = {q: np.array([], dtype=np.complex128) for q in qubits}

    # now loop through steps and push on the pulse shape, or an Id operation if none specified
    for step in seq:
        stepLength = step.maxPts
        for q in qubits:
            if q in step.pulses.keys():
                concatShapes[q] = np.append(concatShapes[q], step.pulses[q].shape)
            else:
                concatShapes[q] = np.append(concatShapes[q], np.zeros(stepLength))
    
    # plot
    for (ct,q) in enumerate(qubits):
        plt.subplot(len(qubits),1,ct+1)
        waveformToPlot = concatShapes[q]
        p = plt.plot(np.linspace(0,len(waveformToPlot)/AWGFreq,len(waveformToPlot)), np.real(waveformToPlot), 'r')
        p = plt.plot(np.linspace(0,len(waveformToPlot)/AWGFreq,len(waveformToPlot)), np.imag(waveformToPlot), 'b')
        plt.ylim((-1.05,1.05))
    plt.show(p)


if __name__ == '__main__':
    pass
