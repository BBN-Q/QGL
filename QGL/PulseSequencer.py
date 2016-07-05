# coding: utf-8
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

from __future__ import unicode_literals
from copy import copy
import json
import numpy as np
import operator
from collections import OrderedDict
from functools import reduce
from builtins import str

from . import ChannelLibrary, PulseShapes


class Pulse(object):
    '''
    A single channel pulse object
        label - name of the pulse
        channel - logical channel the pulse acts upon
        shape - numpy array pulse shape
        frameChange - accumulated phase from the pulse
    '''
    def __init__(self,
                 label,
                 channel,
                 shapeParams,
                 amp=1.0,
                 phase=0,
                 frameChange=0,
                 ignoredStrParams=[]):
        self.label = label
        if isinstance(channel, (list, tuple)):
            # with more than one qubit, need to look up the channel
            self.channel = ChannelLibrary.QubitFactory(reduce(operator.add, [c.label for c in channel]))
        else:
            self.channel = channel
        self.shapeParams = shapeParams
        self.phase = phase
        self.amp = amp
        self.frameChange = frameChange
        self.isTimeAmp = False
        if hasattr(self.channel, 'frequency'):
            self.frequency = self.channel.frequency
        else:
            self.frequency = 0
        self.ignoredStrParams = ignoredStrParams
        requiredParams = ['length', 'shapeFun']
        for param in requiredParams:
            if param not in shapeParams.keys():
                raise NameError("shapeParams must include {0}".format(param))
        if (self.shapeParams['shapeFun'] == PulseShapes.constant or
            self.shapeParams['shapeFun'] == PulseShapes.square):
            self.isTimeAmp = True

    def __repr__(self):
        return "Pulse({0}, {1}, {2}, {3}, {4}, {5}, {6})".format(
            self.label,
            self.channel,
            self.shapeParams,
            self.amp,
            self.phase,
            self.frameChange,
            self.ignoredStrParams)

    def __str__(self):
        kwvals = []
        # object parameters outside of shapeParams
        for param in ["amp", "phase", "frameChange"]:
            if param not in self.ignoredStrParams:
                kwvals.append("{0}={1}".format(param, getattr(self, param)))
        # parameters inside shapeParams
        for n, v in self.shapeParams.items():
            if (n not in self.ignoredStrParams and
                n in self.channel.pulseParams and
                self.channel.pulseParams[n] != v):
                kwvals.append("{0}={1}".format(n, v))
        if kwvals:
            kwstr = ", " + ", ".join(kwvals)
        else:
            kwstr = ""
        return '{0}({1}{2})'.format(self.label, self.channel.label, kwstr)

    def _repr_pretty_(self, p, cycle):
        p.text(str(self))

    # adding pulses concatenates the pulse shapes
    def __add__(self, other):
        if self.channel != other.channel:
            raise NameError("Can only concatenate pulses acting on the same channel")
        return CompositePulse([self, other])

    # unary negation inverts the pulse amplitude and frame change
    def __neg__(self):
        return Pulse(self.label, self.channel, copy(self.shapeParams), -self.amp, self.phase, -self.frameChange)

    def __mul__(self, other):
        return self.promote()*other.promote()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self == other

    def hashshape(self):
        return hash(frozenset(self.shapeParams.items()))

    def promote(self):
        # promote a Pulse to a PulseBlock
        pb =  PulseBlock()
        pb.pulses[self.channel] = self
        return pb

    @property
    def length(self):
        return self.shapeParams['length']

    @length.setter
    def length(self, value):
        self.shapeParams['length'] = value
        return value

    @property
    def isZero(self):
        return self.amp == 0

    @property
    def shape(self):
        params = copy(self.shapeParams)
        params['samplingRate'] = self.channel.physChan.samplingRate
        params.pop('shapeFun')
        return self.shapeParams['shapeFun'](**params)

def TAPulse(label, channel, length, amp, phase=0, frameChange=0, ignoredStrParams=[]):
    '''
    Creates a time/amplitude pulse with the given pulse length and amplitude
    '''
    params = {'length': length, 'shapeFun': PulseShapes.constant}
    ignoredStrParams.append('shapeFun')
    p = Pulse(label, channel, params, amp, phase, frameChange, ignoredStrParams)
    return p

class CompositePulse(object):
    '''
    A sequential series of pulses that reside within one time bin of a pulse block
    '''
    def __init__(self, pulses, label=""):
        self.pulses = pulses
        self.label = label

    def __repr__(self):
        return "CompositePulse({0}, {1})".format(self.pulses, self.label)

    def __str__(self):
        if self.label != "":
            return self.label
        else:
            return "+".join([str(p) for p in self.pulses])

    def _repr_pretty_(self, p, cycle):
        p.text(str(self))

    def __add__(self, other):
        if self.channel != other.channel:
            raise NameError("Can only concatenate pulses acting on the same channel")
        if hasattr(other, 'pulses'):
            return CompositePulse(self.pulses + other.pulses)
        else:
            return CompositePulse(self.pulses + [other])

    def __mul__(self, other):
        return self.promote()*other.promote()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self == other

    def promote(self):
        # promote a CompositePulse to a PulseBlock
        pb =  PulseBlock()
        pb.pulses[self.channel] = self
        return pb

    @property
    def channel(self):
        # Assume that the first pulse in the composite contains the channel information
        return self.pulses[0].channel

    @property
    def length(self):
        return sum(p.length for p in self.pulses)

    @property
    def frameChange(self):
        return sum(p.frameChange for p in self.pulses)

    @property
    def isZero(self):
        return all(p.isZero for p in self.pulses)

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
        self.pulses = OrderedDict()
        self.label = None

    def __repr__(self):
        return "Pulses " + ";".join([str(pulse) for pulse in self.pulses.values()]) + " alignment: {0}".format(self.alignment)

    def __str__(self):
        labelPart = "{0}: ".format(self.label) if self.label else ""
        return labelPart + "*".join([str(pulse) for pulse in self.pulses.values()])

    def _repr_pretty_(self, p, cycle):
        labelPart = "{0}: ".format(self.label) if self.label else ""
        p.text(labelPart + "⊗ ".join([str(pulse) for pulse in self.pulses.values()]))

    #Overload the multiplication operator to combine pulse blocks
    def __mul__(self, rhs):
        # make sure RHS is a PulseBlock
        rhs = rhs.promote()
        # we need to go one layer deep in the copy so that we can manipulate self.pulses and self.channel w/o affecting the original object
        # should bundle this behavior into a __copy__ method
        result = copy(self)
        result.pulses = copy(self.pulses)

        for (k, v) in rhs.pulses.items():
            if k in result.pulses.keys():
                raise NameError("Attempted to multiply pulses acting on the same space")
            else:
                result.pulses[k] = v
        return result

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            # ignore label in equality testing
            mydict = self.__dict__.copy()
            otherdict = other.__dict__.copy()
            mydict.pop('label')
            otherdict.pop('label')
            return mydict == otherdict
        return False

    def __ne__(self, other):
        return not self == other

    #PulseBlocks don't need to be promoted, so just return self
    def promote(self):
        return self

    @property
    def channel(self):
        return self.pulses.keys()

    #The maximum number of points needed for any channel on this block
    @property
    def length(self):
        return max([p.length for p in self.pulses.values()])

def align(pulseBlock, mode="center"):
    # make sure we have a PulseBlock
    pulseBlock = pulseBlock.promote()
    pulseBlock.alignment = mode
    return pulseBlock
