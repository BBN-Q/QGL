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

from collections import namedtuple

class Pulse(namedtuple("Pulse", ["label", "channel", "length", "amp", "phase", "frequency",
                                 "frameChange", "shapeParams", "isTimeAmp",
                                 "isZero", "ignoredStrParams"])):
    __slots__ = ()

    def __new__(cls, label, channel, shapeParams, amp=1.0, phase=0, frameChange=0, ignoredStrParams=[]):
        if hasattr(channel, 'frequency'):
            frequency = channel.frequency
        else:
            frequency = 0
        requiredParams = ['length', 'shapeFun']
        for param in requiredParams:
            if param not in shapeParams.keys():
                raise NameError("shapeParams must include {0}".format(param))
        isTimeAmp = (shapeParams['shapeFun'] == PulseShapes.constant)
        isZero = (amp == 0)
        return super(cls, Pulse).__new__(cls, label, channel,
                                         shapeParams['length'], amp, phase,
                                         frequency, frameChange, shapeParams,
                                         isTimeAmp, isZero, ignoredStrParams)

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

    def hashshape(self):
        return hash(frozenset(self.shapeParams.items()))

    def __add__(self, other):
        if self.channel != other.channel:
            raise NameError(
                "Can only concatenate pulses acting on the same channel")
        return CompositePulse("", [self, other])

    # unary negation inverts the pulse amplitude and frame change
    def __neg__(self):
        return Pulse(self.label, self.channel, copy(self.shapeParams),
                     -self.amp, self.phase, -self.frameChange)

    def __mul__(self, other):
        """ Overload multiplication of Pulses as a "tensor" operator"""
        ptype = promote_type(self, other)
        return self.promote(ptype) * other.promote(ptype)

    def promote(self, ptype):
        # promote a Pulse to a PulseBlock
        return ptype(self)

    @property
    def shape(self):
        params = copy(self.shapeParams)
        params['samplingRate'] = self.channel.physChan.samplingRate
        params.pop('shapeFun')
        return self.shapeParams['shapeFun'](**params)


def TAPulse(label,
            channel,
            length,
            amp,
            phase=0,
            frameChange=0,
            ignoredStrParams=None):
    """
    Creates a time/amplitude pulse with the given pulse length and amplitude
    """
    params = {'length': length, 'shapeFun': PulseShapes.constant}
    if ignoredStrParams:
        if 'shapeFun' not in ignoredStrParams:
            ignoredStrParams.append('shapeFun')
    else:
        ignoredStrParams = ['shapeFun']
    return Pulse(label, channel, params, amp, phase, frameChange, ignoredStrParams)


class CompositePulse(namedtuple("CompositePulse", ["label", "pulses"])):
    '''
    A sequential series of pulses that reside within one time bin of a pulse block
    '''
    __slots__ = ()

    def __str__(self):
        if self.label != "":
            return '{0}({1})'.format(self.label, self.channel.label)
        else:
            return "+".join([str(p) for p in self.pulses])

    def _repr_pretty_(self, p, cycle):
        p.text(str(self))

    def __add__(self, other):
        if self.channel != other.channel:
            raise NameError(
                "Can only concatenate pulses acting on the same channel")
        if hasattr(other, 'pulses'):
            return CompositePulse("", self.pulses + other.pulses)
        else:
            return CompositePulse("", self.pulses + [other])

    def __mul__(self, other):
        ptype = promote_type(self, other)
        return self.promote(ptype) * other.promote(ptype)

    def promote(self, ptype):
        # promote a CompositePulse to a PulseBlock
        return ptype(self)

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
    The basic building block for pulse sequences. This is what we can concatenate
    together to make sequences. We overload the * operator so that we can combine
    pulse blocks on different channels.
    '''

    def __init__(self, *pulses):
        self.alignment = 'left'
        self.pulses = OrderedDict([(pulse.channel, pulse) for pulse in pulses])
        # The maximum length for any channel on this block
        self.length = max(p.length for p in self.pulses.values())
        self.label = None

    def __repr__(self):
        return "Pulses " + ";".join([
            str(pulse) for pulse in self.pulses.values()
        ]) + " alignment: {0}".format(self.alignment)

    def __str__(self):
        labelPart = "{0}: ".format(self.label) if self.label else ""
        return labelPart + "*".join(
            [str(pulse) for pulse in self.pulses.values()])

    def _repr_pretty_(self, p, cycle):
        labelPart = "{0}: ".format(self.label) if self.label else ""
        p.text(labelPart + "⊗ ".join([str(pulse)
                                      for pulse in self.pulses.values()]))

    # Overload the multiplication operator to combine pulse blocks
    def __mul__(self, rhs):
        ptype = promote_type(self, rhs)
        if ptype == CompoundGate:
            return rhs * self
        # otherwise, we are promoting to a PulseBlock
        rhs = rhs.promote(ptype)

        # copy PulseBlock so we don't modify other references
        result = copy(self)
        result.pulses = copy(self.pulses)
        for (k, v) in rhs.pulses.items():
            if k in result.pulses.keys():
                raise NameError(
                    "Attempted to multiply pulses acting on the same space")
            else:
                result.pulses[k] = v
        result.length = max(self.length, rhs.length)
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

    @property
    def channel(self):
        return self.pulses.keys()

    def promote(self, ptype):
        if ptype == CompoundGate:
            return CompoundGate(self)
        elif ptype == PulseBlock:
            # promoting to PulseBlock, so we can just return self
            return self
        else:
            ptype(self)

def align(pulseBlock, mode="center"):
    # make sure we have a PulseBlock
    pulseBlock = pulseBlock.promote(PulseBlock)
    pulseBlock.alignment = mode
    return pulseBlock

class CompoundGate(object):
    '''
    A wrapper around a python list to allow us to define '*' on lists.
    Used by multi-pulse structures like CNOT_CR so that we can retain the
    "boundaries" of the oepration.
    '''
    def __init__(self, seq, label=None):
        if isinstance(seq, list):
            self.seq = seq
        else:
            self.seq = [seq]
        self.label = label

    def __str__(self):
        txt = ', '.join(str(s) for s in self.seq)
        if self.label:
            return "{0}([{1}])".format(self.label, txt)
        else:
            return "CompoundGate([" + txt + "])"

    def _repr_pretty_(self, p, cycle):
        p.text(str(self))

    def __iter__(self):
        for n in range(len(self.seq)):
            yield self.seq[n]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return all(s1 == s2 for s1, s2 in zip(self.seq, other.seq))
        return False

    def __mul__(self, other):
        if isinstance(other, CompoundGate):
            other_seq = other.seq
        else:
            other_seq = [other]
        new_seq = []
        # zip together the sequences
        for s1, s2 in zip(self.seq, other_seq):
            new_seq.append(s1 * s2)

        # append any dangling section
        if len(self.seq) < len(other_seq):
            new_seq += other_seq[len(self.seq):]
        elif len(self.seq) > len(other_seq):
            new_seq += self.seq[len(other_seq):]
        return CompoundGate(new_seq)

    def __add__(self, other):
        if isinstance(other, CompoundGate):
            return CompoundGate(self.seq + other.seq)
        else:
            new_seq = copy(self.seq)
            new_seq.append(other)
            return CompoundGate(new_seq)

    @property
    def channel(self):
        # FIXME, should probably iterate over self.seq to build the effective
        # channel set
        return self.seq[0].channel

    def promote(self, ptype):
        # CompoundGates cannot be promoted
        return self

def promote_type(lhs, rhs):
    '''
    Returns the appropriate type for the '*' operator given operands 'lhs' and
    'rhs'.
    '''
    if isinstance(lhs, CompoundGate) or isinstance(rhs, CompoundGate):
        return CompoundGate
    else:
        return PulseBlock
