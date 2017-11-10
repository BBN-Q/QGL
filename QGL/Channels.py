'''
Channels is where we store information for mapping virtual (qubit) channel to real channels.

Created on Jan 19, 2012

Original Author: Colm Ryan

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

from . import PulseShapes
import numpy as np

from math import tan, cos, pi
# Python 2/3 compatibility: use the py3 meaning of 'str'
from builtins import str

from atom.api import Atom, Str, Float, Instance, \
    Dict, Enum, Bool, Typed, Int

from copy import deepcopy


class Channel(Atom):
    '''
    Every channel has a label and some printers.
    '''
    label = Str()
    enabled = Bool(True)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "{0}('{1}')".format(self.__class__.__name__, self.label)

    def json_encode(self):
        jsonDict = self.__getstate__()

        #Turn objects into labels
        for member in ["phys_chan", "gate_chan", "trig_chan", "receiver_chan", "source", "target"]:
            if member in jsonDict and not isinstance(jsonDict[member], str):
                obj = jsonDict.pop(member)
                if obj:
                    jsonDict[member] = obj.label

        #We want the name of shape functions
        if "pulse_params" in jsonDict:
            pulse_params = deepcopy(jsonDict.pop("pulse_params"))
            if "shape_fun" in pulse_params:
                pulse_params["shape_fun"] = pulse_params["shape_fun"].__name__
            jsonDict["pulse_params"] = pulse_params

        return jsonDict


class PhysicalChannel(Channel):
    '''
    The main class for actual AWG channels.
    '''
    instrument = Str() # i.e. the AWG or receiver
    translator = Str()
    generator = Str()
    sampling_rate = Float(default=1.2e9)
    delay = Float()


class LogicalChannel(Channel):
    '''
    The main class from which we will generate sequences.
    At some point it needs to be assigned to a physical channel.
    '''
    #During initilization we may just have a string reference to the channel
    phys_chan = Instance((str, PhysicalChannel))

    def __init__(self, **kwargs):
        super(LogicalChannel, self).__init__(**kwargs)
        if self.phys_chan is None:
            self.phys_chan = PhysicalChannel(label=kwargs['label'] + '-phys')


class PhysicalMarkerChannel(PhysicalChannel):
    '''
    A digital output channel on an AWG.
    '''
    gate_buffer = Float(0.0).tag(
        desc="How much extra time should be added onto the beginning of a gating pulse")
    gate_min_width = Float(0.0).tag(desc="The minimum marker pulse width")


class PhysicalQuadratureChannel(PhysicalChannel):
    '''
    Something used to implement a standard qubit channel with two analog channels and a microwave gating channel.
    '''
    I_channel = Str()
    Q_channel = Str()
    #During initilization we may just have a string reference to the channel
    amp_factor = Float(1.0)
    phase_skew = Float(0.0)

    @property
    def correctionT(self):
        return np.array(
            [[self.amp_factor, self.amp_factor * tan(self.phase_skew * pi / 180)],
             [0, 1 / cos(self.phase_skew * pi / 180)]])

class ReceiverChannel(PhysicalChannel):
    '''
    A trigger input on a receiver.
    '''
    channel = Str()

class LogicalMarkerChannel(LogicalChannel):
    '''
    A class for digital channels for gating sources or triggering other things.
    '''
    pulse_params = Dict(default={'shape_fun': PulseShapes.constant,
                                'length': 10e-9})


class Qubit(LogicalChannel):
    '''
    The main class for generating qubit pulses.  Effectively a logical "QuadratureChannel".
    '''
    pulse_params = Dict(default={'length': 20e-9,
                                'piAmp': 1.0,
                                'pi2Amp': 0.5,
                                'shape_fun': PulseShapes.gaussian,
                                'cutoff': 2,
                                'drag_scaling': 0,
                                'sigma': 5e-9})
    gate_chan = Instance((str, LogicalMarkerChannel))
    frequency = Float(0.0).tag(
        desc='modulation frequency of the channel (can be positive or negative)')

    def __init__(self, **kwargs):
        super(Qubit, self).__init__(**kwargs)


class Measurement(LogicalChannel):
    '''
    A class for measurement channels.
    Measurements are special because they can be different types:
    autodyne which needs an IQ pair or hetero/homodyne which needs just a marker channel.
    '''
    meas_type = Enum('autodyne', 'homodyne').tag(
        desc='Type of measurement (autodyne, homodyne)')

    autodyne_freq = Float(0.0).tag(
        desc='use to bake the modulation into the pulse, so that it has constant phase')
    frequency = Float(0.0).tag(
        desc='use frequency to asssociate modulation with the channel')
    pulse_params = Dict(default={'length': 100e-9,
                                'amp': 1.0,
                                'shape_fun': PulseShapes.tanh,
                                'cutoff': 2,
                                'sigma': 1e-9})
    gate_chan = Instance((str, LogicalMarkerChannel))
    trig_chan = Instance((str, LogicalMarkerChannel))
    receiver_chan = Instance((str, ReceiverChannel))

    def __init__(self, **kwargs):
        super(Measurement, self).__init__(**kwargs)
        if self.trig_chan is None:
            self.trig_chan = LogicalMarkerChannel(label='digitizerTrig')


class Edge(LogicalChannel):
    '''
    Defines an arc/directed edge between qubit vertices. If a device supports bi-directional
    connectivity, that is represented with two independent Edges.

    An Edge is also effectively an abstract channel, so it carries the same properties as a
    Qubit channel.
    '''
    # allow string in source and target so that we can store a label or an object
    source = Instance((str, Qubit))
    target = Instance((str, Qubit))
    pulse_params = Dict(default={'length': 20e-9,
                                'amp': 1.0,
                                'phase': 0.0,
                                'shape_fun': PulseShapes.gaussian,
                                'cutoff': 2,
                                'drag_scaling': 0,
                                'sigma': 5e-9,
                                'riseFall': 20e-9})
    gate_chan = Instance((str, LogicalMarkerChannel))
    frequency = Float(0.0).tag(
        desc='modulation frequency of the channel (can be positive or negative)')

    def __init__(self, **kwargs):
        super(Edge, self).__init__(**kwargs)

    def isforward(self, source, target):
        ''' Test whether (source, target) matches the directionality of the edge. '''
        nodes = (self.source, self.target)
        if (source not in nodes) or (target not in nodes):
            raise ValueError('One of {0} is not a node in the edge'.format((
                source, target)))
        if (self.source, self.target) == (source, target):
            return True
        else:
            return False


NewLogicalChannelList = [Qubit, Edge, LogicalMarkerChannel, Measurement]
NewPhysicalChannelList = [
    PhysicalMarkerChannel,
    PhysicalQuadratureChannel,
    ReceiverChannel
]
