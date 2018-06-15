'''
Channels is where we store information for mapping virtual (qubit) channel to real channels.

Created on Jan 19, 2012

Original Author: Colm Ryan
Modified By: Graham Rowlands

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

from . import config
from . import PulseShapes
import numpy as np
from math import tan, cos, pi
from pony.orm import *
from copy import deepcopy
import datetime

# Get these in global scope for module imports
Channel = None
PhysicalChannel = None
LogicalChannel = None
PhysicalQuadratureChannel = None
PhysicalMarkerChannel = None
LogicalMarkerChannel = None
ReceiverChannel = None
Measurement = None
Qubit = None
Edge = None
MicrowaveSource = None
ChannelDatabase = None
Digitizer = None
AWG = None

def define_entities(db):

    class ChannelDatabase(db.Entity):
        label      = Required(str)
        channels   = Set("Channel", cascade_delete=True)
        sources    = Set("MicrowaveSource", cascade_delete=True)
        awgs       = Set("AWG", cascade_delete=True)
        digitizers = Set("Digitizer", cascade_delete=True)
        time       = Optional(datetime.datetime)

    class MicrowaveSource(db.Entity):
        label           = Required(str)
        source_type     = Required(str)
        address         = Optional(str)
        power           = Optional(float)
        logical_channel = Optional("PhysicalChannel")
        channel_db      = Optional("ChannelDatabase")

    class Digitizer(db.Entity):
        label           = Required(str)
        address         = Optional(str)
        # stream_types    = Set(str)
        channels        = Set("ReceiverChannel")
        trigger_source  = Required(str, default="External", py_check=lambda x: x in ['External', 'Internal'])
        channel_db      = Optional("ChannelDatabase")
        # nbr_segments     = Required(int, default=1) # This should be automatic
        # nbr_round_robins = Required(int, default=100) # This should be automatic
        # acquire_mode     = Required(str,default="digitizer", py_check=lambda x: x in ['digitizer', 'averager'])
        
        def get_chan(self, name):
            return self.channels.select(lambda x: x.label.endswith(name)).first()

    class AWG(db.Entity):
        label            = Required(str)
        address          = Optional(str)
        channels         = Set("PhysicalChannel", reverse="awg")
        trigger_interval = Required(float, default=100e-6)
        trigger_source   = Required(str, default="External", py_check=lambda x: x in ['External', 'Internal'])
        delay            = Required(float, default=0.0)
        master           = Required(bool, default=False)
        channel_db       = Optional("ChannelDatabase")

        def get_chan(self, name):
            return self.channels.select(lambda x: x.label.endswith(name)).first()

    class Channel(db.Entity):
        '''
        Every channel has a label and some printers.
        '''
        label      = Required(str)
        channel_db = Optional("ChannelDatabase")

        def __repr__(self):
            return str(self)
        def __str__(self):
            return "{0}('{1}')".format(self.__class__.__name__, self.label)

    class PhysicalChannel(Channel):
        '''
        The main class for actual AWG channels.
        '''
        instrument      = Optional(str) # i.e. the AWG or receiver
        translator      = Optional(str)
        generator       = Optional(MicrowaveSource, reverse="logical_channel")
        sampling_rate   = Required(float, default=1.2e9)
        delay           = Required(float, default=0.0)

        # Required reverse connections
        logical_channel = Optional("LogicalChannel")
        # quad_channel_I  = Optional("PhysicalQuadratureChannel", reverse="I_channel")
        # quad_channel_Q  = Optional("PhysicalQuadratureChannel", reverse="Q_channel")
        # marker_channel  = Optional("PhysicalMarkerChannel")
        awg             = Optional(AWG)

    class LogicalChannel(Channel):
        '''
        The main class from which we will generate sequences.
        At some point it needs to be assigned to a physical channel.
            frequency: modulation frequency of the channel (can be positive or negative)
        '''
        #During initilization we may just have a string reference to the channel
        phys_chan     = Optional(PhysicalChannel)
        frequency     = Required(float, default=0.0)
        pulse_params  = Optional(Json, default={}) 
        gate_chan     = Optional("LogicalMarkerChannel")
        receiver_chan = Optional("ReceiverChannel")

    class PhysicalMarkerChannel(PhysicalChannel):
        '''
        A digital output channel on an AWG.
            gate_buffer: How much extra time should be added onto the beginning of a gating pulse
            gate_min_width: The minimum marker pulse width
        '''
        gate_buffer    = Required(float, default=0.0)
        gate_min_width = Required(float, default=0.0)
        # phys_channel   = Optional(PhysicalChannel)

    class PhysicalQuadratureChannel(PhysicalChannel):
        '''
        Something used to implement a standard qubit channel with two analog channels and a microwave gating channel.
        '''
        amp_factor           = Required(float, default=1.0)
        phase_skew           = Required(float, default=0.0)
        I_channel_offset     = Required(float, default=0.0)
        Q_channel_offset     = Required(float, default=0.0)
        I_channel_amp_factor = Required(float, default=1.0)
        Q_channel_amp_factor = Required(float, default=1.0)

    class ReceiverChannel(PhysicalChannel):
        '''
        A trigger input on a receiver.
        '''
        triggering_channel = Optional(LogicalChannel)
        channel            = Optional(int)
        stream_type        = Required(str, default="Raw", py_check=lambda x: x.lower() in["raw", "demodulated", "integrated", "averaged"])
        if_freq            = Required(float, default=0.0)
        kernel             = Optional(str)
        kernel_bias        = Required(float, default=0.0)
        threshold          = Required(float, default=0.0)
        threshold_invert   = Required(bool, default=False)
        digitizer          = Optional(Digitizer)

    def pulse_check(name):
        return name in ["constant", "gaussian", "drag", "gaussOn", "gaussOff", "dragGaussOn", "dragGaussOff",
                       "tanh", "exp_decay", "autodyne", "arb_axis_drag"]
        
    class LogicalMarkerChannel(LogicalChannel):
        '''
        A class for digital channels for gating sources or triggering other things.
        '''
        meas_channel = Optional(LogicalChannel)
        trig_channel = Optional("Measurement")

        def __init__(self, **kwargs):
            if "pulse_params" not in kwargs.keys():
                kwargs["pulse_params"] =  {'shape_fun': "constant",'length': 10e-9}
            super(LogicalMarkerChannel, self).__init__(**kwargs)

    class Qubit(LogicalChannel):
        '''
        The main class for generating qubit pulses.  Effectively a logical "QuadratureChannel".
            frequency: modulation frequency of the channel (can be positive or negative)
        '''
        edge_source = Optional("Edge", reverse="source")
        edge_target = Optional("Edge", reverse="target")

        def __init__(self, **kwargs):
            if "pulse_params" not in kwargs.keys():
                kwargs["pulse_params"] =  {'length': 20e-9,
                                'piAmp': 1.0,
                                'pi2Amp': 0.5,
                                'shape_fun': "gaussian",
                                'cutoff': 2,
                                'drag_scaling': 0,
                                'sigma': 5e-9}
            super(Qubit, self).__init__(**kwargs)

    class Measurement(LogicalChannel):
        '''
        A class for measurement channels.
        Measurements are special because they can be different types:
        autodyne which needs an IQ pair or hetero/homodyne which needs just a marker channel.
            meas_type: Type of measurement (autodyne, homodyne)
            autodyne_freq: use to bake the modulation into the pulse, so that it has constant phase
            frequency: use to asssociate modulation with the channel
        '''
        meas_type     = Required(str, default='autodyne', py_check=lambda x: x in ['autodyne', 'homodyne'])
        autodyne_freq = Required(float, default=0.0)
        trig_chan     = Optional(LogicalMarkerChannel)

        def __init__(self, **kwargs):
            if "pulse_params" not in kwargs.keys():
                kwargs["pulse_params"] =  {'length': 100e-9,
                                    'amp': 1.0,
                                    'shape_fun': "tanh",
                                    'cutoff': 2,
                                    'sigma': 1e-9}
            super(Measurement, self).__init__(**kwargs)

    class Edge(LogicalChannel):
        '''
        Defines an arc/directed edge between qubit vertices. If a device supports bi-directional
        connectivity, that is represented with two independent Edges.

        An Edge is also effectively an abstract channel, so it carries the same properties as a
        Qubit channel.
        '''
        # allow string in source and target so that we can store a label or an object
        source = Required(Qubit)
        target = Required(Qubit)
        
        def __init__(self, **kwargs):
            if "pulse_params" not in kwargs.keys():
                kwargs["pulse_params"] =   {'length': 20e-9,
                                            'amp': 1.0,
                                            'phase': 0.0,
                                            'shape_fun': "gaussian",
                                            'cutoff': 2,
                                            'drag_scaling': 0,
                                            'sigma': 5e-9,
                                            'riseFall': 20e-9}
            super(Edge, self).__init__(**kwargs)

    globals()["Channel"] = Channel
    globals()["PhysicalChannel"] = PhysicalChannel
    globals()["LogicalChannel"] = LogicalChannel
    globals()["PhysicalQuadratureChannel"] = PhysicalQuadratureChannel
    globals()["PhysicalMarkerChannel"] = PhysicalMarkerChannel
    globals()["LogicalMarkerChannel"] = LogicalMarkerChannel
    globals()["ReceiverChannel"] = ReceiverChannel
    globals()["Measurement"] = Measurement
    globals()["Qubit"] = Qubit
    globals()["Edge"] = Edge
    globals()["MicrowaveSource"] = MicrowaveSource
    globals()["ChannelDatabase"] = ChannelDatabase
    globals()["Digitizer"] = Digitizer
    globals()["AWG"] = AWG