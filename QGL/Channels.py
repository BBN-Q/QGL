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
from bbndb import define_entities as def_ent
from bbndb import qgl

# Get these in global scope for module imports
Channel                   = None
PhysicalChannel           = None
LogicalChannel            = None
PhysicalQuadratureChannel = None
PhysicalMarkerChannel     = None
LogicalMarkerChannel      = None
ReceiverChannel           = None
Measurement               = None
Qubit                     = None
Edge                      = None
MicrowaveSource           = None
ChannelDatabase           = None
Digitizer                 = None
AWG                       = None

# The main definitions have been moved to bbndb
# to keep the schema consistent across Auspex and QGL
# We retain this local framework in order to maintain
# the paths for QGL primitives

def define_entities(db):
    def_ent(db)

    globals()["Channel"]                   = qgl.Channel
    globals()["PhysicalChannel"]           = qgl.PhysicalChannel
    globals()["LogicalChannel"]            = qgl.LogicalChannel
    globals()["PhysicalQuadratureChannel"] = qgl.PhysicalQuadratureChannel
    globals()["PhysicalMarkerChannel"]     = qgl.PhysicalMarkerChannel
    globals()["LogicalMarkerChannel"]      = qgl.LogicalMarkerChannel
    globals()["ReceiverChannel"]           = qgl.ReceiverChannel
    globals()["Measurement"]               = qgl.Measurement
    globals()["Qubit"]                     = qgl.Qubit
    globals()["Edge"]                      = qgl.Edge
    globals()["MicrowaveSource"]           = qgl.MicrowaveSource
    globals()["ChannelDatabase"]           = qgl.ChannelDatabase
    globals()["Digitizer"]                 = qgl.Digitizer
    globals()["AWG"]                       = qgl.AWG