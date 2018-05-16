'''
Various tools to interface with pyGSTi for running GST experiments.

Created on May 16, 2018

Original Author: Guilhem Ribeill

Copyright 2018 Raytheon BBN Technologies

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
from .PulsePrimitives import *
from pygsti.objects import GateString
from itertools import chain

#Default mapping from pyGSTi naming convention to QGL gates.
gst_gate_map = {"Gx": X,
                "Gy": Y,
                "Gi": Id}

def gst_map_1Q(gst_list, qubit, qgl_map=gst_gate_map, append_meas=True):
    """
    Helper function that takes an arbitrarily nested list of pygsti gatestrings
    and converts them into QGL sequences, keeping the same nesting of lists.

    Inputs:
        gst_list: GateString to convert, or possibly nested list of pyGSTi GateStrings.
        qubit: QGL qubit to apply the sequence to
        qgl_map: Dictionary that maps between pyGSTi "Gx" string to QGL pulse
        append_meas: Append a measurement to each sequence.
    Returns:
        QGL sequences, preserving the input list nesting (as a generator)
    """
    if isinstace(gst_list, GateString):
        gst_list = [gst_list]
    for item in gst_list:
        if isinstance(item, GateString):
            mapped = map(lambda x: qgl_map[x](qubit), item.tup)
            if append_meas:
                yield list(chain(mapped, [MEAS(qubit)]))
            else:
                yield list(mapped)
        elif isinstance(item, list):
            yield list(gst_map(item, qubit, qgl_map=qgl_map, append_meas=append_meas))

def create_gst_sequence_from_pygsti(gst_list, qubit, gate_map=gst_gate_map):
    """ Returns list of QGL sequences from a pyGSTi GateString list. See gst_map_1Q.
    """
    return list(gst_map_1Q(gst_list, qubit, qgl_map=gate_map, append_meas=True))
