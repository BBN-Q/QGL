'''
Channels is where we store information for mapping virtual (qubit) channel to
real channels.

Split from Channels.py on Jan 14, 2016.

Original Author: Colm Ryan

Copyright 2016 Raytheon BBN Technologies

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Include modification to yaml loader (MIT License) from
https://gist.github.com/joshbode/569627ced3076931b02f

Scientific notation fix for yaml from
https://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number
'''

import sys
import os
import re
import traceback
import datetime
import importlib
from pony.orm import *
import networkx as nx
import yaml

# FSEvents observer in watchdog cannot have multiple watchers of the same path
# use kqueue instead
if sys.platform == 'darwin':
    from watchdog.observers.kqueue import KqueueObserver as Observer
else:
    from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time

from . import config
from . import Channels
from . import PulseShapes

channelLib = None

def set_from_dict(obj, settings):
    for prop_name in obj.to_dict().keys():
        if prop_name in settings.keys():
            try:
                setattr(obj, prop_name, settings[prop_name])
            except Exception as e:
                print(f"{obj.label}: Error loading {prop_name} from config")

class ChannelLibrary(object):

    def __init__(self, library_file=None, blank=False, channelDict={}, **kwargs):
        """Create the channel library. We assume that the user wants the config file in the 
        usual locations specified in the config files."""

        # Load the basic config options from the yaml
        self.library_file = config.load_config(library_file)

        self.connectivityG = nx.DiGraph()
        
        self.channelDict = {c.label: c for  c in select(c for c in Channels.Channel)}

        # Update the global reference
        global channelLib
        channelLib = self

    #Dictionary methods
    def __getitem__(self, key):
        return self.channelDict[key]

    def __setitem__(self, key, value):
        self.channelDict[key] = value

    def __delitem__(self, key):
        del self.channelDict[key]

    def __contains__(self, key):
        return key in self.channelDict

    def keys(self):
        return self.channelDict.keys()

    def values(self):
        return self.channelDict.values()

    def build_connectivity_graph(self):
        # build connectivity graph
<<<<<<< HEAD
        self.connectivityG.clear()
        for chan in self.channelDict.values():
            if isinstance(chan,
                          Channels.Qubit) and chan not in self.connectivityG:
                self.connectivityG.add_node(chan)
        for chan in self.channelDict.values():
            if isinstance(chan, Channels.Edge):
                self.connectivityG.add_edge(chan.source, chan.target)
                self.connectivityG[chan.source][chan.target]['channel'] = chan

    def load_from_library(self, return_only=False):
        """Loads the YAML library, creates the QGL objects, and returns a list of the visited filenames
        for the filewatcher."""
        if not self.library_file:
            return
        try:
            with open(self.library_file, 'r') as FID:
                loader = config.Loader(FID)
                try:
                    tmpLib = loader.get_single_data()
                    filenames = loader.filenames
                finally:
                    loader.dispose()

            # Check to see if we have the mandatory sections
            for section in ['instruments', 'qubits']: #, 'filters']:
                if section not in tmpLib.keys():
                    raise ValueError("{} section not present in config file {}.".format(section, self.library_file))

            instr_dict   = tmpLib['instruments']
            qubit_dict   = tmpLib['qubits']
            # filter_dict  = tmpLib['filters']
            trigger_dict = tmpLib.get('markers', {}) # This section is optional
            edge_dict    = tmpLib.get('edges', {}) # This section is optional
            master_awgs  = []

            # Construct the channel library
            channel_dict = {}
            marker_lens = {}

            for name, instr in instr_dict.items():
                if "tx_channels" in instr.keys():
                    for chan_name, channel in instr["tx_channels"].items():
                        if channel is None:
                            params = {}
                        else:
                            params = {k: v for k,v in channel.items() if k in Channels.PhysicalQuadratureChannel.__atom_members__.keys()}
                        params["label"] = name + "-" + chan_name
                        params["instrument"] = name
                        params["translator"] = instr["type"] + "Pattern"
                        params["__module__"] = "QGL.Channels"
                        params["__class__"]  = "PhysicalQuadratureChannel"
                        channel_dict[params["label"]] = params
                if "rx_channels" in instr.keys():
                    for chan_name, channel in instr["rx_channels"].items():
                        if channel is None:
                            params = {}
                        else:
                            params = {k: v for k,v in channel.items() if k in Channels.PhysicalMarkerChannel.__atom_members__.keys()}
                        params["label"] = name + "-" + chan_name
                        params["instrument"] = name
                        params["translator"] = instr["type"] + "Pattern"
                        params["__module__"] = "QGL.Channels"
                        params["__class__"]  = "PhysicalMarkerChannel"
                        channel_dict[params["label"]] = params
                if "markers" in instr.keys():
                    for mark_name, marker in instr["markers"].items():
                        if marker is None:
                            params = {}
                        else:
                            params = {k: v for k,v in marker.items() if k in Channels.PhysicalMarkerChannel.__atom_members__.keys()}
                        params["label"] = name + "-" + mark_name
                        params["instrument"] = name
                        params["translator"] = instr["type"] + "Pattern"
                        params["__module__"] = "QGL.Channels"
                        params["__class__"]  = "PhysicalMarkerChannel"
                        if "length" in marker.keys():
                            marker_lens[params["label"]] = marker["length"]
                        channel_dict[params["label"]] = params
                if "master" in instr.keys() and instr["master"]:
                    if instr['type'] != 'TDM':
                        slave_chan = instr["slave_trig"] if "slave_trig" in instr.keys() else "slave"
                        master_awgs.append(name + "-" + slave_chan)
                    else:
                        master_awgs.append(name)
                # Eventually we should support multiple masters...
                if "slave_trig" in instr.keys():
                    params = {}
                    params["label"]        = "slave_trig"
                    params["phys_chan"]    = name + "-" + instr["slave_trig"]
                    if params["phys_chan"] in marker_lens.keys():
                        length = marker_lens[params["phys_chan"]]
                    else:
                        length = 1e-7
                    params["pulse_params"] = {"length": length, "shape_fun": "constant"}
                    params["__module__"]   = "QGL.Channels"
                    params["__class__"]    = "LogicalMarkerChannel"
                    channel_dict[params["label"]] = params

            # Establish the slave trigger, assuming for now that we have a single
            # APS master. This might change later.
            if len(master_awgs) > 1:
                raise ValueError("More than one AWG is marked as master.")
            # elif len(master_awgs) == 1  and instr_dict[master_awgs[0].split('-')[0]]['type'] != 'TDM':
            #     params = {}
            #     params["label"]       = "slave_trig"
            #     params["phys_chan"]    = master_awgs[0]
            #     if params["phys_chan"] in marker_lens.keys():
            #         length = marker_lens[params["phys_chan"]]
            #     else:
            #         length = 1e-7
            #     params["pulse_params"] = {"length": length, "shape_fun": "constant"}
            #     params["__module__"]  = "QGL.Channels"
            #     params["__class__"]   = "LogicalMarkerChannel"
            #     channel_dict[params["label"]] = params

            # for name, filt in filter_dict.items():
            #     if "StreamSelector" in filt["type"]:
            #         params = {k: v for k,v in filt.items() if k in Channels.ReceiverChannel.__atom_members__.keys()}
            #         params["label"]      = "RecvChan-" + name # instr_dict[filt["instrument"]]["name"] + "-" + name
            #         params["channel"]    = str(params["channel"]) # Convert to a string
            #         params["instrument"] = filt["source"]
            #         params["__module__"] = "QGL.Channels"
            #         params["__class__"]  = "ReceiverChannel"
            #         if "source" not in filt.keys():
            #             raise ValueError("No instrument (source) specified for Stream Selector")
            #         if filt["source"] not in instr_dict.keys() and filt["source"] not in channel_dict.keys():
            #             raise ValueError("Stream Selector source {} not found among list of instruments.".format(filt["source"]))
            #         params["instrument"] = filt["source"]

            #         channel_dict[params["label"]] = params

            for name, qubit in qubit_dict.items():
                # Create a stream selector
                rcv_inst, rcv_chan, rcv_stream = qubit["measure"]["receiver"].split()
                rcv_params = {}
                rcv_params["label"]      = "RecvChan-" + name + "-SS"
                rcv_params["channel"]    = str(rcv_chan)
                rcv_params["instrument"] = rcv_inst
                rcv_params["__module__"] = "QGL.Channels"
                rcv_params["__class__"]  = "ReceiverChannel"
                channel_dict[rcv_params["label"]] = rcv_params

                # Create the Qubits
                if len(qubit["control"]["AWG"].split()) != 2:
                    print("Control AWG specification for {} ({}) must have a device, channel".format(name, qubit["control"]["AWG"]))
                    raise ValueError("Control AWG specification for {} ({}) must have a device, channel".format(name, qubit["control"]["AWG"]))
                ctrl_instr, ctrl_chan = qubit["control"]["AWG"].split()
                params = {k: v for k,v in qubit["control"].items() if k in Channels.Qubit.__atom_members__.keys()}
                params["label"]      = name
                params["phys_chan"]   = ctrl_instr + "-" + ctrl_chan 
                params["__module__"] = "QGL.Channels"
                params["__class__"]  = "Qubit"
                channel_dict[params["label"]] = params
                if 'generator' in qubit["control"].keys():
                    channel_dict[params["phys_chan"]]["generator"] = qubit["control"]["generator"]

                # Create the measurements
                if len(qubit["measure"]["AWG"].split()) != 2:
                    print("Measurement AWG specification for {} ({}) must have a device, channel".format(name, qubit["measure"]["AWG"]))
                    raise ValueError("Measurement AWG specification for {} ({}) must have a device, channel".format(name, qubit["measure"]["AWG"]))
                meas_instr, meas_chan = qubit["measure"]["AWG"].split()
                params = {k: v for k,v in qubit["measure"].items() if k in Channels.Measurement.__atom_members__.keys()}
                params["label"]        = "M-{}".format(name)
                # parse the digitizer trigger from the marker dictionary, if available. If not, expected in the form Instr Ch
                dig_trig = trigger_dict.get(qubit["measure"]["trigger"], qubit["measure"]["trigger"]) if trigger_dict else qubit["measure"]["trigger"]
                params["trig_chan"]     = "digTrig-" + dig_trig
                params["phys_chan"]     = meas_instr + "-" + meas_chan
                params["meas_type"]     = "autodyne"
                params["receiver_chan"] = rcv_params["label"]
                params["__module__"]   = "QGL.Channels"
                params["__class__"]    = "Measurement"
                channel_dict[params["label"]] = params
                if 'generator' in qubit["measure"].keys():
                    channel_dict[params["phys_chan"]]["generator"] = qubit["measure"]["generator"]

                # Create the receiver channels
                if "receiver" in qubit["measure"].keys():
                    if len(qubit["measure"]["receiver"].split()) != 3:
                        print("Receiver specification for {} ({}) must have an instrument name, physical channel, and stream".format(name, qubit["measure"]["receiver"]))
                        raise ValueError("Receiver specification for {} ({}) must have a stream selector".format(name, qubit["measure"]["receiver"]))
                    phys_instr, phys_marker = dig_trig.split()
                    params = {}
                    params["label"]        = "digTrig-" + dig_trig
                    params["phys_chan"]     = phys_instr + "-" + phys_marker
                    if params["phys_chan"] in marker_lens.keys():
                        length = marker_lens[params["phys_chan"]]
                    else:
                        length = 3e-7
                    params["pulse_params"]  = {"length": length, "shape_fun": "constant"}
                    params["__module__"]   = "QGL.Channels"
                    params["__class__"]    = "LogicalMarkerChannel"
                    # Don't duplicate triggers to the same digitizer
                    if params["label"] not in channel_dict.keys():
                        channel_dict[params["label"]] = params
=======
        for chan in select(q for q in Channels.Qubit if q not in self.connectivityG):
            self.connectivityG.add_node(chan)
        for chan in select(e for e in Channels.Edge):
            self.connectivityG.add_edge(chan.source, chan.target)
            self.connectivityG[chan.source][chan.target]['channel'] = chan
>>>>>>> Ditch atom, move to Pony.orm for all channel library objects.


def QubitFactory(label, **kwargs):
    ''' Return a saved qubit channel or create a new one. '''
    thing = select(el for el in Channels.Qubit if el.label==label).first()
    if thing:
        return thing
    else:
        return Channels.Qubit(label=label, **kwargs)
    
def MeasFactory(label, **kwargs):
    ''' Return a saved measurement channel or create a new one. '''
    thing = select(el for el in Channels.Measurement if el.label==label).first()
    if thing:
        return thing
    else:
        return Channels.Measurement(label=label, **kwargs)

def MarkerFactory(label, **kwargs):
    ''' Return a saved Marker channel or create a new one. '''
    thing = select(el for el in Channels.LogicalMarkerChannel if el.label==label).first()
    if thing:
        return thing
    else:
        return Channels.LogicalMarkerChannel(label=label, **kwargs)

def EdgeFactory(source, target):
    if channelLib.connectivityG.has_edge(source, target):
        return channelLib.connectivityG[source][target]['channel']
    elif channelLib.connectivityG.has_edge(target, source):
        return channelLib.connectivityG[target][source]['channel']
    else:
        raise ValueError('Edge {0} not found in connectivity graph'.format((
            source, target)))

