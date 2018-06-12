'''
Channels is where we store information for mapping virtual (qubit) channel to
real channels.

Split from Channels.py on Jan 14, 2016.
Moved to pony ORM from atom June 1, 2018

Original Author: Colm Ryan
Modified By: Graham Rowlands

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
import datetime
import traceback
import datetime
import importlib
import inspect
from pony.orm import *
import networkx as nx

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

def copy_objs(chans, srcs, new_channel_db):
    new_chans       = []
    new_srcs        = []
    old_to_new_chan = {}
    old_to_new_src  = {}

    for chan in chans:
        c = copy_entity(chan, new_channel_db)
        new_chans.append(c)
        old_to_new_chan[chan] = c

    for src in srcs:
        c = copy_entity(src, new_channel_db)
        new_srcs.append(c)
        old_to_new_src[src] = c

    # Fix links... pony updates the relationships symmetriacally so we get some for free
    for thing in new_chans + new_srcs:
        for attr in thing._attrs_:
            if attr:
                if isinstance(getattr(thing, attr.name), Channels.Channel):
                    if getattr(thing, attr.name) in old_to_new_chan.keys():
                        setattr(thing, attr.name, old_to_new_chan[getattr(thing, attr.name)])
                elif isinstance(getattr(thing, attr.name), Channels.MicrowaveSource):
                    if getattr(thing, attr.name) in old_to_new_src.keys():
                        setattr(thing, attr.name, old_to_new_src[getattr(thing, attr.name)])

    return new_chans, new_srcs

def copy_entity(obj, new_channel_db):
    """Copy a pony entity instance"""
    kwargs = {a.name: getattr(obj, a.name) for a in obj.__class__._attrs_ if a.name not in ["id", "classtype", "pulse_params"]}
    # if "pulse_params" in kwargs.keys():
    #     kwargs["pulse_params"] = dict(kwargs["pulse_params"])
    kwargs["channel_db"] = new_channel_db
    return obj.__class__(**kwargs)

# def copy_entity(obj, new_channel_db):
#     """Copy a pony entity instance"""
#     attr_names = [a.name for a in obj.__class__._attrs_]
#     skip = ["id", "classtype", "channel_db"]
#     print(obj)
#     kwargs = {"channel_db": new_channel_db}
#     for a in attr_names: #["label", "source_type"]:
#         print("\t", a)
#         if a in dir(obj):
#             val = getattr(obj, a)
#             if a == "pulse_params":
#                 val = dict(val)
#             if val is not None and a not in skip:
#                 print("\t\t", getattr(obj, a))
#                 kwargs[a] = val
#     return obj.__class__(**kwargs)

class ChannelLibrary(object):

    def __init__(self, channel_db_name=None, database_file=None, channelDict={}, **kwargs):
        """Create the channel library."""

        global channelLib
        if channelLib is not None:
            channelLib.db.disconnect()

        config.load_db()
        if database_file:
            self.database_file = database_file
        elif config.db_file:
            self.database_file = config.db_file
        else:
            self.database_file = ":memory:"

        self.db = Database()
        Channels.define_entities(self.db)
        self.db.bind('sqlite', filename=self.database_file, create_db=True)
        self.db.generate_mapping(create_tables=True)

        # Dirty trick: push the correct entity defs to the calling context
        for var in ["Measurement","Qubit","Edge"]:
            inspect.stack()[1][0].f_globals[var] = getattr(Channels, var)

        self.connectivityG = nx.DiGraph()
        
        # This is still somewhere legacy QGL behavior. Massage db into dict for lookup.
        self.channelDict = {}
        self.channels = []
        self.sources = []
        self.channelDatabase = Channels.ChannelDatabase(label="__temp__", time=datetime.datetime.now())
        self.channel_db_name = channel_db_name if channel_db_name else "temp"

        config.load_config()

        # self.load_most_recent()
        # config.load_config()

        # Update the global reference
        channelLib = self

    def get_current_channels(self):
        return list(select(c for c in Channels.Channel if c.channel_db is None)) + list(select(c for c in Channels.MicrowaveSource if c.channel_db is None))

    def update_channelDict(self):
        self.channelDict = {c.label: c for c in self.get_current_channels()}

    def list(self):
        select((c.label, c.time, c.id) for c in Channels.ChannelDatabase).show()

    def load_by_id(self, id_num):
        obj = select(c for c in Channels.ChannelDatabase if c.id==id_num).first()
        self.load(obj)

    def clear(self):
        select(c for c in Channels.Channel if c.channel_db == self.channelDatabase).delete(bulk=True)
        select(c for c in Channels.MicrowaveSource if c.channel_db == self.channelDatabase).delete(bulk=True)
        self.channelDatabase.time  = datetime.datetime.now()

    def load(self, obj): #, delete=True):
        self.clear()

        chans = list(obj.channels)
        srcs  = list(obj.sources)

        # self.channelDatabase = Channels.ChannelDatabase(label="__temp__", time=datetime.datetime.now())
        new_chans, new_srcs = copy_objs(chans, srcs, self.channelDatabase)

        self.channels = new_chans
        self.sources = new_srcs
        self.channel_db_name = obj.label
        # self.channelDatabase = None

    # def load_most_recent(self, name=None):
    #     if name is None:
    #         name = self.channel_db_name
    #     mrcd = Channels.ChannelDatabase.select(lambda d: d.label==name).order_by(desc(Channels.ChannelDatabase.time)).first()
    #     if mrcd:
    #         self.load(mrcd)

    # def new(self, name):
    #     # self.channelDatabase.delete()
    #     self.clear()
    #     commit()

    #     # self.channelDatabase = Channels.ChannelDatabase(label="__temp__", time=datetime.datetime.now())

    #     # self.channelDatabase = None
    #     self.channel_db_name = name
    #     self.channels = []
    #     self.sources = []

    def save(self):
        self.save_as(self.channel_db_name)

    def save_as(self, name):
        # self.channelDatabase.label = name
        # self.channelDatabase.time  = datetime.datetime.now()
        # cd = self.channelDatabase
        # self.channelDatabase = None
        # commit()
        # self.load(cd, delete=False)
        
        # Get channels that are part of the currently active db
        # chans = list(select(c for c in Channels.Channel if c.channel_db is None))
        # srcs  = list(select(c for c in Channels.MicrowaveSource if c.channel_db is None))
        chans = list(self.channelDatabase.channels)
        srcs  = list(self.channelDatabase.sources)
        cd    = Channels.ChannelDatabase(label=name, time=datetime.datetime.now(), channels=chans, sources=srcs)
        new_chans, new_srcs = copy_objs(chans, srcs, cd)

        # self.channels = new_chans
        # self.sources = new_srcs
        # self.channelDatabase = None
        commit()
        # self.channel_db_name = name

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

# Convenience functions for generating and linking channels
# TODO: move these to a shim layer shared by Auspex/QGL

class APS2(object):
    def __init__(self, label, address=None, delay=0.0):
        self.chan12 = Channels.PhysicalQuadratureChannel(label=f"{label}-12", instrument=label, translator="APS2Pattern", channel_db=channelLib.channelDatabase)
        self.m1     = Channels.PhysicalMarkerChannel(label=f"{label}-12m1", instrument=label, translator="APS2Pattern", channel_db=channelLib.channelDatabase)
        self.m2     = Channels.PhysicalMarkerChannel(label=f"{label}-12m2", instrument=label, translator="APS2Pattern", channel_db=channelLib.channelDatabase)
        self.m3     = Channels.PhysicalMarkerChannel(label=f"{label}-12m3", instrument=label, translator="APS2Pattern", channel_db=channelLib.channelDatabase)
        self.m4     = Channels.PhysicalMarkerChannel(label=f"{label}-12m4", instrument=label, translator="APS2Pattern", channel_db=channelLib.channelDatabase)
        
        self.trigger_interval = None
        self.trigger_source   = "External"
        self.address          = address
        self.delay            = delay
        self.master           = False

class X6(object):
    def __init__(self, label, address=None):
        self.chan1 = Channels.ReceiverChannel(label=f"RecvChan-{label}-1", channel_db=channelLib.channelDatabase)
        self.chan2 = Channels.ReceiverChannel(label=f"RecvChan-{label}-2", channel_db=channelLib.channelDatabase)
        
        self.address          = address
        self.reference        = "external"
        self.nbr_segments     = 1
        self.nbr_round_robins = 100
        self.acquire_mode     = "digitizer"

def new_qubit(label):
    return Channels.Qubit(label=label, channel_db=channelLib.channelDatabase)

def new_source(label, source_type, address, power=-30.0):
    return Channels.MicrowaveSource(label=label, source_type=source_type, address=address, power=power, channel_db=channelLib.channelDatabase)

def set_control(qubit, awg, generator=None):
    qubit.phys_chan = awg.chan12
    if generator:
        qubit.phys_chan.generator = generator
    
def set_measure(qubit, awg, dig, generator=None, dig_channel=1, trig_channel=1, gate=False, gate_channel=2, trigger_length=1e-7):
    meas = Channels.Measurement(label=f"M-{qubit.label}", channel_db=channelLib.channelDatabase)
    meas.phys_chan     = awg.chan12
    
    meas.trig_chan              = Channels.LogicalMarkerChannel(label=f"digTrig-{qubit.label}", channel_db=channelLib.channelDatabase)
    meas.trig_chan.phys_chan    = getattr(awg, f"m{trig_channel}")
    meas.trig_chan.pulse_params = {"length": trigger_length, "shape_fun": "constant"}
    meas.receiver_chan          = getattr(dig, f"chan{dig_channel}")

    if generator:
        meas.phys_chan.generator = generator

    if gate:
        meas.gate_chan           = Channels.LogicalMarkerChannel(label=f"M-{qubit.label}-gate", channel_db=channelLib.channelDatabase)
        meas.gate_chan.phys_chan = getattr(awg, f"m{gate_channel}")
        
def set_master(awg, trig_channel=2, pulse_length=1e-7):
    st = Channels.LogicalMarkerChannel(label="slave_trig", channel_db=channelLib.channelDatabase)
    st.phys_chan = getattr(awg, f"m{trig_channel}")
    st.pulse_params = {"length": pulse_length, "shape_fun": "constant"}
    awg.master = True
    awg.trigger_source = "Internal"

def QubitFactory(label, **kwargs):
    ''' Return a saved qubit channel or create a new one. '''
    # TODO: this will just get the first entry in the whole damned DB!
    # thing = select(el for el in Channels.Qubit if el.label==label).first()
    thing = {c.label: c for c in channelLib.get_current_channels()}[label]
    if thing:
        return thing
    else:
        return Channels.Qubit(label=label, **kwargs)
    
def MeasFactory(label, **kwargs):
    ''' Return a saved measurement channel or create a new one. '''
    thing = {c.label: c for c in channelLib.get_current_channels()}[label]
    if thing:
        return thing
    else:
        return Channels.Measurement(label=label, **kwargs)

def MarkerFactory(label, **kwargs):
    ''' Return a saved Marker channel or create a new one. '''
    thing = {c.label: c for c in channelLib.get_current_channels()}[label]
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

