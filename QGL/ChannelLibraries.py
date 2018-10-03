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
from functools import wraps
from pony.orm import *
import numpy as np
import networkx as nx

from . import config
from . import Channels
from . import PulseShapes
from .PulsePrimitives import clear_pulse_cache
import bbndb

channelLib = None
db = None

def localize_db_objects(f):
    """Since we can't mix db objects from separate sessions, re-fetch entities by their unique IDs"""
    @wraps(f)
    def wrapper(*args, **kwargs):
        global db
        with db_session:
            args_info    = [(type(arg), arg.id) if isinstance(arg, db.Entity) else (arg, None) for arg in args]
            kwargs_info  = {k: (type(v), v.id) if isinstance(v, db.Entity) else (v, None) for k, v in kwargs.items()}
        with db_session:
            new_args   = [c[i] if i else c for c, i in args_info]
            new_kwargs = {k: v[0][v[1]] if v[1] else v[0] for k,v in kwargs_info.items()}
            return f(*new_args, **new_kwargs)
    return wrapper

def set_from_dict(obj, settings):
    for prop_name in obj.to_dict().keys():
        if prop_name in settings.keys():
            try:
                setattr(obj, prop_name, settings[prop_name])
            except Exception as e:
                print(f"{obj.label}: Error loading {prop_name} from config")

@db_session
def copy_objs(*entities, new_channel_db):
    # Entities is a list of lists of entities of specific types
    new_entities    = []
    old_to_new      = {}
    links_to_change = {}

    for ent in entities:
        new_ents = []
        for obj in ent:
            c, links = copy_entity(obj, new_channel_db)
            new_ents.append(c)
            links_to_change[c] = links
            old_to_new[c.label] = c
        new_entities.append(new_ents)

    for chan, link_info in links_to_change.items():
        for attr_name, link_name in link_info.items():
            if isinstance(link_name, pony.orm.core.Multiset):
                new = [old_to_new[ln] for ln in link_name]
            else:
                new = old_to_new[link_name]
            setattr(chan, attr_name, new)

    return new_entities

@db_session
def copy_entity(obj, new_channel_db):
    """Copy a pony entity instance"""
    kwargs = {a.name: getattr(obj, a.name) for a in obj._attrs_ if a.name not in ["id", "classtype"]}

    # Extract any links to other entities
    links = {}
    for attr in obj._attrs_:
        if attr.name not in ["channel_db"]:
            obj_attr = getattr(obj, attr.name)
            if hasattr(obj_attr, "id"):
                kwargs.pop(attr.name)
                links[attr.name] = obj_attr.label

    kwargs["channel_db"] = new_channel_db
    return obj.__class__(**kwargs), links

class ChannelLibrary(object):

    def __init__(self, database_file=None, channelDict={}, **kwargs):
        """Create the channel library."""

        global channelLib, db
        if channelLib is not None:
            channelLib.db.disconnect()

        config.load_db()
        self.database_provider = "sqlite"
        if database_file:
            self.database_file = database_file
        elif config.db_file:
            self.database_file = config.db_file
        else:
            self.database_file = ":memory:"

        db = Database()
        Channels.define_entities(db, cache_callback=clear_pulse_cache)
        db.bind(self.database_provider, filename=self.database_file, create_db=True)
        db.generate_mapping(create_tables=True)
        bbndb.database = db

        # Dirty trick: push the correct entity defs to the calling context
        for var in ["Measurement","Qubit","Edge"]:
            inspect.stack()[1][0].f_globals[var] = getattr(Channels, var)

        self.connectivityG = nx.DiGraph()

        # This is still somewhere legacy QGL behavior. Massage db into dict for lookup.
        self.channelDict = {}

        # Check to see whether there is already a temp database
        with db_session:
            w_dbs = list(select(d for d in Channels.ChannelDatabase if d.label == "working"))
            if len(w_dbs) > 1:
                # self.clear(channel_db=cdb, create_new=False)
                raise Exception("More than one working database exists!")
            elif len(w_dbs) == 1:
                self.channelDatabase = w_dbs[0]
                commit()
                self.update_channelDict()
            elif len(w_dbs) == 0:
                self.channelDatabase = Channels.ChannelDatabase(label="working", time=datetime.datetime.now())
            # commit()

        config.load_config()

        # Update the global reference
        channelLib = self

    @db_session
    def get_current_channels(self):
        cdb = Channels.ChannelDatabase[self.channelDatabase.id] # Can't use external object
        return list(cdb.channels) + list(cdb.sources)

    @db_session
    def update_channelDict(self):
        commit()
        self.channelDict = {c.label: c for c in self.get_current_channels()}

    @db_session
    def ls(self):
        select((c.label, c.time, c.id) for c in Channels.ChannelDatabase).sort_by(1, 2).show()

    @db_session
    def ent_by_type_name(self, name, show=False):
        q = select(c for c in getattr(bbndb.qgl,name) if c.channel_db.label == "working")
        if show:
            select(c.label for c in getattr(bbndb.qgl,name) if c.channel_db.label == "working").sort_by(1).show()
        else:
            return {el.label: el for el in q}

    def receivers(self):
        return self.ent_by_type_name("Receiver")

    def transmitter(self):
        return self.ent_by_type_name("Transmitter")

    def qubits(self):
        return self.ent_by_type_name("Qubit")

    def meas(self):
        return self.ent_by_type_name("Measurement")

    def ls_receivers(self):
        return self.ent_by_type_name("Receiver", show=True)

    def ls_transmitters(self):
        return self.ent_by_type_name("Transmitter", show=True)

    def ls_qubits(self):
        return self.ent_by_type_name("Qubit", show=True)

    def ls_meas(self):
        return self.ent_by_type_name("Measurement", show=True)

    @db_session
    def load(self, name, index=1):
        """Load the latest instance for a particular name. Specifying index = 2 will select the second most recent instance """
        obj = list(select(c for c in Channels.ChannelDatabase if c.label==name).sort_by(desc(Channels.ChannelDatabase.time)))
        self.load_obj(obj[-index])

    @db_session
    def load_by_id(self, id_num):
        obj = select(c for c in Channels.ChannelDatabase if c.id==id_num).first()
        self.load_obj(obj)

    @db_session
    def clear(self, channel_db=None, create_new=True):
        # If no database is specified, clear self.database
        channel_db = channel_db if channel_db else self.channelDatabase
        # First clear items that don't have Sets of other items
        for ent in [Channels.MicrowaveSource, Channels.Channel, Channels.Transmitter, Channels.Receiver, Channels.Transceiver]:
            select(c for c in ent if c.channel_db == channel_db).delete(bulk=True)
            commit()
        # Now clear items that do potentially have sets of items (which should be deleted)
        for ent in [Channels.ChannelDatabase]:
            select(d for d in ent if d.label == "working").delete(bulk=True)
            commit()
        if create_new:
            self.channelDatabase = Channels.ChannelDatabase(label="working", time=datetime.datetime.now())
            commit()

    @db_session
    def load_obj(self, obj):
        commit()
        self.clear()
        chans, srcs, d2as, a2ds, trans = map(list, [obj.channels, obj.sources, obj.transmitters, obj.receivers, obj.transceivers])
        copy_objs(chans, srcs, d2as, a2ds, trans, new_channel_db=self.channelDatabase)
        commit()
        self.update_channelDict()

    @db_session
    def save_as(self, name):
        chans, srcs, d2as, a2dsm, trans = map(list, [self.channelDatabase.channels, self.channelDatabase.sources,
                                            self.channelDatabase.transmitters, self.channelDatabase.receivers,
                                            self.channelDatabase.transceivers])
        commit()
        cd = Channels.ChannelDatabase(label=name, time=datetime.datetime.now())
        new_chans, new_srcs, new_d2as, new_a2ds, new_trans = copy_objs(chans, srcs, d2as, a2ds, trans, new_channel_db=cd)
        cd.channels, cd.sources, cd.transmitters, cd.receivers, cd.transceivers = new_chans, new_srcs, new_d2as, new_a2ds, new_trans
        commit()

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

    @db_session
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

@db_session
def new_APS2(label, address):
    cdb    = Channels.ChannelDatabase[channelLib.channelDatabase.id] # Can't use external object
    chan12 = Channels.PhysicalQuadratureChannel(label=f"{label}-12", instrument=label, translator="APS2Pattern", channel_db=cdb)
    m1     = Channels.PhysicalMarkerChannel(label=f"{label}-12m1", instrument=label, translator="APS2Pattern", channel_db=cdb)
    m2     = Channels.PhysicalMarkerChannel(label=f"{label}-12m2", instrument=label, translator="APS2Pattern", channel_db=cdb)
    m3     = Channels.PhysicalMarkerChannel(label=f"{label}-12m3", instrument=label, translator="APS2Pattern", channel_db=cdb)
    m4     = Channels.PhysicalMarkerChannel(label=f"{label}-12m4", instrument=label, translator="APS2Pattern", channel_db=cdb)

    this_transmitter = Channels.Transmitter(label=label, model="APS2", address=address, channels=[chan12, m1, m2, m3, m4], channel_db=cdb)
    this_transmitter.trigger_source = "External"
    this_transmitter.address        = address

    commit()
    return this_transmitter

@db_session
def new_APS2_rack(label, num, start_address):
    cdb              = Channels.ChannelDatabase[channelLib.channelDatabase.id] # Can't use external object
    address_start    = ".".join(start_address.split(".")[:3])
    address_end      = int(start_address.split(".")[-1])
    transmitters     = [new_APS2(f"{label}_U{i}", f"{address_start}.{address_end+i}") for i in range(1,num+1)]
    this_transceiver = Channels.Transceiver(label=label, model="APS2Rack", transmitters=transmitters, channel_db=cdb)

    commit()
    return this_transceiver

@db_session
def new_X6(label, address, dsp_channel=0, record_length=1024):
    cdb   = Channels.ChannelDatabase[channelLib.channelDatabase.id] # Can't use external object
    chan1 = Channels.ReceiverChannel(label=f"RecvChan-{label}-1", channel=1, dsp_channel=dsp_channel, channel_db=cdb)
    chan2 = Channels.ReceiverChannel(label=f"RecvChan-{label}-2", channel=2, dsp_channel=dsp_channel, channel_db=cdb)

    this_receiver = Channels.Receiver(label=label, model="X6-1000M", address=address, channels=[chan1, chan2],
                                  record_length=record_length, channel_db=cdb)
    this_receiver.trigger_source = "External"
    this_receiver.stream_types   = "raw, demodulated, integrated"
    this_receiver.address        = address

    # Add a default kernel
    chan1.kernel = np.ones(record_length, dtype=np.complex).tobytes()
    chan2.kernel = np.ones(record_length, dtype=np.complex).tobytes()

    commit()
    return this_receiver

@db_session
def new_Alazar(label, address, record_length=1024):
    cdb   = Channels.ChannelDatabase[channelLib.channelDatabase.id] # Can't use external object
    chan1 = Channels.ReceiverChannel(label=f"RecvChan-{label}-1", channel=1, channel_db=cdb)
    chan2 = Channels.ReceiverChannel(label=f"RecvChan-{label}-2", channel=2, channel_db=cdb)

    this_receiver = Channels.Receiver(label=label, model="AlazarATS9870", address=address, channels=[chan1, chan2],
                                  record_length=record_length, channel_db=cdb)
    this_receiver.trigger_source = "External"
    this_receiver.stream_types   = "raw"
    this_receiver.address        = address

    commit()
    return this_receiver

@db_session
def new_qubit(label):
    cdb   = Channels.ChannelDatabase[channelLib.channelDatabase.id] # Can't use external object
    thing = Channels.Qubit(label=label, channel_db=cdb)
    return thing

@db_session
def new_source(label, model, address, power=-30.0, frequency=5.0e9):
    cdb   = Channels.ChannelDatabase[channelLib.channelDatabase.id] # Can't use external object
    thing = Channels.MicrowaveSource(label=label, model=model, address=address, power=power, frequency=frequency, channel_db=cdb)
    return thing

@localize_db_objects
def set_control(qubit, transmitter, generator=None):
    quads   = [c for c in transmitter.channels if isinstance(c, Channels.PhysicalQuadratureChannel)]
    markers = [c for c in transmitter.channels if isinstance(c, Channels.PhysicalMarkerChannel)]

    if isinstance(transmitter, Channels.Transmitter) and len(quads) > 1:
        raise ValueError("In set_control the Transmitter must have a single quadrature channel or a specific channel must be passed instead")
    elif isinstance(transmitter, Channels.Transmitter) and len(quads) == 1:
        phys_chan = quads[0]
    elif isinstance(transmitter, Channels.PhysicalQuadratureChannel):
        phys_chan = transmitter
    else:
        raise ValueError("In set_control the Transmitter must have a single quadrature channel or a specific channel must be passed instead")

    qubit.phys_chan = phys_chan
    if generator:
        qubit.phys_chan.generator = generator

@localize_db_objects
def set_measure(qubit, transmitter, receivers, generator=None, receivers_channel=1, trig_channel=None, gate=False, gate_channel=None, trigger_length=1e-7):
    cdb     = Channels.ChannelDatabase[channelLib.channelDatabase.id] # Can't use external object
    quads   = [c for c in transmitter.channels if isinstance(c, Channels.PhysicalQuadratureChannel)]
    markers = [c for c in transmitter.channels if isinstance(c, Channels.PhysicalMarkerChannel)]

    if isinstance(transmitter, Channels.Transmitter) and len(quads) > 1:
        raise ValueError("In set_measure the Transmitter must have a single quadrature channel or a specific channel must be passed instead")
    elif isinstance(transmitter, Channels.Transmitter) and len(quads) == 1:
        phys_chan = quads[0]
    elif isinstance(transmitter, Channels.PhysicalQuadratureChannel):
        phys_chan = transmitter
    else:
        raise ValueError("In set_measure the Transmitter must have a single quadrature channel or a specific channel must be passed instead")

    meas = Channels.Measurement(label=f"M-{qubit.label}", channel_db=cdb)
    meas.phys_chan = phys_chan
    if generator:
        meas.phys_chan.generator = generator

    phys_trig_channel = trig_channel if trig_channel else transmitter.get_chan("12m1")

    trig_chan              = Channels.LogicalMarkerChannel(label=f"receiversTrig-{qubit.label}", channel_db=cdb)
    trig_chan.phys_chan    = phys_trig_channel
    trig_chan.pulse_params = {"length": trigger_length, "shape_fun": "constant"}
    meas.trig_chan         = trig_chan

    if isinstance(receivers, Channels.Receiver) and len(receivers.channels) > 1:
        raise ValueError("In set_measure the Receiver must have a single receiver channel or a specific channel must be passed instead")
    elif isinstance(receivers, Channels.Receiver) and len(receivers.channels) == 1:
        rcv_chan = receivers.channels[0]
    elif isinstance(receivers, Channels.ReceiverChannel):
        rcv_chan = receivers
    else:
        raise ValueError("In set_measure the Transmitter must have a single quadrature channel or a specific channel must be passed instead")

    meas.receiver_chan = rcv_chan

    if gate:
        phys_gate_channel   = gate_channel if gate_channel else transmitter.get_chan("12m2")
        gate_chan           = Channels.LogicalMarkerChannel(label=f"M-{qubit.label}-gate", channel_db=cdb)
        gate_chan.phys_chan = phys_gate_channel
        meas.gate_chan      = gate_chan

@localize_db_objects
def set_master(transmitter, trig_channel, pulse_length=1e-7):
    if not isinstance(trig_channel, Channels.PhysicalMarkerChannel):
        raise ValueError("In set_master the trigger channel must be an instance of PhysicalMarkerChannel")

    cdb = Channels.ChannelDatabase[channelLib.channelDatabase.id] # Can't use external object
    st = Channels.LogicalMarkerChannel(label="slave_trig", channel_db=cdb)
    st.phys_chan = trig_channel
    st.pulse_params = {"length": pulse_length, "shape_fun": "constant"}
    transmitter.master = True
    transmitter.trigger_source = "Internal"

@db_session
def QubitFactory(label, **kwargs):
    ''' Return a saved qubit channel or create a new one. '''
    # TODO: this will just get the first entry in the whole damned DB!
    # thing = select(el for el in Channels.Qubit if el.label==label).first()
    thing = {c.label: c for c in channelLib.get_current_channels() if isinstance(c, Channels.Qubit)}[label]
    if thing:
        return thing
    else:
        return Channels.Qubit(label=label, **kwargs)

@db_session
def MeasFactory(label, **kwargs):
    ''' Return a saved measurement channel or create a new one. '''
    thing = {c.label: c for c in channelLib.get_current_channels() if isinstance(c, Channels.Measurement)}[label]
    if thing:
        return thing
    else:
        return Channels.Measurement(label=label, **kwargs)

@db_session
def MarkerFactory(label, **kwargs):
    ''' Return a saved Marker channel or create a new one. '''
    thing = {c.label: c for c in channelLib.get_current_channels()}[label]
    if thing:
        return thing
    else:
        return Channels.LogicalMarkerChannel(label=label, **kwargs)

@db_session
def EdgeFactory(source, target):
    if channelLib.connectivityG.has_edge(source, target):
        return channelLib.connectivityG[source][target]['channel']
    elif channelLib.connectivityG.has_edge(target, source):
        return channelLib.connectivityG[target][source]['channel']
    else:
        raise ValueError('Edge {0} not found in connectivity graph'.format((
            source, target)))
