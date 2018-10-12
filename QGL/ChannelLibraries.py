'''
Channels is where we store information for mapping virtual (qubit) channel to
real channels.

Split from Channels.py on Jan 14, 2016.
Moved to SQLAlchemy ORM from atom 2018

Original Author: Colm Ryan
Modified By: Graham Rowlands

Copyright 2016-2018 Raytheon BBN Technologies

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
import numpy as np
import networkx as nx

import bbndb

from . import config
from . import Channels
from . import PulseShapes
from .PulsePrimitives import clear_pulse_cache

from sqlalchemy.orm.session import make_transient
from sqlalchemy.pool import StaticPool

channelLib = None

def check_session_dirty(f):
    """Since we can't mix db objects from separate sessions, re-fetch entities by their unique IDs"""
    @wraps(f)
    def wrapper(cls, *args, **kwargs):
        if 'force' in kwargs and kwargs['force']:
            return f(cls, *args, **kwargs)
        elif len(cls.session.dirty | cls.session.new) != 0:
            raise Exception("Uncommitted transactions for working database. Either use force=True or commit/revert your changes.")
    return wrapper

class ChannelLibrary(object):

    def __init__(self, db_resource_name=None):
        """Create the channel library."""

        global channelLib

        self.db_provider = "sqlite"
        self.db_resource_name = ":memory:"

        if bbndb.engine:
            # Use current db
            self.db = bbndb.engine
        else:

            if db_resource_name:
                self.db_resource_name = db_resource_name
            elif config.load_db():
                self.db_resource_name = config.load_db()

            self.db = bbndb.engine = bbndb.create_engine(f'{self.db_provider}:///{self.db_resource_name}',
                                                            connect_args={'check_same_thread':False},
                                                            poolclass=StaticPool,
                                                            echo=False)

        bbndb.Base.metadata.create_all(bbndb.engine)
        bbndb.Session.configure(bind=bbndb.engine)
        self.Session = bbndb.Session
        self.session = self.Session()

        self.connectivityG = nx.DiGraph()

        # Check to see whether there is already a temp database
        working_dbs = self.query(Channels.ChannelDatabase, label="working").all()
        if len(working_dbs) > 1:
            raise Exception("More than one working database exists!")
        elif len(working_dbs) == 1:
            self.channelDatabase = working_dbs[0]
        elif len(working_dbs) == 0:
            self.channelDatabase = Channels.ChannelDatabase(label="working", time=datetime.datetime.now())
            self.session.add(self.channelDatabase)
            self.session.commit()

        self.update_channelDict()

        # Update the global reference
        channelLib = self

    def query(self, obj_type, **kwargs):
        return self.session.query(obj_type).filter_by(**kwargs)

    def get_current_channels(self):
        return self.channelDatabase.channels + self.channelDatabase.generators

    def update_channelDict(self):
        self.channelDict = {c.label: c for c in self.get_current_channels()}

    def ls(self):
        cdb = Channels.ChannelDatabase
        q = self.session.query(cdb.label, cdb.time, cdb.id).\
            order_by(Channels.ChannelDatabase.id, Channels.ChannelDatabase.label).all()
        for i, (label, time, id) in enumerate(q):
                t = time.strftime("(%Y) %b. %d @ %I:%M:%S %p")
                print(f"[{id}] {t} -> {label}")

    def ent_by_type(self, obj_type, show=False):
        q = self.session.query(obj_type).filter(obj_type.channel_db.has(label="working")).order_by(obj_type.label).all()
        if show:
            for i, el in enumerate(q):
                print(f"[{i}] -> {el.label}")
        else:
            return q

    def receivers(self):
        return self.ent_by_type(Channels.Receiver)

    def transmitters(self):
        return self.ent_by_type(Channels.Transmitter)

    def transceivers(self):
        return self.ent_by_type(Channels.Transceiver)

    def qubits(self):
        return self.ent_by_type(Channels.Qubit)

    def meas(self):
        return self.ent_by_type(Channels.Measurement)

    def ls_receivers(self):
        return self.ent_by_type(Channels.Receiver, show=True)

    def ls_transmitters(self):
        return self.ent_by_type(Channels.Transmitter, show=True)

    def ls_qubits(self):
        return self.ent_by_type(Channels.Qubit, show=True)

    def ls_measurements(self):
        return self.ent_by_type(Channels.Measurement, show=True)

    @check_session_dirty
    def load(self, name, index=1):
        """Load the latest instance for a particular name. Specifying index = 2 will select the second most recent instance """
        cdb = Channels.ChannelDatabase
        items = self.session.query(cdb).filter(cdb.label==name).order_by(cdb.time.desc()).all()
        self.load_obj(items[-index])

    @check_session_dirty
    def load_by_id(self, id_num):
        cdb = Channels.ChannelDatabase
        item = self.session.query(cdb).filter(cdb.id==id_num).first()
        self.load_obj(item)

    def clear(self, channel_db=None, create_new=True):
        # If no database is specified, clear self.database
        channel_db = channel_db if channel_db else self.channelDatabase

        self.session.delete(channel_db)
        self.session.commit()

        if create_new:
            self.channelDatabase = Channels.ChannelDatabase(label="working", time=datetime.datetime.now())
            self.session.add(self.channelDatabase)
            self.session.commit()
        channelLib = self

    def load_obj(self, obj):
        self.clear(create_new=False)
        self.channelDatabase = bbndb.deepcopy_sqla_object(obj, self.session)
        self.channelDatabase.label = "working"
        self.update_channelDict()

    def commit(self):
        self.session.commit()

    def revert(self):
        self.session.rollback()

    def save_as(self, name):
        if name == "working":
            raise ValueError("Cannot save as `working` since that is the default working environment name...")
        self.commit()
        new_channelDatabase = bbndb.deepcopy_sqla_object(self.channelDatabase, self.session)
        new_channelDatabase.label = name
        new_channelDatabase.time = datetime.datetime.now()
        self.commit()

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
        for chan in select(q for q in Channels.Qubit if q not in self.connectivityG):
=======
        for chan in self.session.query(Channels.Qubit).filter(Channels.Qubit not in self.connectivityG).all():
            self.connectivityG.add_node(chan)
        for chan in self.session.query(Channels.Edge): #select(e for e in Channels.Edge):
            self.connectivityG.add_edge(chan.source, chan.target)
            self.connectivityG[chan.source][chan.target]['channel'] = chan
>>>>>>> Ditch atom, move to Pony.orm for all channel library objects.

    def new_APS2(self, label, address):
        chan12 = Channels.PhysicalQuadratureChannel(label=f"{label}-12", instrument=label, translator="APS2Pattern", channel_db=self.channelDatabase)
        m1     = Channels.PhysicalMarkerChannel(label=f"{label}-12m1", instrument=label, translator="APS2Pattern", channel_db=self.channelDatabase)
        m2     = Channels.PhysicalMarkerChannel(label=f"{label}-12m2", instrument=label, translator="APS2Pattern", channel_db=self.channelDatabase)
        m3     = Channels.PhysicalMarkerChannel(label=f"{label}-12m3", instrument=label, translator="APS2Pattern", channel_db=self.channelDatabase)
        m4     = Channels.PhysicalMarkerChannel(label=f"{label}-12m4", instrument=label, translator="APS2Pattern", channel_db=self.channelDatabase)

        this_transmitter = Channels.Transmitter(label=label, model="APS2", address=address, channels=[chan12, m1, m2, m3, m4], channel_db=self.channelDatabase)
        this_transmitter.trigger_source = "external"
        this_transmitter.address        = address

        self.session.add(this_transmitter)
        return this_transmitter

    def new_APS2_rack(self, label, num, start_address):
        address_start    = ".".join(start_address.split(".")[:3])
        address_end      = int(start_address.split(".")[-1])
        transmitters     = [new_APS2(f"{label}_U{i}", f"{address_start}.{address_end+i}") for i in range(1,num+1)]
        this_transceiver = Channels.Transceiver(label=label, model="APS2Rack", transmitters=transmitters, channel_db=self.channelDatabase)

        self.session.add(this_transceiver)
        return this_transceiver

    def new_X6(self, label, address, dsp_channel=0, record_length=1024):
        chan1 = Channels.ReceiverChannel(label=f"RecvChan-{label}-1", channel=1, dsp_channel=dsp_channel, channel_db=self.channelDatabase)
        chan2 = Channels.ReceiverChannel(label=f"RecvChan-{label}-2", channel=2, dsp_channel=dsp_channel, channel_db=self.channelDatabase)

        this_receiver = Channels.Receiver(label=label, model="X6-1000M", address=address, channels=[chan1, chan2],
                                      record_length=record_length, channel_db=self.channelDatabase)
        this_receiver.trigger_source = "external"
        this_receiver.stream_types   = "raw, demodulated, integrated"
        this_receiver.address        = address

        # Add a default kernel
        chan1.kernel = np.ones(record_length, dtype=np.complex).tobytes()
        chan2.kernel = np.ones(record_length, dtype=np.complex).tobytes()

        self.session.add(this_receiver)
        return this_receiver

    def new_Alazar(self, label, address, record_length=1024):
        chan1 = Channels.ReceiverChannel(label=f"RecvChan-{label}-1", channel=1, channel_db=self.channelDatabase)
        chan2 = Channels.ReceiverChannel(label=f"RecvChan-{label}-2", channel=2, channel_db=self.channelDatabase)

        this_receiver = Channels.Receiver(label=label, model="AlazarATS9870", address=address, channels=[chan1, chan2],
                                      record_length=record_length, channel_db=self.channelDatabase)
        this_receiver.trigger_source = "external"
        this_receiver.stream_types   = "raw"
        this_receiver.address        = address

        self.session.add(this_receiver)
        return this_receiver

    def new_qubit(self, label):
        thing = Channels.Qubit(label=label, channel_db=self.channelDatabase)
        self.session.add(thing)
        return thing

    def new_source(self, label, model, address, power=-30.0, frequency=5.0e9):
        thing = Channels.Generator(label=label, model=model, address=address, power=power, frequency=frequency, channel_db=self.channelDatabase)
        self.session.add(thing)
        return thing

    def set_control(self, qubit, transmitter, generator=None):
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

    def set_measure(self, qubit, transmitter, receivers, generator=None, trig_channel=None, gate=False, gate_channel=None, trigger_length=1e-7):
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

        meas = Channels.Measurement(label=f"M-{qubit.label}", channel_db=self.channelDatabase)
        meas.phys_chan = phys_chan
        if generator:
            meas.phys_chan.generator = generator

        phys_trig_channel = trig_channel if trig_channel else transmitter.get_chan("12m1")

        trig_chan              = Channels.LogicalMarkerChannel(label=f"receiversTrig-{qubit.label}", channel_db=self.channelDatabase)
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
            gate_chan           = Channels.LogicalMarkerChannel(label=f"M-{qubit.label}-gate", channel_db=self.channelDatabase)
            gate_chan.phys_chan = phys_gate_channel
            meas.gate_chan      = gate_chan

    def set_master(self, transmitter, trig_channel, pulse_length=1e-7):
        if not isinstance(trig_channel, Channels.PhysicalMarkerChannel):
            raise ValueError("In set_master the trigger channel must be an instance of PhysicalMarkerChannel")

        st = Channels.LogicalMarkerChannel(label="slave_trig", channel_db=self.channelDatabase)
        st.phys_chan = trig_channel
        st.pulse_params = {"length": pulse_length, "shape_fun": "constant"}
        transmitter.master = True
        transmitter.trigger_source = "internal"

def QubitFactory(label, **kwargs):
    ''' Return a saved qubit channel or create a new one. '''
    q = channelLib.session.query(Channels.Qubit).filter(Channels.Qubit.label==label).all()
    if len(q) == 1:
        return q[0]
    else:
        c = Channels.Qubit(label=label, channel_db=channelLib.channelDatabase, **kwargs)
        channelLib.session.add(c)
        return c

def MeasFactory(label, **kwargs):
    ''' Return a saved measurement channel or create a new one. '''
    q = channelLib.session.query(Channels.Measurement).filter(Channels.Measurement.label==label).all()
    if len(q) == 1:
        return q[0]
    else:
        c = Channels.Measurement(label=label, channel_db=channelLib.channelDatabase, **kwargs)
        channelLib.session.add(c)
        return c

def MarkerFactory(label, **kwargs):
    ''' Return a saved Marker channel or create a new one. '''
    q = channelLib.session.query(Channels.LogicalMarkerChannel).filter(Channels.LogicalMarkerChannel.label==label).all()
    if len(q) == 1:
        return q[0]
    else:
        c = Channels.LogicalMarkerChannel(label=label, channel_db=channelLib.channelDatabase, **kwargs)
        channelLib.session.add(c)
        return c

def EdgeFactory(source, target):
    if channelLib.connectivityG.has_edge(source, target):
        return channelLib.connectivityG[source][target]['channel']
    elif channelLib.connectivityG.has_edge(target, source):
        return channelLib.connectivityG[target][source]['channel']
    else:
        raise ValueError('Edge {0} not found in connectivity graph'.format((
            source, target)))
