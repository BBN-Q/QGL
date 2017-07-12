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
import importlib
from atom.api import Atom, Str, Int, Typed
import networkx as nx
import yaml

from . import Channels
from . import PulseShapes
from JSONLibraryUtils import FileWatcher
from . import config

class LoaderMeta(type):
    def __new__(metacls, __name__, __bases__, __dict__):
        """Add include constructer to class."""
        # register the include constructor on the class
        cls = super().__new__(metacls, __name__, __bases__, __dict__)
        cls.add_constructor('!include', cls.construct_include)
        return cls
class Loader(yaml.Loader, metaclass=LoaderMeta):
    """YAML Loader with `!include` constructor."""
    def __init__(self, stream):
        """Initialise Loader."""
        try:
            self._root = os.path.split(stream.name)[0]
        except AttributeError:
            self._root = os.path.curdir
        super().__init__(stream)
        self.add_implicit_resolver(
            u'tag:yaml.org,2002:float',
            re.compile(u'''^(?:
             [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
            |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
            |\\.[0-9_]+(?:[eE][-+][0-9]+)?
            |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
            |[-+]?\\.(?:inf|Inf|INF)
            |\\.(?:nan|NaN|NAN))$''', re.X),
            list(u'-+0123456789.'))
        self.filenames = [os.path.abspath(stream.name)]
    def construct_include(self, node):
        """Include file referenced at node."""
        filename = os.path.abspath(os.path.join(
            self._root, self.construct_scalar(node)
        ))
        extension = os.path.splitext(filename)[1].lstrip('.')
        self.filenames.append(filename)
        with open(filename, 'r') as f:
            if extension in ('yaml', 'yml'):
                return yaml.load(f, Loader)
            else:
                return ''.join(f.readlines())

class ChannelLibrary(Atom):
    # channelDict = Dict(Str, Channel)
    channelDict = Typed(dict)
    connectivityG = Typed(nx.DiGraph)
    libFile = Str()
    fileWatcher = Typed(FileWatcher.LibraryFileWatcher)
    version = Int(5)

    specialParams = ['physChan', 'gateChan', 'trigChan', 'receiverChan',
                     'source', 'target']

    def __init__(self, channelDict={}, **kwargs):
        super(ChannelLibrary, self).__init__(channelDict=channelDict, **kwargs)
        self.connectivityG = nx.DiGraph()
        yaml_filenames = self.load_from_library()
        if self.libFile and yaml_filenames:
            self.fileWatcher = FileWatcher.LibraryFileWatcher(
                self.libFile, self.update_from_file)

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
        self.connectivityG.clear()
        for chan in self.channelDict.values():
            if isinstance(chan,
                          Channels.Qubit) and chan not in self.connectivityG:
                self.connectivityG.add_node(chan)
        for chan in self.channelDict.values():
            if isinstance(chan, Channels.Edge):
                self.connectivityG.add_edge(chan.source, chan.target)
                self.connectivityG[chan.source][chan.target]['channel'] = chan

    def load_from_library(self):
        """Loads the YAML library, creates the QGL objects, and returns a list of the visited filenames 
        for the filewatcher."""
        if not self.libFile:
            return
        try:
            with open(self.libFile, 'r') as FID:
                loader = Loader(FID)
                try:
                    tmpLib = loader.get_single_data()
                    filenames = loader.filenames
                finally:
                    loader.dispose()

            instr_dict  = tmpLib['instruments']
            qubit_dict  = tmpLib['qubits']
            filter_dict = tmpLib['filters']
            master_awg  = []

            # Construct the channel library
            channel_dict = {}

            for name, instr in instr_dict.items():
                if "channels" in instr.keys():
                    for channel in instr["channels"]:
                        params = {k: v for k,v in channel.items() if k in Channels.PhysicalQuadratureChannel.__atom_members__.keys()}
                        params["label"] = name + "-" + channel["name"]
                        params["instrument"] = name
                        params["translator"] = instr["type"] + "Pattern"
                        params["__module__"] = "QGL.Channels"
                        params["__class__"]  = "PhysicalQuadratureChannel"
                        channel_dict[params["label"]] = params
                if "markers" in instr.keys():
                    for marker in instr["markers"]:
                        params = {k: v for k,v in marker.items() if k in Channels.PhysicalMarkerChannel.__atom_members__.keys()}
                        params["label"] = name + "-" + marker["name"]
                        params["instrument"] = name
                        params["translator"] = instr["type"] + "Pattern"
                        params["__module__"] = "QGL.Channels"
                        params["__class__"]  = "PhysicalMarkerChannel"
                        channel_dict[params["label"]] = params
                if "master" in instr.keys():
                    if instr["master"]:
                        master_awg.append(name + "-" + instr["slaveTrig"])

            # Establish the slave trigger
            if len(master_awg) != 1:
                raise ValueError("More (or less) than one AWG is marked as master.")
            else:
                params = {}
                params["label"]       = "slaveTrig"
                params["physChan"]    = master_awg[0]
                params["pulseParams"] = {"length": 1e-7, "shapeFun": "constant"}
                params["__module__"]  = "QGL.Channels"
                params["__class__"]   = "LogicalMarkerChannel"
                channel_dict[params["label"]] = params

            for name, filt in filter_dict.items():
                if "StreamSelector" in filt["type"]:
                    params = {k: v for k,v in filt.items() if k in Channels.ReceiverChannel.__atom_members__.keys()}
                    params["label"]      = "RecvChan-" + name # instr_dict[filt["instrument"]]["name"] + "-" + name
                    params["__module__"] = "QGL.Channels"
                    params["__class__"]  = "ReceiverChannel"

                    if "source" not in filt.keys():
                        raise ValueError("No instrument (source) specified for Stream Selector")
                    if filt["source"] not in instr_dict.keys():
                        raise ValueError("Stream Selector source {} not found among list of instruments.".format(filt["source"]))
                    params["instrument"] = filt["source"]

                    channel_dict[params["label"]] = params

            for name, qubit in qubit_dict.items():
                # Create the Qubits
                if len(qubit["control"]["AWG"].split()) != 2:
                    print("Control AWG specification for {} ({}) must have a device, channel".format(name, qubit["control"]["AWG"]))
                    raise ValueError("Control AWG specification for {} ({}) must have a device, channel".format(name, qubit["control"]["AWG"]))
                ctrl_instr, ctrl_chan = qubit["control"]["AWG"].split()
                params = {k: v for k,v in qubit["control"].items() if k in Channels.Qubit.__atom_members__.keys()}
                params["label"]      = name
                params["physChan"]   = ctrl_instr + "-" + ctrl_chan
                params["__module__"] = "QGL.Channels"
                params["__class__"]  = "Qubit"
                channel_dict[params["label"]] = params
                if 'generator' in qubit["control"].keys():
                    channel_dict[params["physChan"]]["generator"] = qubit["control"]["generator"]

                # Create the measurements
                if len(qubit["measure"]["AWG"].split()) != 2:
                    print("Measurement AWG specification for {} ({}) must have a device, channel".format(name, qubit["measure"]["AWG"]))
                    raise ValueError("Measurement AWG specification for {} ({}) must have a device, channel".format(name, qubit["measure"]["AWG"]))
                meas_instr, meas_chan = qubit["measure"]["AWG"].split()
                params = {k: v for k,v in qubit["measure"].items() if k in Channels.Measurement.__atom_members__.keys()}
                params["label"]        = "M-{}".format(name)
                params["trigChan"]     = "digTrig-" + meas_instr
                params["physChan"]     = meas_instr + "-" + meas_chan
                params["measType"]     = "autodyne"
                params["receiverChan"] = "RecvChan-" + qubit["measure"]["receiver"]
                params["__module__"]   = "QGL.Channels"
                params["__class__"]    = "Measurement"
                channel_dict[params["label"]] = params
                if 'generator' in qubit["measure"].keys():
                    channel_dict[params["physChan"]]["generator"] = qubit["measure"]["generator"]

                # Create the receiver channels
                if "receiver" in qubit["measure"].keys():
                    if len(qubit["measure"]["receiver"].split()) != 1:
                        print("Receiver specification for {} ({}) must have a stream selector".format(name, qubit["measure"]["receiver"]))
                        raise ValueError("Receiver specification for {} ({}) must have a stream selector".format(name, qubit["measure"]["receiver"]))
                    phys_instr, phys_marker = qubit["measure"]["trigger"].split()
                    params = {k: v for k,v in marker.items() if k in Channels.LogicalMarkerChannel.__atom_members__.keys()}
                    params["label"]        = "digTrig-" + phys_instr
                    params["physChan"]     = phys_instr + "-" + phys_marker
                    params["pulseParams"]  = {"length": 3e-7, "shapeFun": "constant"}
                    params["__module__"]   = "QGL.Channels"
                    params["__class__"]    = "LogicalMarkerChannel"
                    # Don't duplicate triggers to the same digitizer
                    if params["label"] not in channel_dict.keys():
                        channel_dict[params["label"]] = params

                # Create the measurement gate chan:
                if "gate" in qubit["measure"].keys():
                    phys_instr, phys_marker = qubit["measure"]["gate"].split()
                    params = {}
                    params["label"]      = "M-{}-gate".format(name)
                    params["physChan"]   = phys_instr + "-" + phys_marker
                    params["__module__"] = "QGL.Channels"
                    params["__class__"]  = "LogicalMarkerChannel"
                    channel_dict[params["label"]] = params

                # Create the control gate chan:
                if "gate" in qubit["control"].keys():
                    phys_instr, phys_marker = qubit["control"]["gate"].split()
                    params = {}
                    params["label"]      = "{}-gate".format(name)
                    params["physChan"]   = phys_instr + "-" + phys_marker
                    params["__module__"] = "QGL.Channels"
                    params["__class__"]  = "LogicalMarkerChannel"
                    channel_dict[params["label"]] = params

            def instantiate(paramDict):
                if 'pulseParams' in paramDict:
                    if 'shapeFun' in paramDict['pulseParams']:
                        shapeFun = paramDict['pulseParams']['shapeFun']
                        paramDict['pulseParams']['shapeFun'] = getattr(PulseShapes, shapeFun)
                if '__class__' in paramDict:
                    className  = paramDict.pop('__class__')
                    moduleName = paramDict.pop('__module__')
                    __import__(moduleName)
                    return getattr(sys.modules[moduleName], className)(**paramDict)

            channel_dict = {k: instantiate(v) for k,v in channel_dict.items()}
            # connect objects labeled by strings
            for chan in channel_dict.values():
                for param in self.specialParams:
                    if hasattr(chan, param):
                        setattr(chan, param,
                                channel_dict.get(getattr(chan, param), None)
                                )
            self.channelDict.update(channel_dict)
            self.build_connectivity_graph()
            return filenames

        except IOError:
            print('No channel library found.')
        except ValueError:
            print('Failed to load channel library.')

    def update_from_file(self):
        if not self.libFile:
            return
        try:
            self.load_from_library()
        except:
            print('Failed to update channel library from file. Probably is just half-written.')
            return

        # reset pulse cache
        from . import PulsePrimitives
        PulsePrimitives._memoize.cache.clear()

    def on_awg_change(self, oldName, newName):
        print("Change AWG", oldName, newName)
        for chName in self.channelDict:
            if isinstance(self.channelDict[chName],
                          (Channels.PhysicalMarkerChannel,
                           Channels.PhysicalQuadratureChannel)):
                awgName, awgChannel = chName.rsplit('-', 1)
                if awgName == oldName:
                    newLabel = "{0}-{1}".format(newName, awgChannel)
                    print("Changing {0} to {1}".format(chName, newLabel))
                    self.physicalChannelManager.name_changed(chName, newLabel)

def QubitFactory(label, **kwargs):
    ''' Return a saved qubit channel or create a new one. '''
    if channelLib and label in channelLib and isinstance(channelLib[label],
                                                         Channels.Qubit):
        return channelLib[label]
    else:
        return Channels.Qubit(label=label, **kwargs)


def MeasFactory(label, measType='autodyne', **kwargs):
    ''' Return a saved measurement channel or create a new one. '''
    if channelLib and label in channelLib and isinstance(channelLib[label],
                                                         Channels.Measurement):
        return channelLib[label]
    else:
        return Channels.Measurement(label=label, measType=measType, **kwargs)


def EdgeFactory(source, target):
    if not channelLib:
        raise ValueError('Connectivity graph not found')
    if channelLib.connectivityG.has_edge(source, target):
        return channelLib.connectivityG[source][target]['channel']
    elif channelLib.connectivityG.has_edge(target, source):
        return channelLib.connectivityG[target][source]['channel']
    else:
        raise ValueError('Edge {0} not found in connectivity graph'.format((
            source, target)))

channelLib = ChannelLibrary(libFile=config.configFile)
