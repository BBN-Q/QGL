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
from atom.api import Atom, Str, Int, Typed
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

from . import Channels
from . import PulseShapes
from . import config

channelLib = None

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

class MyEventHandler(FileSystemEventHandler):

    def __init__(self, file_paths, callback):
        super(MyEventHandler, self).__init__()
        self.file_paths = [os.path.normpath(fn) for fn in file_paths]
        self.callback = callback
        self.paused = True

        # The spotlight indexer in MacOSX retriggers events... maybe we should hash the files?
        self.grace_period = 3.0 if sys.platform == 'darwin' else 1.0 
        self.last_library_update = datetime.datetime.now()

    def on_modified(self, event):
        try:
            if any([os.path.samefile(event.src_path, fp) for fp in self.file_paths]):
                if not self.paused:
                    # Build in some sanity checking since we seem to get multiple
                    # events firing in a number of situations.
                    now = datetime.datetime.now()
                    
                    if (now-self.last_library_update).total_seconds() > (self.grace_period):
                        self.last_library_update = now
                        """
                        Hold off for half a second
                        If the event is from the file being opened to be written this gives
                        time for it to be written.
                        """
                        time.sleep(0.5)
                        self.callback()
        except FileNotFoundError:
            #Temporary settings files generated using yaml_dump get deleted
            #faster than the above code can catch it.
            pass

class LibraryFileWatcher(object):
    def __init__(self, main_path, callback):
        super(LibraryFileWatcher, self).__init__()
        
        self.main_path = os.path.normpath(main_path)
        self.callback = callback

        # Perform a preliminary loading to find all of the connected files...
        # TODO: modularity
        with open(os.path.abspath(self.main_path), 'r') as FID:
            loader = Loader(FID)
            try:
                tmpLib = loader.get_single_data()
                self.filenames = [os.path.normpath(lf) for lf in loader.filenames]
            finally:
                loader.dispose()

        self.eventHandler = MyEventHandler(self.filenames, self.callback)
        self.observer = Observer()
        self.observer.schedule(self.eventHandler, path=os.path.dirname(os.path.abspath(main_path)))

        self.observer.start()
        self.resume()

    def __del__(self):
        self.observer.stop()
        self.observer.join()

    def pause(self):
        self.eventHandler.paused = True

    def resume(self):
        self.eventHandler.paused = False

class ChannelLibrary(Atom):
    # channelDict = Dict(Str, Channel)
    channelDict = Typed(dict)
    connectivityG = Typed(nx.DiGraph)
    library_file = Str()
    fileWatcher = Typed(LibraryFileWatcher)
    version = Int(5)
    last_library_update = Str()

    specialParams = ['phys_chan', 'gate_chan', 'trig_chan', 'receiver_chan',
                     'source', 'target']

    def __init__(self, library_file=config.configFile, channelDict={}, **kwargs):
        """Create the channel library. We assume that the user wants the config file in the 
        usual locations specified in the config files."""
        if library_file:
            super(ChannelLibrary, self).__init__(channelDict=channelDict, library_file=library_file, **kwargs)
            self.connectivityG = nx.DiGraph()
            yaml_filenames = self.load_from_library()
            if self.library_file and yaml_filenames:
                self.fileWatcher = LibraryFileWatcher(self.library_file, self.update_from_file)
        else: # we want a blank library if library_file is none
            super(ChannelLibrary, self).__init__(channelDict={})
            self.connectivityG = nx.DiGraph()

        # Update the global reference
        global channelLib
        if channelLib:
            # Don't let the 
            channelLib.fileWatcher = None
        channelLib = self

        self.last_library_update = str(datetime.datetime.now())

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

    def load_from_library(self, return_only=False):
        """Loads the YAML library, creates the QGL objects, and returns a list of the visited filenames
        for the filewatcher."""
        if not self.library_file:
            return
        try:
            with open(self.library_file, 'r') as FID:
                loader = Loader(FID)
                try:
                    tmpLib = loader.get_single_data()
                    filenames = loader.filenames
                finally:
                    loader.dispose()

            # Check to see if we have the mandatory sections
            for section in ['instruments', 'qubits', 'filters']:
                if section not in tmpLib.keys():
                    raise ValueError("{} section not present in config file {}.".format(section, self.library_file))

            instr_dict   = tmpLib['instruments']
            qubit_dict   = tmpLib['qubits']
            filter_dict  = tmpLib['filters']
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
                    slave_chan = instr["slave_trig"] if "slave_trig" in instr.keys() else "slave"
                    master_awgs.append(name + "-" + slave_chan)
                # Eventually we should support multiple masters...
                # if "slave_trig" in instr.keys():
                #     params = {}
                #     params["label"]        = name + "_slave_trig"
                #     params["phys_chan"]    = name + "-" + instr["slave_trig"]
                #     params["pulse_params"] = {"length": 1e-7, "shape_fun": "constant"}
                #     params["__module__"]   = "QGL.Channels"
                #     params["__class__"]    = "LogicalMarkerChannel"
                #     print(params["label"], "***")
                #     channel_dict[params["label"]] = params

            # Establish the slave trigger, assuming for now that we have a single
            # APS master. This might change later.
            if len(master_awgs) > 1:
                raise ValueError("More than one AWG is marked as master.")
            elif len(master_awgs) == 1:
                params = {}
                params["label"]       = "slave_trig"
                params["phys_chan"]    = master_awgs[0]
                if params["phys_chan"] in marker_lens.keys():
                    length = marker_lens[params["phys_chan"]]
                else:
                    length = 1e-7
                params["pulse_params"] = {"length": length, "shape_fun": "constant"}
                params["__module__"]  = "QGL.Channels"
                params["__class__"]   = "LogicalMarkerChannel"
                channel_dict[params["label"]] = params

            for name, filt in filter_dict.items():
                if "StreamSelector" in filt["type"]:
                    params = {k: v for k,v in filt.items() if k in Channels.ReceiverChannel.__atom_members__.keys()}
                    params["label"]      = "RecvChan-" + name # instr_dict[filt["instrument"]]["name"] + "-" + name
                    params["channel"]    = str(params["channel"]) # Convert to a string
                    params["instrument"] = filt["source"]
                    params["__module__"] = "QGL.Channels"
                    params["__class__"]  = "ReceiverChannel"
                    if "source" not in filt.keys():
                        raise ValueError("No instrument (source) specified for Stream Selector")
                    if filt["source"] not in instr_dict.keys() and filt["source"] not in channel_dict.keys():
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
                params["trig_chan"]     = "digTrig-" + qubit["measure"]["trigger"]
                params["phys_chan"]     = meas_instr + "-" + meas_chan
                params["meas_type"]     = "autodyne"
                params["receiver_chan"] = "RecvChan-" + qubit["measure"]["receiver"]
                params["__module__"]   = "QGL.Channels"
                params["__class__"]    = "Measurement"
                channel_dict[params["label"]] = params
                if 'generator' in qubit["measure"].keys():
                    channel_dict[params["phys_chan"]]["generator"] = qubit["measure"]["generator"]

                # Create the receiver channels
                if "receiver" in qubit["measure"].keys():
                    if len(qubit["measure"]["receiver"].split()) != 1:
                        print("Receiver specification for {} ({}) must have a stream selector".format(name, qubit["measure"]["receiver"]))
                        raise ValueError("Receiver specification for {} ({}) must have a stream selector".format(name, qubit["measure"]["receiver"]))
                    phys_instr, phys_marker = qubit["measure"]["trigger"].split()
                    params = {}
                    params["label"]        = "digTrig-" + qubit["measure"]["trigger"]
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

                # Create the measurement gate chan:
                if "gate" in qubit["measure"].keys():
                    phys_instr, phys_marker = qubit["measure"]["gate"].split()
                    params = {}
                    params["label"]      = "M-{}-gate".format(name)
                    params["phys_chan"]   = phys_instr + "-" + phys_marker
                    params["__module__"] = "QGL.Channels"
                    params["__class__"]  = "LogicalMarkerChannel"
                    channel_dict[params["label"]] = params
                    channel_dict["M-{}".format(name)]["gate_chan"] = params["label"]

                # Create the control gate chan:
                if "gate" in qubit["control"].keys():
                    phys_instr, phys_marker = qubit["control"]["gate"].split()
                    params = {}
                    params["label"]      = "{}-gate".format(name)
                    params["phys_chan"]   = phys_instr + "-" + phys_marker
                    params["__module__"] = "QGL.Channels"
                    params["__class__"]  = "LogicalMarkerChannel"
                    channel_dict[params["label"]] = params
                    channel_dict[name]["gate_chan"] = params["label"]


            for trig_name, trigger in trigger_dict.items():
                phys_instr, phys_marker = trigger.split()
                params = {}
                params["label"]      = trig_name
                params["phys_chan"]   = phys_instr + "-" + phys_marker
                if params["phys_chan"] in marker_lens.keys():
                    length = marker_lens[params["phys_chan"]]
                else:
                    length = 1e-7
                params["__module__"] = "QGL.Channels"
                params["__class__"]  = "LogicalMarkerChannel"
                channel_dict[params["label"]] = params

            for name, edge in edge_dict.items():
                # Create the Edges
                if len(edge["AWG"].split()) != 2:
                    print("Control AWG specification for {} ({}) must have a device, channel".format(name, edge["AWG"]))
                    raise ValueError("Control AWG specification for {} ({}) must have a device, channel".format(name, edge["AWG"]))
                ctrl_instr, ctrl_chan = edge["AWG"].split()
                params = {k: v for k,v in edge.items() if k in Channels.Edge.__atom_members__.keys()}
                params["label"]      = name
                params["phys_chan"]   = ctrl_instr + "-" + ctrl_chan
                params["__module__"] = "QGL.Channels"
                params["__class__"]  = "Edge"
                channel_dict[params["label"]] = params
                if 'generator' in edge.keys():
                    channel_dict[params["phys_chan"]]["generator"] = edge["generator"]

                # Create the edge gate chan:
                if "gate" in edge.keys():
                    phys_instr, phys_marker = edge["gate"].split()
                    params = {}
                    params["label"]      = "{}-gate".format(name)
                    params["phys_chan"]   = phys_instr + "-" + phys_marker
                    params["__module__"] = "QGL.Channels"
                    params["__class__"]  = "LogicalMarkerChannel"
                    channel_dict[params["label"]] = params
                    channel_dict[name]["gate_chan"] = params["label"]

            if return_only:
                return channel_dict
            else:
                channel_dict = {k: self.instantiate(v) for k,v in channel_dict.items()}
                # connect objects labeled by strings
                for chan in channel_dict.values():
                    for param in self.specialParams:
                        if hasattr(chan, param) and getattr(chan, param) is not None:
                            chan_to_find = channel_dict.get(getattr(chan, param), None)
                            if not chan_to_find:
                                print("Couldn't find {} of {} in the channel_dict!".format(param, chan))
                            setattr(chan, param, chan_to_find)

                self.channelDict.update(channel_dict)
                self.build_connectivity_graph()
                return filenames

        except IOError:
            print('No channel library found.')
        except Exception as e:
            print('Failed to load channel library: received exception', e)
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_tb(exc_traceback, limit=4, file=sys.stdout)

    def instantiate(self, paramDict):
        if 'pulse_params' in paramDict:
            if 'shape_fun' in paramDict['pulse_params']:
                shape_fun = paramDict['pulse_params']['shape_fun']
                paramDict['pulse_params']['shape_fun'] = getattr(PulseShapes, shape_fun)
        if '__class__' in paramDict:
            className  = paramDict.pop('__class__')
            moduleName = paramDict.pop('__module__')
            __import__(moduleName)
            return getattr(sys.modules[moduleName], className)(**paramDict)


    def update_from_file(self):
        """
        Only update relevant parameters
        Helps avoid both stale references from replacing whole channel objects (as in load_from_library)
        and the overhead of recreating everything.
        """

        if not self.library_file:
            return
        try:
            all_params = self.load_from_library(return_only=True)

            # update & insert
            for chName, chParams in all_params.items():
                if chName in self.channelDict:
                    self.update_from_json(chName, chParams)
                else:
                    self.channelDict[chName] = self.instantiate(chParams)
                    self.update_from_json(chName, chParams)

            # remove
            for chName in list(self.channelDict.keys()):
                if chName not in all_params:
                    del self.channelDict[chName]

            self.build_connectivity_graph()

            print("Updated library")
            self.last_library_update = str(datetime.datetime.now())
        except:
            print('Failed to update channel library from file. Is there a typo?.')
            return

        # reset pulse cache
        from . import PulsePrimitives
        PulsePrimitives._memoize.cache.clear()


    def update_from_json(self, chName, chParams):
        # connect objects labeled by strings
        if 'pulse_params' in chParams.keys():
            paramDict = {str(k): v for k, v in chParams['pulse_params'].items()}
            shapeFunName = paramDict.pop('shape_fun', None)
            if shapeFunName:
                paramDict['shape_fun'] = getattr(PulseShapes, shapeFunName)
            self.channelDict[chName].pulse_params = paramDict

        for param in self.specialParams:
            if param in chParams.keys():
                setattr(self.channelDict[chName],
                        param,
                        self.channelDict.get(chParams[param], None)
                        )
        # TODO: how do we follow changes to selected AWG or generator?

        # ignored or specially handled parameters
        ignoreList = self.specialParams + ['pulse_params', 'AWG', 'generator', '__class__', '__module__']
        for paramName in chParams:
            if paramName not in ignoreList:
                setattr(self.channelDict[chName], paramName, chParams[paramName])

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

def MarkerFactory(label, **kwargs):
    '''Return a marker channel by name. Must be defined under top-level `markers`
    keyword in measurement configuration YAML.
    '''
    if not channelLib:
        raise ValueError('ChannelLibrary not found, has an instance of ChannelLibrary been created?')
    if label in channelLib and isinstance(channelLib[label], Channels.LogicalMarkerChannel):
        return channelLib[label]
    else:
        raise ValueError("Marker channel {} not found in channel library.".format(label))

def QubitFactory(label, **kwargs):
    ''' Return a saved qubit channel or create a new one. '''
    if not channelLib:
        raise ValueError('ChannelLibrary not found, has an instance of ChannelLibrary been created?')
    if label in channelLib and isinstance(channelLib[label], Channels.Qubit):
        return channelLib[label]
    else:
        return Channels.Qubit(label=label, **kwargs)

def MeasFactory(label, meas_type='autodyne', **kwargs):
    ''' Return a saved measurement channel or create a new one. '''
    if not channelLib:
        raise ValueError('ChannelLibrary not found, has an instance of ChannelLibrary been created?')
    if label in channelLib and isinstance(channelLib[label], Channels.Measurement):
        return channelLib[label]
    else:
        return Channels.Measurement(label=label, meas_type=meas_type, **kwargs)

def EdgeFactory(source, target):
    if not channelLib:
        raise ValueError('Connectivity graph not found. Has a ChannelLibrary has been created?')
    if channelLib.connectivityG.has_edge(source, target):
        return channelLib.connectivityG[source][target]['channel']
    elif channelLib.connectivityG.has_edge(target, source):
        return channelLib.connectivityG[target][source]['channel']
    else:
        raise ValueError('Edge {0} not found in connectivity graph'.format((
            source, target)))
