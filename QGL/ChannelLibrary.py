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
'''

import sys
import json
import importlib
from atom.api import Atom, Str, Int, Typed
import networkx as nx

from . import Channels
from . import PulseShapes
from JSONLibraryUtils import LibraryCoders, FileWatcher, JSONMigrators
from . import config


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
        self.load_from_library()
        if self.libFile:
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

    def write_to_file(self, fileName=None):
        libFileName = fileName if fileName != None else self.libFile
        if libFileName:
            #Pause the file watcher to stop cicular updating insanity
            if self.fileWatcher:
                self.fileWatcher.pause()
            with open(libFileName, 'w') as FID:
                json.dump(self,
                          FID,
                          cls=LibraryCoders.LibraryEncoder,
                          indent=2,
                          sort_keys=True)
            if self.fileWatcher:
                self.fileWatcher.resume()

    def json_encode(self, matlabCompatible=False):
        return {"channelDict": self.channelDict, "version": self.version}

    def load_from_library(self):
        if not self.libFile:
            return
        try:
            with open(self.libFile, 'r') as FID:
                tmpLib = json.load(FID, cls=ChannelDecoder)
            if not isinstance(tmpLib, ChannelLibrary):
                raise ValueError('Failed to load channel library')

            # connect objects labeled by strings
            for chan in tmpLib.channelDict.values():
                for param in self.specialParams:
                    if hasattr(chan, param):
                        setattr(chan, param,
                                tmpLib.channelDict.get(getattr(chan, param), None)
                                )
            self.channelDict.update(tmpLib.channelDict)
            # grab library version
            self.version = tmpLib.version
            self.build_connectivity_graph()
        except IOError:
            print('No channel library found.')
        except ValueError:
            print('Failed to load channel library.')

    def update_from_file(self):
        """
        Only update relevant parameters
        Helps avoid both stale references from replacing whole channel objects (as in load_from_library)
        and the overhead of recreating everything.
        """

        if not self.libFile:
            return
        try:
            with open(self.libFile, 'r') as FID:
                allParams = json.load(FID)['channelDict']

            # update & insert
            for chName, chParams in allParams.items():
                if chName in self.channelDict:
                    self.update_from_json(chName, chParams)
                else:
                    # load class from name and update from json
                    className = chParams['x__class__']
                    moduleName = chParams['x__module__']

                    mod = importlib.import_module(moduleName)
                    cls = getattr(mod, className)
                    self.channelDict[chName] = cls()
                    self.update_from_json(chName, chParams)

            # remove
            for chName in list(self.channelDict.keys()):
                if chName not in allParams:
                    del self.channelDict[chName]

            self.build_connectivity_graph()
        except:
            print('Failed to update channel library from file. Probably is just half-written.')
            return

        # reset pulse cache
        from . import PulsePrimitives
        PulsePrimitives._memoize.cache.clear()

    def update_from_json(self, chName, chParams):
        # connect objects labeled by strings
        if 'pulseParams' in chParams.keys():
            paramDict = {str(k): v for k, v in chParams['pulseParams'].items()}
            shapeFunName = paramDict.pop('shapeFun', None)
            if shapeFunName:
                paramDict['shapeFun'] = getattr(PulseShapes, shapeFunName)
            self.channelDict[chName].pulseParams = paramDict

        for param in self.specialParams:
            if param in chParams.keys():
                setattr(self.channelDict[chName],
                        param,
                        self.channelDict.get(chParams[param], None)
                        )
        # TODO: how do we follow changes to selected AWG or generator?

        # ignored or specially handled parameters
        ignoreList = self.specialParams + ['pulseParams', 'AWG', 'generator',
                     'x__class__', 'x__module__']
        for paramName in chParams:
            if paramName not in ignoreList:
                setattr(self.channelDict[chName], paramName,
                        chParams[paramName])

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


class ChannelDecoder(json.JSONDecoder):
    def __init__(self, **kwargs):
        super(ChannelDecoder, self).__init__(object_hook=self.dict_to_obj,
                                             **kwargs)

    def dict_to_obj(self, jsonDict):
        import QGL.PulseShapes
        if 'x__class__' in jsonDict or '__class__' in jsonDict:
            #Pop the class and module
            className = jsonDict.pop('x__class__', None)
            if not className:
                className = jsonDict.pop('__class__')
            moduleName = jsonDict.pop('x__module__', None)
            if not moduleName:
                moduleName = jsonDict.pop('__module__')

            __import__(moduleName)

            #Re-encode the strings as ascii (this should go away in Python 3)
            jsonDict = {str(k): v for k, v in jsonDict.items()}

            # instantiate the object
            inst = getattr(sys.modules[moduleName], className)(**jsonDict)

            return inst
        else:
            #Re-encode the strings as ascii (this should go away in Python 3)
            jsonDict = {str(k): v for k, v in jsonDict.items()}
            shapeFun = jsonDict.pop('shapeFun', None)
            if shapeFun:
                jsonDict['shapeFun'] = getattr(QGL.PulseShapes, shapeFun)
            return jsonDict


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

# global channel library
migrator = JSONMigrators.ChannelMigrator(config.channelLibFile)
migrationMsg = migrator.migrate()
for msg in migrationMsg:
    print(msg)

channelLib = ChannelLibrary(libFile=config.channelLibFile)
