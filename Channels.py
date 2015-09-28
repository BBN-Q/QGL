'''
Channels is where we store information for mapping virtual (qubit) channel to real channels.

Created on Jan 19, 2012

Original Author: Colm Ryan

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

import json
import PulseShapes
import Compiler
import numpy as np
import networkx as nx

from math import tan,cos,pi

from instruments.AWGs import AWG
from instruments.MicrowaveSources import MicrowaveSource
from DictManager import DictManager

from atom.api import Atom, Str, Unicode, Float, Instance, Property, cached_property, \
                        Dict, Enum, Bool, Typed, observe, Int

import FileWatcher

import importlib

from copy import deepcopy

class Channel(Atom):
    '''
    Every channel has a label and some printers.
    '''
    label = Str()
    enabled = Bool(True)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "{0}: {1}".format(self.__class__.__name__, self.label)

    def json_encode(self):
        jsonDict = self.__getstate__()

        #Strip out pass-through properties
        for k,m in self.members().items():
            if isinstance(m, Property):
                del jsonDict[k]

        #Turn instruments back into unicode labels
        for member in ["AWG", "generator", "physChan", "gateChan", "trigChan", "source", "target"]:
            if member in jsonDict:
                obj = jsonDict.pop(member)
                if obj:
                    jsonDict[member] = obj.label

        #We want the name of shape functions
        if "pulseParams" in jsonDict:
            pulseParams = deepcopy(jsonDict.pop("pulseParams"))
            if "shapeFun" in pulseParams:
                pulseParams["shapeFun"] = pulseParams["shapeFun"].__name__
            jsonDict["pulseParams"] = pulseParams

        return jsonDict

class PhysicalChannel(Channel):
    '''
    The main class for actual AWG channels.
    '''
    AWG = Typed(AWG)
    generator = Typed(MicrowaveSource)
    samplingRate = Property()
    delay = Float()

    def _get_samplingRate(self):
        return self.AWG.samplingRate

    def _default_AWG(self):
        return AWG()


class LogicalChannel(Channel):
    '''
    The main class from which we will generate sequences.
    At some point it needs to be assigned to a physical channel.
    '''
    #During initilization we may just have a string reference to the channel
    physChan = Instance((unicode,PhysicalChannel))

    def __init__(self, **kwargs):
        super(LogicalChannel, self).__init__(**kwargs)
        if self.physChan is None:
            self.physChan = PhysicalChannel(label=kwargs['label']+'-phys')

class PhysicalMarkerChannel(PhysicalChannel):
    '''
    An digital output channel on an AWG.
    '''
    gateBuffer = Float(0.0).tag(desc="How much extra time should be added onto the beginning of a gating pulse")
    gateMinWidth = Float(0.0).tag(desc="The minimum marker pulse width")

class PhysicalQuadratureChannel(PhysicalChannel):
    '''
    Something used to implement a standard qubit channel with two analog channels and a microwave gating channel.
    '''
    IChannel = Str()
    QChannel = Str()
    #During initilization we may just have a string reference to the channel
    ampFactor = Float(1.0)
    phaseSkew = Float(0.0)

    @cached_property
    def correctionT(self):
        return np.array([[self.ampFactor, self.ampFactor*tan(self.phaseSkew*pi/180)], [0, 1/cos(self.phaseSkew*pi/180)]])

    @observe('ampFactor', 'phaseSkew')
    def _reset_correctionT(self, change):
        if change['type'] == 'update':
            self.get_member('correctionT').reset(self)

class LogicalMarkerChannel(LogicalChannel):
    '''
    A class for digital channels for gating sources or triggering other things.
    '''
    pulseParams = Dict(default={'shapeFun': PulseShapes.square, 'length':100e-9})

class Qubit(LogicalChannel):
    '''
    The main class for generating qubit pulses.  Effectively a logical "QuadratureChannel".
    '''
    pulseParams = Dict(default={'length':20e-9, 'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.gaussian, 'cutoff':2, 'dragScaling':0, 'sigma':5e-9})
    gateChan = Instance((unicode, LogicalMarkerChannel))
    frequency = Float(0.0).tag(desc='modulation frequency of the channel (can be positive or negative)')

    def __init__(self, **kwargs):
        super(Qubit, self).__init__(**kwargs)
        if self.gateChan is None:
            self.gateChan = LogicalMarkerChannel(label=kwargs['label']+'-gate')

class Measurement(LogicalChannel):
    '''
    A class for measurement channels.
    Measurements are special because they can be different types:
    autodyne which needs an IQ pair or hetero/homodyne which needs just a marker channel.
    '''
    measType = Enum('autodyne','homodyne').tag(desc='Type of measurment (autodyne, homodyne)')
    autodyneFreq = Float(0.0).tag(desc='use to bake the modulation into the pulse, so that it has constant phase')
    frequency = Float(0.0).tag(desc='use frequency to asssociate modulation with the channel')
    pulseParams = Dict(default={'length':100e-9, 'amp':1.0, 'shapeFun':PulseShapes.tanh, 'cutoff':2, 'sigma':1e-9})
    gateChan = Instance((unicode, LogicalMarkerChannel))
    trigChan = Instance((unicode, LogicalMarkerChannel))
    
    def __init__(self, **kwargs):
        super(Measurement, self).__init__(**kwargs)
        if self.gateChan is None:
            self.gateChan = LogicalMarkerChannel(label=kwargs['label']+'-gate')
        if self.trigChan is None:
            self.trigChan = LogicalMarkerChannel(label='digitizerTrig')

class Edge(LogicalChannel):
    '''
    Defines an arc/directed edge between qubit vertices. If a device supports bi-directional
    connectivity, that is represented with two independent Edges.

    An Edge is also effectively an abstract channel, so it carries the same properties as a 
    Qubit channel.
    '''
    # allow unicode in source and target so that we can store a label or an object
    source = Instance((unicode, Qubit))
    target = Instance((unicode, Qubit))
    pulseParams = Dict(default={'length':20e-9, 'amp':1.0, 'phase':0.0, 'shapeFun':PulseShapes.gaussian, 'cutoff':2, 'dragScaling':0, 'sigma':5e-9, 'riseFall': 20e-9})
    gateChan = Instance((unicode, LogicalMarkerChannel))
    frequency = Float(0.0).tag(desc='modulation frequency of the channel (can be positive or negative)')

    def __init__(self, **kwargs):
        super(Edge, self).__init__(**kwargs)
        if self.gateChan is None:
            self.gateChan = LogicalMarkerChannel(label=kwargs['label']+'-gate')

    def isforward(self, source, target):
        ''' Test whether (source, target) matches the directionality of the edge. '''
        nodes = (self.source, self.target)
        if (source not in nodes) or (target not in nodes):
            raise ValueError('One of {0} is not a node in the edge'.format((source, target)))
        if (self.source, self.target) == (source, target):
            return True
        else:
            return False

def QubitFactory(label, **kwargs):
    ''' Return a saved qubit channel or create a new one. '''
    if Compiler.channelLib and label in Compiler.channelLib and isinstance(Compiler.channelLib[label], Qubit):
        return Compiler.channelLib[label]
    else:
        return Qubit(label=label, **kwargs)

def MeasFactory(label, measType='autodyne', **kwargs):
    ''' Return a saved measurement channel or create a new one. '''
    if Compiler.channelLib and label in Compiler.channelLib and isinstance(Compiler.channelLib[label], Measurement):
        return Compiler.channelLib[label]
    else:
        return Measurement(label=label, measType=measType, **kwargs)

def EdgeFactory(source, target):
    if not Compiler.channelLib:
        raise ValueError('Connectivity graph not found')
    if Compiler.channelLib.connectivityG.has_edge(source, target):
        return Compiler.channelLib.connectivityG[source][target]['channel']
    elif Compiler.channelLib.connectivityG.has_edge(target, source):
        return Compiler.channelLib.connectivityG[target][source]['channel']
    else:
        raise ValueError('Edge {0} not found in connectivity graph'.format((source, target)))

class ChannelLibrary(Atom):
    # channelDict = Dict(Str, Channel)
    channelDict = Typed(dict)
    connectivityG = Typed(nx.DiGraph)
    logicalChannelManager = Typed(DictManager)
    physicalChannelManager = Typed(DictManager)
    libFile = Str()
    fileWatcher = Typed(FileWatcher.LibraryFileWatcher)
    version = Int(0)

    def __init__(self, channelDict={}, **kwargs):
        super(ChannelLibrary, self).__init__(channelDict=channelDict, **kwargs)
        self.connectivityG = nx.DiGraph()
        self.load_from_library()
        self.logicalChannelManager = DictManager(itemDict=self.channelDict,
                                                 displayFilter=lambda x : isinstance(x, LogicalChannel),
                                                 possibleItems=NewLogicalChannelList)
        self.physicalChannelManager = DictManager(itemDict=self.channelDict,
                                                  displayFilter=lambda x : isinstance(x, PhysicalChannel),
                                                  possibleItems=NewPhysicalChannelList)

        if self.libFile:
            self.fileWatcher = FileWatcher.LibraryFileWatcher(self.libFile, self.update_from_file)

    #Dictionary methods
    def __getitem__(self, key):
        return self.channelDict[key]

    def __setitem__(self, key, value):
        self.channelDict[key] = value

    def __delitem__(self, key):
        del self.channelDict[key]

    # def __len__(self):
    #     return len(self.channelDict)

    # def __iter__(self):
    #     for x in self.channelDict:
    #         yield x

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
            if isinstance(chan, Qubit) and chan not in self.connectivityG:
                self.connectivityG.add_node(chan)
        for chan in self.channelDict.values():
            if isinstance(chan, Edge):
                self.connectivityG.add_edge(chan.source, chan.target)
                self.connectivityG[chan.source][chan.target]['channel'] = chan

    def write_to_file(self,fileName=None):
        import JSONHelpers
        libFileName = fileName if fileName != None else self.libFile
        if libFileName:
            #Pause the file watcher to stop cicular updating insanity
            if self.fileWatcher:
                self.fileWatcher.pause()
            with open(libFileName, 'w') as FID:
                json.dump(self, FID, cls=JSONHelpers.LibraryEncoder, indent=2, sort_keys=True)
            if self.fileWatcher:
                self.fileWatcher.resume()

    def json_encode(self, matlabCompatible=False):
        return {
            "channelDict": self.channelDict,
            "version": self.version
        }

    def load_from_library(self):
        import JSONHelpers
        if self.libFile:
            try:
                with open(self.libFile, 'r') as FID:
                    tmpLib = json.load(FID, cls=JSONHelpers.ChannelDecoder)
                    if not isinstance(tmpLib, ChannelLibrary):
                        raise ValueError('Failed to load channel library')

                    # connect objects labeled by strings
                    for chan in tmpLib.channelDict.values():
                        if hasattr(chan, 'physChan'):
                            chan.physChan = tmpLib[chan.physChan] if chan.physChan in tmpLib.channelDict else None
                        if hasattr(chan, 'gateChan'):
                            chan.gateChan = tmpLib[chan.gateChan] if chan.gateChan in tmpLib.channelDict else None
                        if hasattr(chan, 'trigChan'):
                            chan.trigChan = tmpLib[chan.trigChan] if chan.trigChan in tmpLib.channelDict else None
                        if hasattr(chan, 'source'):
                            chan.source = tmpLib[chan.source] if chan.source in tmpLib.channelDict else None
                        if hasattr(chan, 'target'):
                            chan.target = tmpLib[chan.target] if chan.target in tmpLib.channelDict else None
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

        if self.libFile:
            with open(self.libFile, 'r') as FID:
                try:
                    allParams = json.load(FID)['channelDict']
                except ValueError:
                    print('Failed to update channel library from file. Probably is just half-written.')
                    return

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
                        self.channelDict[chName]  = cls()
                        self.update_from_json(chName, chParams)

                # remove
                for chName in self.channelDict.keys():
                    if chName not in allParams:
                        del self.channelDict[chName]

                self.build_connectivity_graph()

    def update_from_json(self,chName, chParams):
        # ignored or specially handled parameters
        ignoreList = ['pulseParams', 'physChan', 'gateChan', 'trigChan', 'source', 'target', 'AWG', 'generator', 'x__class__', 'x__module__']
        if 'pulseParams' in chParams.keys():
            paramDict = {k.encode('ascii'):v for k,v in chParams['pulseParams'].items()}
            shapeFunName = paramDict.pop('shapeFun', None)
            if shapeFunName:
                paramDict['shapeFun'] = getattr(PulseShapes, shapeFunName)
            self.channelDict[chName].pulseParams = paramDict
        if 'physChan' in chParams.keys():
            self.channelDict[chName].physChan = self.channelDict[chParams['physChan']] if chParams['physChan'] in self.channelDict else None
        if 'gateChan' in chParams.keys():
            self.channelDict[chName].gateChan = self.channelDict[chParams['gateChan']] if chParams['gateChan'] in self.channelDict else None
        if 'trigChan' in chParams.keys():
            self.channelDict[chName].trigChan = self.channelDict[chParams['trigChan']] if chParams['trigChan'] in self.channelDict else None
        if 'source' in chParams.keys():
            self.channelDict[chName].source = self.channelDict[chParams['source']] if chParams['source'] in self.channelDict else None
        if 'target' in chParams.keys():
            self.channelDict[chName].target = self.channelDict[chParams['target']] if chParams['target'] in self.channelDict else None
        # TODO: how do we follow changes to selected AWG or generator?

        for paramName in chParams:
            if paramName not in ignoreList:
                setattr(self.channelDict[chName], paramName, chParams[paramName])

    def on_awg_change(self, oldName, newName):
        print "Change AWG", oldName, newName
        for chName in self.channelDict:
            if (isinstance(self.channelDict[chName], PhysicalMarkerChannel) or
               isinstance(self.channelDict[chName], PhysicalQuadratureChannel)):
                awgName, awgChannel = chName.rsplit('-',1)
                if awgName == oldName:
                    newLabel = "{0}-{1}".format(newName,awgChannel)
                    print "Changing {0} to {1}".format(chName, newLabel)
                    self.physicalChannelManager.name_changed(chName, newLabel)


NewLogicalChannelList = [Qubit, Edge, LogicalMarkerChannel, Measurement]
NewPhysicalChannelList = [PhysicalMarkerChannel, PhysicalQuadratureChannel]

if __name__ == '__main__':
    import os.path
    import sys
    #Load the libraries
    execfile(os.path.join(os.path.dirname(sys.argv[0]), os.path.pardir, 'startup.py'))

    # create a channel params file
    import enaml
    from enaml.qt.qt_application import QtApplication

    with enaml.imports():
        from ChannelsViews import ChannelLibraryWindow

    app = QtApplication()
    view = ChannelLibraryWindow(channelLib=Compiler.channelLib, instrumentLib=instrumentLib)
    view.show()

    app.start()
