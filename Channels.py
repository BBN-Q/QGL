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

from math import tan,cos,pi

from instruments.AWGs import AWG
from instruments.MicrowaveSources import MicrowaveSource
from DictManager import DictManager

from atom.api import Atom, Str, Unicode, Float, Instance, Property, cached_property, \
                        Dict, Enum, Bool, Typed, observe

import FileWatcher

from copy import deepcopy

class Channel(Atom):
    '''
    Every channel has a label and some printers.
    '''
    label = Str()
    enabled = Bool(True)

    def __repr__(self):
        return "{0}: {1}".format(self.__class__.__name__, self.label)

    def __str__(self):
        return json.dumps(self, sort_keys=True, indent=2, default=json_serializer)

    def json_encode(self):
        jsonDict = self.__getstate__()

        #Strip out pass-through properties
        for k,m in self.members().items():
            if isinstance(m, Property):
                del jsonDict[k]

        #Turn instruments back into unicode labels
        for member in ["AWG", "generator", "physChan", "gateChan"]:
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
    AWG = Property()

    def _get_AWG(self):
        return self.physChan.AWG

class PhysicalMarkerChannel(PhysicalChannel):
    '''
    An digital output channel on an AWG.
    '''
        
class PhysicalQuadratureChannel(PhysicalChannel):
    '''
    Something used to implement a standard qubit channel with two analog channels and a microwave gating channel.
    '''
    IChannel = Str()
    QChannel = Str()
    #During initilization we may just have a string reference to the channel
    gateChan = Instance((unicode, PhysicalMarkerChannel))
    ampFactor = Float(1.0)
    phaseSkew = Float(0.0)
    SSBFreq = Float(0.0)

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
    pulseParams = Dict(default={'length':20e-9, 'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.gaussian, 'buffer':0.0, 'cutoff':2, 'dragScaling':0, 'sigma':5e-9})

    def __init__(self, **kwargs):
        super(Qubit, self).__init__(**kwargs)
        if self.physChan is None:
            self.physChan = PhysicalQuadratureChannel(label=kwargs['label']+'-phys')

class Measurement(LogicalChannel):
    '''
    A class for measurement channels. 
    Measurments are special because they can be different types:
    autodyne which needs an IQ pair or hetero/homodyne which needs just a marker channel. 
    '''
    measType = Enum('autodyne','homodyne').tag(desc='Type of measurment (autodyne, homodyne)')
    autodyneFreq = Float()
    pulseParams = Dict(default={'length':100e-9, 'amp':1.0, 'shapeFun':PulseShapes.tanh, 'buffer':0.0, 'cutoff':2, 'sigma':1e-9})


def QubitFactory(label, **kwargs):
    ''' Return a saved qubit channel or create a new one. '''
    if Compiler.channelLib and label in Compiler.channelLib.channelDict and isinstance(Compiler.channelLib[label], Qubit):
        return Compiler.channelLib[label]
    else:
        return Qubit(label=label, **kwargs)

def MeasFactory(label, measType='autodyne', **kwargs):
    ''' Return a saved measurment channel or create a new one. '''
    if Compiler.channelLib and label in Compiler.channelLib.channelDict and isinstance(Compiler.channelLib[label], Measurement):
        return Compiler.channelLib[label]
    else:
        return Measurement(label=label, measType = measType, **kwargs)

class ChannelLibrary(Atom):
    # channelDict = Dict(Str, Channel)
    channelDict = Typed(dict)
    logicalChannelManager = Typed(DictManager)
    physicalChannelManager = Typed(DictManager)
    libFile = Str()
    fileWatcher = Typed(FileWatcher.LibraryFileWatcher)

    def __init__(self, channelDict={}, **kwargs):
        super(ChannelLibrary, self).__init__(channelDict=channelDict, **kwargs)
        self.load_from_library()
        self.logicalChannelManager = DictManager(itemDict=self.channelDict,
                                                 displayFilter=lambda x : isinstance(x, LogicalChannel),
                                                 possibleItems=NewLogicalChannelList)
        self.physicalChannelManager = DictManager(itemDict=self.channelDict,
                                                  displayFilter=lambda x : isinstance(x, PhysicalChannel),
                                                  possibleItems=NewPhysicalChannelList)

        if self.libFile:
            self.fileWatcher = FileWatcher.LibraryFileWatcher(self.libFile, self.update_from_file)

    #Overload [] to allow direct pulling of channel info
    def __getitem__(self, chanLabel):
        return self.channelDict[chanLabel]

    def write_to_file(self):
        import JSONHelpers
        if self.libFile:
            #Pause the file watcher to stop cicular updating insanity
            if self.fileWatcher:
                self.fileWatcher.pause()
            with open(self.libFile, 'w') as FID:
                json.dump(self, FID, cls=JSONHelpers.LibraryEncoder, indent=2, sort_keys=True)
            if self.fileWatcher:
                self.fileWatcher.resume()

    def json_encode(self, matlabCompatible=False):
        return {"channelDict":self.channelDict}

    def load_from_library(self):
        import JSONHelpers
        if self.libFile:
            try:
                with open(self.libFile, 'r') as FID:
                    tmpLib = json.load(FID, cls=JSONHelpers.ChannelDecoder)
                    if isinstance(tmpLib, ChannelLibrary):
                        for chan in tmpLib.channelDict.values():
                            if isinstance(chan, LogicalChannel):
                                chan.physChan = tmpLib[chan.physChan] if chan.physChan in tmpLib.channelDict else None
                            elif isinstance(chan, PhysicalQuadratureChannel):
                                chan.gateChan = tmpLib[chan.gateChan] if chan.gateChan in tmpLib.channelDict else None
                        self.channelDict.update(tmpLib.channelDict)

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
        updateList = ['ampFactor', 'phaseSkew', 'SSBFreq', 'delay', 'pulseParams']
        if self.libFile:
            with open(self.libFile, 'r') as FID:
                try:
                    allParams = json.load(FID)['channelDict']
                except ValueError:
                    print('Failed to update channel library from file. Probably is just half-written.')
                    return
                for chName, chParams in allParams.items():
                    if chName in self.channelDict:
                        for paramName in updateList:
                            if paramName in chParams:
                                #Deal with unicode/string difference
                                if paramName == 'pulseParams':
                                    paramDict = {k.encode('ascii'):v for k,v in chParams['pulseParams'].items()}
                                    shapeFunName = paramDict.pop('shapeFun', None)
                                    if shapeFunName:
                                        paramDict['shapeFun'] = getattr(PulseShapes, shapeFunName)
                                    setattr(self.channelDict[chName], 'pulseParams', paramDict)
                                else:
                                    setattr(self.channelDict[chName], paramName, chParams[paramName])


NewLogicalChannelList = [Qubit, LogicalMarkerChannel, Measurement]
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
