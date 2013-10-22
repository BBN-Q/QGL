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

from traits.api import HasTraits, Str, Float, Instance, DelegatesTo, Property, cached_property, \
                        DictStrAny, Dict, Either, Enum, Bool, on_trait_change, Any

import FileWatcher

class Channel(HasTraits):
    '''
    Every channel has a name and some printers.
    '''
    name = Str()
    enabled = Bool(True)

    def __repr__(self):
        return "{0}: {1}".format(self.__class__.__name__, self.name)

    def __str__(self):
        return json.dumps(self, sort_keys=True, indent=2, default=json_serializer)

class PhysicalChannel(Channel):
    '''
    The main class for actual AWG channels.
    '''
    AWG = Instance(AWG)
    generator = Instance(MicrowaveSource)
    samplingRate = DelegatesTo('AWG')
    delay = Float()

    def __init__(self, **kwargs):
        super(PhysicalChannel, self).__init__(**kwargs)
        if self.AWG is None:
            self.AWG = AWG()


class LogicalChannel(Channel):
    '''
    The main class from which we will generate sequences. 
    At some point it needs to be assigned to a physical channel.
    '''
    #During initilization we may just have a string reference to the channel
    physChan = Either(Str, Instance(PhysicalChannel))
    AWG = DelegatesTo('physChan')

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
    gateChan = Either(Str, Instance(PhysicalMarkerChannel))
    ampFactor = Float(1.0)
    phaseSkew = Float(0.0)
    SSBFreq = Float(0.0)
    correctionT = Property(depends_on=['ampFactor', 'phaseSkew'])

    @cached_property
    def _get_correctionT(self):
        return np.array([[self.ampFactor, self.ampFactor*tan(self.phaseSkew*pi/180)], [0, 1/cos(self.phaseSkew*pi/180)]])
                
class LogicalMarkerChannel(LogicalChannel):
    '''
    A class for digital channels for gating sources or triggering other things.
    '''
    pulseParams = DictStrAny({'shapeFun': PulseShapes.square, 'length':100e-9})

class Qubit(LogicalChannel):
    '''
    The main class for generating qubit pulses.  Effectively a logical "QuadratureChannel".
    '''
    pulseParams = DictStrAny({'length':20e-9, 'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.gaussian, 'buffer':0.0, 'cutoff':2, 'dragScaling':0, 'sigma':5e-9})

    def __init__(self, **kwargs):
        super(Qubit, self).__init__(**kwargs)
        if self.physChan is None:
            self.physChan = PhysicalQuadratureChannel(name=kwargs['name']+'-phys')

class Measurement(LogicalChannel):
    '''
    A class for measurement channels. 
    Measurments are special because they can be different types:
    autodyne which needs an IQ pair or hetero/homodyne which needs just a marker channel. 
    '''
    measType = Enum('autodyne','homodyne', desc='Type of measurment (autodyne, homodyne)')
    autodyneFreq = Float
    pulseParams = DictStrAny({'length':100e-9, 'amp':1.0, 'shapeFun':PulseShapes.tanh, 'buffer':0.0, 'cutoff':2, 'sigma':1e-9})


def QubitFactory(name, **kwargs):
    ''' Return a saved qubit channel or create a new one. '''
    if Compiler.channelLib and name in Compiler.channelLib.channelDict and isinstance(Compiler.channelLib[name], Qubit):
        return Compiler.channelLib[name]
    else:
        return Qubit(name=name, **kwargs)

def MeasFactory(name, measType='autodyne', **kwargs):
    ''' Return a saved measurment channel or create a new one. '''
    if Compiler.channelLib and name in Compiler.channelLib.channelDict and isinstance(Compiler.channelLib[name], Measurement):
        return Compiler.channelLib[name]
    else:
        return Measurement(measType = measType)

class ChannelLibrary(HasTraits):
    channelDict = Dict(Str, Channel)
    libFile = Str(transient=True)
    fileWatcher = Any(None, transient=True)

    def __init__(self, **kwargs):
        super(ChannelLibrary, self).__init__(**kwargs)
        self.load_from_library()
        if self.libFile:
            self.fileWatcher = FileWatcher.LibraryFileWatcher(self.libFile, self.update_from_file)

    #Overload [] to allow direct pulling of channel info
    def __getitem__(self, chanName):
        return self.channelDict[chanName]

    @on_trait_change('channelDict.anytrait')
    def write_to_library(self):
        import JSONHelpers
        if self.libFile:
            #Pause the file watcher to stop cicular updating insanity
            if self.fileWatcher:
                self.fileWatcher.pause()
            with open(self.libFile, 'w') as FID:
                json.dump(self, FID, cls=JSONHelpers.LibraryEncoder, indent=2, sort_keys=True)
            if self.fileWatcher:
                self.fileWatcher.resume()

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
    from enaml.stdlib.sessions import show_simple_view

    with enaml.imports():
        from ChannelsViews import ChannelLibraryWindow
    show_simple_view(ChannelLibraryWindow(channelLib=Compiler.channelLib, instrumentLib=instrumentLib))



    
    