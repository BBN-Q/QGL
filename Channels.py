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
import config
import numpy as np

from math import tan,cos,pi

from instruments.AWGs import AWG
from instruments.MicrowaveSources import MicrowaveSource

from traits.api import HasTraits, Str, Float, Instance, DelegatesTo, Property, cached_property, \
                        DictStrAny, Dict, Either, Enum, Bool

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
    ampFactor = Float()
    phaseSkew = Float()
    SSBFreq = Float()
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
    pulseParams = DictStrAny({'length':20e-9, 'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.gaussian, 'buffer':0.0, 'cutoff':2, 'dragScaling':0})

class Measurement(LogicalChannel):
    '''
    A class for measurement channels. 
    Measurments are special because they can be different types:
    autodyne which needs an IQ pair or hetero/homodyne which needs just a marker channel. 
    '''
    measType = Enum('autodyne','homodyne', desc='Type of measurment (autodyne, homodyne)')
    autodyneFreq = Float
    pulseParams = DictStrAny({'length':20e-9, 'amp':1.0, 'shapeFun':PulseShapes.tanh, 'buffer':0.0, 'cutoff':2})


def QubitFactory(name, **kwargs):
    #delayed import to avoid circular imports
    from Libraries import channelLib
    ''' Return a saved qubit channel or create a new one. '''
    if name in channelLib.channelDict and isinstance(channelLib[name], Qubit):
        return channelLib[name]
    else:
        return Qubit(name=name, **kwargs)

def MeasFactory(name, measType='autodyne', **kwargs):
    ''' Return a saved measurment channel or create a new one. '''
    from Libraries import channelLib
    if name in channelLib.channelDict and isinstance(channelLib[name], Measurement):
        return channelLib[name]
    else:
        if measType == 'autodyne':
            return Measurement()

class ChannelLibrary(HasTraits):
    channelDict = Dict(Str, Channel)
    libFile = Str(transient=True)

    def __init__(self, **kwargs):
        super(ChannelLibrary, self).__init__(**kwargs)
        self.load_from_library()

    #Overload [] to allow direct pulling of channel info
    def __getitem__(self, chanName):
        return self.channelDict[chanName]

    @on_trait_change('channelDict.anytrait')
    def write_to_library(self):
        import JSONHelpers
        if self.libFile:
            with open(self.libFile, 'w') as FID:
                json.dump(self, FID, cls=JSONHelpers.LibraryEncoder, indent=2, sort_keys=True)

    def load_from_library(self):
        import JSONHelpers
        if self.libFile:
            try:
                with open(self.libFile, 'r') as FID:
                    tmpLib = json.load(FID, cls=JSONHelpers.ChannelDecoder)
                    if isinstance(tmpLib, ChannelLibrary):
                        for chan in tmpLib.channelDict.values():
                            if isinstance(chan, LogicalChannel):
                                chan.physChan = tmpLib[chan.physChan]
                            elif isinstance(chan, PhysicalQuadratureChannel):
                                chan.gateChan = tmpLib[chan.gateChan]
                        self.channelDict.update(tmpLib.channelDict)

            except IOError:
                print('No channel library found.')
            except ValueError:
                print('Failed to load channel library.')

NewLogicalChannelList = [Qubit, LogicalMarkerChannel, Measurement]
NewPhysicalChannelList = [PhysicalMarkerChannel, PhysicalQuadratureChannel]

if __name__ == '__main__':
    # create a channel params file
    import QGL.Channels
    from Libraries import instrumentLib, channelLib
    import enaml
    from enaml.stdlib.sessions import show_simple_view

    with enaml.imports():
        from ChannelsViews import ChannelLibraryWindow
    show_simple_view(ChannelLibraryWindow(channelLib=channelLib, instrumentLib=instrumentLib))

    # channelLib.channelDict['BBNAPS1-1m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS1-1m1', AWG=instrumentLib['BBNAPS1'], delay=+30e-9)
    # channelLib.channelDict['BBNAPS1-2m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS1-2m1', AWG=instrumentLib['BBNAPS1'], delay=-100e-9)
    # channelLib.channelDict['BBNAPS1-3m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS1-3m1', AWG=instrumentLib['BBNAPS1'], delay=0e-9)
    # channelLib.channelDict['BBNAPS1-4m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS1-4m1', AWG=instrumentLib['BBNAPS1'], delay=0e-9)
    # channelLib.channelDict['BBNAPS1-12'] = QGL.Channels.PhysicalQuadratureChannel(name='BBNAPS1-12', AWG=instrumentLib['BBNAPS1'], generator=instrumentLib['Agilent1'], IChannel='ch1', QChannel='ch2', gateChan=channelLib['BBNAPS1-1m1'], ampFactor=0.8942, phaseSkew=-2.00)
    # channelLib.channelDict['BBNAPS1-34'] = QGL.Channels.PhysicalQuadratureChannel(name='BBNAPS1-34', AWG=instrumentLib['BBNAPS1'], generator=instrumentLib['Autodyne1'], IChannel='ch3', QChannel='ch4', gateChan=channelLib['BBNAPS1-3m1'], ampFactor=0.9374, phaseSkew=5.662)

    # channelLib.channelDict['BBNAPS2-1m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS2-1m1', AWG=instrumentLib['BBNAPS2'], delay=0e-9)
    # channelLib.channelDict['BBNAPS2-2m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS2-2m1', AWG=instrumentLib['BBNAPS2'], delay=-100e-9)
    # channelLib.channelDict['BBNAPS2-3m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS2-3m1', AWG=instrumentLib['BBNAPS2'], delay=+30e-9)
    # channelLib.channelDict['BBNAPS2-4m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS2-4m1', AWG=instrumentLib['BBNAPS2'], delay=0e-9)
    # channelLib.channelDict['BBNAPS2-12'] = QGL.Channels.PhysicalQuadratureChannel(name='BBNAPS2-12', AWG=instrumentLib['BBNAPS2'], generator=instrumentLib['Autodyne2'], IChannel='ch1', QChannel='ch2', gateChan=channelLib['BBNAPS2-1m1'], ampFactor=0.9454, phaseSkew=8.7596)
    # channelLib.channelDict['BBNAPS2-34'] = QGL.Channels.PhysicalQuadratureChannel(name='BBNAPS2-34', AWG=instrumentLib['BBNAPS2'], generator=instrumentLib['Agilent2'], IChannel='ch3', QChannel='ch4', gateChan=channelLib['BBNAPS2-3m1'], ampFactor=1.0045, phaseSkew=2.09)

    # channelLib.channelDict['BBNAPS3-1m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS3-1m1', AWG=instrumentLib['BBNAPS3'], delay=0e-9)
    # channelLib.channelDict['BBNAPS3-2m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS3-2m1', AWG=instrumentLib['BBNAPS3'], delay=-100e-9)
    # channelLib.channelDict['BBNAPS3-3m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS3-3m1', AWG=instrumentLib['BBNAPS3'], delay=+30e-9)
    # channelLib.channelDict['BBNAPS3-4m1'] = QGL.Channels.PhysicalMarkerChannel(name='BBNAPS3-4m1', AWG=instrumentLib['BBNAPS3'], delay=-200e-9)
    # channelLib.channelDict['BBNAPS3-12'] = QGL.Channels.PhysicalQuadratureChannel(name='BBNAPS3-12', AWG=instrumentLib['BBNAPS3'], generator=instrumentLib['Autodyne3'], IChannel='ch1', QChannel='ch2', gateChan=channelLib['BBNAPS3-1m1'], ampFactor=0.9454, phaseSkew=8.7596)
    # channelLib.channelDict['BBNAPS3-34'] = QGL.Channels.PhysicalQuadratureChannel(name='BBNAPS3-34', AWG=instrumentLib['BBNAPS3'], generator=instrumentLib['LabBrick1'], IChannel='ch3', QChannel='ch4', gateChan=channelLib['BBNAPS3-3m1'], ampFactor=1.0045, phaseSkew=2.09)

    # channelLib.channelDict['q1'] = QGL.Channels.Qubit(name='q1',  physChan=channelLib['BBNAPS2-34'], pulseParams={'piAmp':0.7179, 'pi2Amp':0.3604, 'shapeFun':QGL.PulseShapes.drag, 'length':40e-9, 'buffer':2e-9, 'dragScaling':0.88})
    # channelLib.channelDict['q2'] = QGL.Channels.Qubit(name='q2', physChan=channelLib['BBNAPS1-12'], pulseParams={'piAmp':0.6583, 'pi2Amp':0.3153, 'shapeFun':QGL.PulseShapes.drag, 'length':40e-9, 'buffer':2e-9, 'dragScaling':2.65})
    # channelLib.channelDict['q3'] = QGL.Channels.Qubit(name='q3', physChan=channelLib['BBNAPS3-34'], pulseParams={'piAmp':0.6583, 'pi2Amp':0.3153, 'shapeFun':QGL.PulseShapes.drag, 'length':40e-9, 'buffer':2e-9, 'dragScaling':2.65})
    # channelLib.channelDict['digitizerTrig'] = QGL.Channels.LogicalMarkerChannel(name='digitizerTrig', physChan=channelLib['BBNAPS1-2m1'], pulseParams={'length':1e-9, 'amp':1.0, 'shapeFun':QGL.PulseShapes.square})
    # channelLib.channelDict['slaveTrig'] = QGL.Channels.LogicalMarkerChannel(name='slaveTrig', physChan=channelLib['BBNAPS1-4m1'], pulseParams={'length':1e-9, 'amp':1.0, 'shapeFun':QGL.PulseShapes.square})

    # channelLib.channelDict['M-q1'] = QGL.Channels.Measurement(name='M-q1', measType='autodyne', autodyneFreq=10e6, physChan=channelLib['BBNAPS2-12'], pulseParams={'amp':0.5, 'shapeFun':PulseShapes.tanh, 'length':3.33333e-6, 'buffer':2e-9})
    # channelLib.channelDict['M-q2'] = QGL.Channels.Measurement(name='M-q2', measType='autodyne', autodyneFreq=20e6, physChan=channelLib['BBNAPS1-34'], pulseParams={'amp':0.5, 'shapeFun':PulseShapes.tanh, 'length':3.33333e-6, 'buffer':2e-9})
    # channelLib.channelDict['M-q3'] = QGL.Channels.Measurement(name='M-q3', measType='autodyne', autodyneFreq=10e6, physChan=channelLib['BBNAPS3-12'], pulseParams={'amp':0.5, 'shapeFun':PulseShapes.tanh, 'length':3.33333e-6, 'buffer':2e-9})
    
    # channelLib.write_to_library()

    # ChannelDict['q1q2'] = Qubit(name='q1q2', physChan='BBNAPS1-34', pulseParams={'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.drag, 'pulseLength':40e-9, 'buffer':2e-9, 'dragScaling':1})
    # ChannelDict['M-q1q2'] = Measurement(name='M-q1q2', measType='autodyne', physChan='BBNAPS2-34', trigChan='digitizerTrig', pulseParams={'amp':1.0, 'shapeFun':PulseShapes.tanh, 'length':1.6e-6, 'buffer':2e-9})
    


    
    