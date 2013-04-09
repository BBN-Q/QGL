'''
Channels is where we store information for mapping virtual (qubit) channel to real channels.

Created on Jan 19, 2012

Original Author: Colm Ryan

'''

import sys
import json
import PulseShapes
import config
import numpy as np

from math import tan,cos,pi
import config

from warnings import warn

from instruments.InstrumentManager import InstrumentLibrary
from instruments.AWGs import AWG
from instruments.MicrowaveSources import MicrowaveSource

from traits.api import HasTraits, Str, Float, Instance, DelegatesTo, Property, cached_property, \
                        DictStrAny, Dict

class Channel(HasTraits):
    '''
    Every channel has a name and some printers.
    '''
    name = Str()

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

    def _set_AWG(self, AWG):
        print('got here')

    def _set_generator(self, gen):
        print('got here')

class LogicalChannel(Channel):
    '''
    The main class from which we will generate sequences. 
    At some point it needs to be assigned to a physical channel.
    '''
    physicalChannel = Instance(PhysicalChannel)

class PhysicalMarkerChannel(PhysicalChannel):
    '''
    An digital output channel on an AWG.
    '''
    markerCoupling = DelegatesTo('AWG')
        
class PhysicalQuadratureChannel(PhysicalChannel):
    '''
    Something used to implement a standard qubit channel with two analog channels and a microwave gating channel.
    '''
    IChannel = Str()
    QChannel = Str()
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
    def __init__(self, name=None, measType='autodyne', physicalChannel=None, trigChan=None, pulseParams=None):
        super(Measurement, self).__init__(name=name, physicalChannel=physicalChannel)
        self.measType = measType
        self._trigChan = trigChan
        defaultPulseParams = {'length':20e-9, 'amp':1.0, 'shapeFun':PulseShapes.tanh, 'buffer':0.0, 'cutoff':2}
        defaultPulseParams.update(pulseParams)
        self.pulseParams = defaultPulseParams
        
    @property
    def trigChan(self):
        if self._trigChan:
            return ChannelDict[self._trigChan]
        else:
            return LogicalMarkerChannel()

def QubitFactory(name, **kwargs):
    ''' Return a saved qubit channel or create a new one. '''
    if name in ChannelDict and isinstance(ChannelDict[name], Qubit):
        return ChannelDict[name]
    else:
        return Qubit(name, **kwargs)

def MeasFactory(name, measType='autodyne', **kwargs):
    ''' Return a saved measurment channel or create a new one. '''
    if name in ChannelDict and isinstance(ChannelDict[name], Measurement):
        return ChannelDict[name]
    else:
        if measType == 'autodyne':
            return Measurment

def update_channel_info(fileName=None):
    '''
    Helper function to load a channel info file into channel and channelInfo objects
    '''
    if fileName is None:
        fileName = config.channelParamsFile
    with open(fileName,'r') as FID:
        try:
            ChannelDict.update(json.load(FID, object_hook=json_deserializer))
        except Exception, e:
            warn('Unable to update channel information.')

class ChannelLibrary(HasTraits):
    channelDict = Dict(Str, Channel)
    libFile = Str(transient=True)
    instrLib = Instance(InstrumentLibrary)

    #Overload [] to allow direct pulling of channel info
    def __getitem__(self, chanName):
        return self.channelDict[chanName]

    def write_to_library(self):
        import JSONHelpers
        if self.libFile:
            with open(self.libFile, 'w') as FID:
                json.dump(self, FID, cls=JSONHelpers.LibraryEncoder, indent=2, sort_keys=True)

        pass

    def load_from_library(self):
        import JSONHelpers
        if self.libFile:
            try:
                with open(self.libFile, 'r') as FID:
                    tmpLib = json.load(FID, cls=JSONHelpers.ChannelDecoder, instrLib=self.instrLib)
                    import pdb; pdb.set_trace()
                    if isinstance(tmpLib, ChannelLibrary):
                        self.channelDict.update(tmpLib.channelDict)
            except IOError:
                print('No sweep library found.')


if __name__ == '__main__':
    # # create a channel params file
    import QGL.Channels
    instrLib = InstrumentLibrary(libFile=config.instrumentLibFile)
    channelLib = QGL.Channels.ChannelLibrary(instrLib=instrLib)

    channelLib.channelDict['BBNAPS1-12'] = QGL.Channels.PhysicalQuadratureChannel(name='BBNAPS1-12', AWG=instrLib.instrDict['BBNAPS1'], generator=instrLib.instrDict['Agilent1'], IChannel='ch1', QChannel='ch2', ampFactor=1.0252, phaseSkew=-4.97)
    channelLib.channelDict['BBNAPS1-34'] = QGL.Channels.PhysicalQuadratureChannel(name='BBNAPS1-34', AWG=instrLib.instrDict['BBNAPS1'], generator=instrLib.instrDict['Agilent2'], IChannel='ch3', QChannel='ch4', ampFactor=1, phaseSkew=0, SSBFreq=31.9e6)
    

    # ChannelDict['BBNAPS2-12'] = PhysicalQuadratureChannel(name='BBNAPS2-12', AWG='BBNAPS2', generator='Agilent2', IChannel='ch1', QChannel='ch2', ampFactor=1, phaseSkew=0)
    # ChannelDict['BBNAPS2-34'] = PhysicalQuadratureChannel(name='BBNAPS2-34', AWG='BBNAPS2', generator='Agilent3', IChannel='ch3', QChannel='ch4', ampFactor=1, phaseSkew=0, SSBFreq=-20e6)
    # ChannelDict['BBNAPS1-1m1'] = PhysicalMarkerChannel(name='BBNAPS1-1m1', AWG='BBNAPS1')
    # ChannelDict['BBNAPS1-2m1'] = PhysicalMarkerChannel(name='BBNAPS1-2m1', AWG='BBNAPS1', delay=-50e-9)
    # ChannelDict['BBNAPS1-3m1'] = PhysicalMarkerChannel(name='BBNAPS1-3m1', AWG='BBNAPS1')
    # ChannelDict['BBNAPS1-4m1'] = PhysicalMarkerChannel(name='BBNAPS1-4m1', AWG='BBNAPS1', delay=-200e-9)
    # ChannelDict['BBNAPS2-1m1'] = PhysicalMarkerChannel(name='BBNAPS2-1m1', AWG='BBNAPS2')
    # ChannelDict['BBNAPS2-2m1'] = PhysicalMarkerChannel(name='BBNAPS2-2m1', AWG='BBNAPS2')
    # ChannelDict['BBNAPS2-3m1'] = PhysicalMarkerChannel(name='BBNAPS2-3m1', AWG='BBNAPS2')
    # ChannelDict['BBNAPS2-4m1'] = PhysicalMarkerChannel(name='BBNAPS2-4m1', AWG='BBNAPS2')

    # ChannelDict['q1'] = Qubit(name='q1',  physicalChannel='BBNAPS1-12', pulseParams={'piAmp':0.7313, 'pi2Amp':0.3648, 'shapeFun':PulseShapes.drag, 'length':26.67e-9, 'buffer':2e-9, 'dragScaling':0.3})
    # ChannelDict['q2'] = Qubit(name='q2', physicalChannel='BBNAPS2-12', pulseParams={'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.drag, 'length':40e-9, 'buffer':2e-9, 'dragScaling':1})
    # ChannelDict['q1q2'] = Qubit(name='q1q2', physicalChannel='BBNAPS1-34', pulseParams={'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.drag, 'pulseLength':40e-9, 'buffer':2e-9, 'dragScaling':1})
    # ChannelDict['M-q1'] = Measurement(name='M-q1', measType='autodyne', physicalChannel='BBNAPS1-34', trigChan='digitizerTrig', pulseParams={'amp':1.0, 'shapeFun':PulseShapes.tanh, 'length':1.6e-6, 'buffer':2e-9})
    # ChannelDict['M-q2'] = Measurement(name='M-q2', measType='autodyne', physicalChannel='BBNAPS2-34', trigChan='digitizerTrig', pulseParams={'amp':1.0, 'shapeFun':PulseShapes.tanh, 'length':1.6e-6, 'buffer':2e-9})
    # ChannelDict['M-q1q2'] = Measurement(name='M-q1q2', measType='autodyne', physicalChannel='BBNAPS2-34', trigChan='digitizerTrig', pulseParams={'amp':1.0, 'shapeFun':PulseShapes.tanh, 'length':1.6e-6, 'buffer':2e-9})
    
    # ChannelDict['digitizerTrig'] = LogicalMarkerChannel(name='digitizerTrig', physicalChannel='BBNAPS1-2m1', pulseParams={'length':40e-9, 'amp':1.0, 'shapeFun':PulseShapes.square})
    # ChannelDict['slaveTrig'] = LogicalMarkerChannel(name='slaveTrig', physicalChannel='BBNAPS1-4m1', pulseParams={'length':1e-9, 'amp':1.0, 'shapeFun':PulseShapes.square})
    # save_channel_info()


    
    