'''
Channels is where we store information for mapping virtual (qubit) channel to real channels.

Created on Jan 19, 2012

@author: cryan


'''

import sys
import json
import PulseShapes
import config
import numpy as np

from operator import itemgetter
from math import tan,cos,pi
import config

from warnings import warn
from types import FunctionType

ChannelDict = {}

class Channel(object):
    '''
    Every channel has a name and some printers.
    '''
    def __init__(self, name=None):
        self.name = name

    def __repr__(self):
        return "{0}: {1}".format(self.__class__.__name__, self.name)

    def __str__(self):
        return json.dumps(self, sort_keys=True, indent=2, default=json_serializer)

class LogicalChannel(Channel):
    '''
    The main class from which we will generate sequences. 
    At some point it needs to be assigned to a physical channel.
    '''
    def __init__(self, name=None, physicalChannel=None):
        super(LogicalChannel, self).__init__(name)
        self._physicalChannel = physicalChannel

    @property
    def physicalChannel(self):
        if self._physicalChannel:
            return ChannelDict[self._physicalChannel]
        else:
            return PhysicalChannel()
    
class PhysicalChannel(Channel):
    '''
    The main class for actual AWG channels.
    '''
    def __init__(self, name=None, AWG=None, generator=None):
        super(PhysicalChannel, self).__init__(name)
        self._AWG = AWG
        self._generator = generator

    @property
    def AWG(self):
        if self._AWG:
            return ChannelDict[self._AWG]
        else:
            # create a default AWG object
            return AWG()

    @property
    def samplingRate(self):
        return self.AWG.samplingRate

    @property
    def generator(self):
        if self._generator:
            return ChannelDict[self._generator]
        else:
            return Generator()

class PhysicalMarkerChannel(PhysicalChannel):
    '''
    An digital output channel on an AWG.
    '''
    def __init__(self, name=None, AWG=None, generator=None, channel=None, delay=0.0, **kwargs):
        super(PhysicalMarkerChannel, self).__init__(name=name, AWG=AWG, generator=generator)
        self.delay = delay
        self.channel = channel
        
class PhysicalQuadratureChannel(PhysicalChannel):
    '''
    Something used to implement a standard qubit channel with two analog channels and a microwave gating channel.
    '''
    def __init__(self, name=None, AWG=None, generator=None, IChannel=None, QChannel=None, delay=0.0, ampFactor=1.0, phaseSkew=0.0, **kwargs):
        super(PhysicalQuadratureChannel, self).__init__(name=name, AWG=AWG, generator=generator)
        self.IChannel = IChannel
        self.QChannel = QChannel
        self.delay = delay
        self.ampFactor = ampFactor
        self.phaseSkew = phaseSkew

    @property
    def correctionT(self):
        return np.array([[self.ampFactor, self.ampFactor*tan(self.phaseSkew*pi/180)], [0, 1/cos(self.phaseSkew*pi/180)]])
                
class LogicalMarkerChannel(LogicalChannel):
    '''
    A class for digital channels for gating sources or triggering other things.
    '''
    def __init__(self, name=None, physicalChannel=None, **kwargs):
        super(LogicalMarkerChannel, self).__init__(name=name, physicalChannel=physicalChannel)        
    
    def gatePulse(self, length, delay=0):
        tmpBlock = PulseSequencer.PulseBlock()
        if delay>0:
            tmpBlock.add_pulse(PatternGen.QId(delay), self)
        tmpBlock.add_pulse(PatternGen.Square(length, amp=1), self)
        return tmpBlock

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

class Qubit(LogicalChannel):
    '''
    The main class for generating qubit pulses.  Effectively a logical "QuadratureChannel".
    '''
    def __init__(self, name=None, physicalChannel=None, pulseParams=None):
        super(Qubit, self).__init__(name=name, physicalChannel=physicalChannel)
        defaultPulseParams = {'length':20e-9, 'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.gaussian, 'buffer':0.0, 'cutoff':2, 'dragScaling':0}
        defaultPulseParams.update(pulseParams)
        self.pulseParams = defaultPulseParams

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

class Generator(Channel):
    '''
    Although not quite a channel, it is tightly linked to channels.
    '''
    def __init__(self, name=None, gateChannel=None, gateBuffer=0.0, gateMinWidth=0.0, gateDelay=0.0):
        super(Generator, self).__init__(name)
        self._gateChannel = gateChannel
        self.gateBuffer = gateBuffer
        self.gateMinWidth = gateMinWidth
        self.gateDelay = gateDelay

    @property
    def gateChannel(self):
        return ChannelDict[self._gateChannel]

class AWG(Channel):
    '''
    Although not quite a channel, it is tightly linked to channels.
    '''
    def __init__(self, name=None, model=None, samplingRate=1.2e9):
        super(AWG, self).__init__(name)
        self.model = model
        self.samplingRate = samplingRate

def save_channel_info(fileName=None):
    '''
    Helper function to save a channelInfo dictionary to a JSON file or string.
    '''
    #Convert the channel into a dictionary
    if fileName is None:
        fileName = config.channelParamsFile
    with open(fileName,'w') as FID:
        json.dump(ChannelDict, FID, sort_keys=True, indent=2, default=json_serializer)
        
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

def json_serializer(obj):
    '''
    Helper function to flatten the channel classes to a dictionary for json serialization.
    We just keep the class name and the properties
    '''
    if isinstance(obj, FunctionType):
        return obj.__name__
    else:
        jsonDict = {'__class__': obj.__class__.__name__}
        #Strip leading underscores off private properties
        newDict = { key.lstrip('_'):value for key,value in obj.__dict__.items()}
        jsonDict.update(newDict)
        return jsonDict

def json_deserializer(jsonDict):
    '''
    Helper function to convert a json representation of a channel back into an object.
    '''
    #Extract the class name from the dictionary
    #If there is no class then assume top level dictionary
    if '__class__' not in jsonDict:
        return jsonDict
    else:
        className = jsonDict.pop('__class__')
        class_ = getattr(sys.modules[__name__], className)
        #Deal with shape functions
        if 'pulseParams' in jsonDict:
            if 'shapeFun' in jsonDict['pulseParams']:
                jsonDict['pulseParams']['shapeFun'] = getattr(PulseShapes, jsonDict['pulseParams']['shapeFun'])
        return class_(**jsonDict)

#Load the ChannelDict on import
update_channel_info()

if __name__ == '__main__':
    # create a channel params file
    ChannelDict['q1'] = Qubit(name='q1',  physicalChannel='BBNAPS1-12', pulseParams={'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.drag, 'pulseLength':40e-9, 'bufferTime':2e-9, 'dragScaling':1})
    ChannelDict['q2'] = Qubit(name='q2', physicalChannel='BBNAPS1-34', pulseParams={'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.drag, 'pulseLength':40e-9, 'bufferTime':2e-9, 'dragScaling':1})
    ChannelDict['q1q2'] = Qubit(name='q1q2', physicalChannel='BBNAPS1-34', pulseParams={'piAmp':1.0, 'pi2Amp':0.5, 'shapeFun':PulseShapes.drag, 'pulseLength':40e-9, 'bufferTime':2e-9, 'dragScaling':1})
    ChannelDict['M-q1'] = Measurement(name='M-q1', measType='autodyne', physicalChannel='BBNAPS2-34', trigChan='digitizerTrig', pulseParams={'amp':1.0, 'shapeFun':PulseShapes.tanh, 'pulseLength':200e-9, 'bufferTime':2e-9})
    ChannelDict['M-q1q2'] = Measurement(name='M-q1q2', measType='autodyne', physicalChannel='BBNAPS2-34', trigChan='digitizerTrig', pulseParams={'amp':1.0, 'shapeFun':PulseShapes.tanh, 'pulseLength':200e-9, 'bufferTime':2e-9})
    
    ChannelDict['digitizerTrig'] = LogicalMarkerChannel(name='digitizerTrig', physicalChannel='BBNAPS1-2m1')

    ChannelDict['BBNAPS1-12'] = PhysicalQuadratureChannel(name='BBNAPS1-12', AWG='BBNAPS1', generator='QPC1-1691', IChannel='ch1', QChannel='ch2', delay=0e-9, ampFactor=1, phaseSkew=0)
    ChannelDict['BBNAPS1-34'] = PhysicalQuadratureChannel(name='BBNAPS1-34', AWG='BBNAPS1', generator='QPC1-1691', IChannel='ch3', QChannel='ch4', delay=0e-9, ampFactor=1, phaseSkew=0)
    ChannelDict['BBNAPS2-12'] = PhysicalQuadratureChannel(name='BBNAPS2-12', AWG='BBNAPS2', generator='Agilent1', IChannel='ch1', QChannel='ch2', delay=0e-9, ampFactor=1, phaseSkew=0)
    ChannelDict['BBNAPS2-34'] = PhysicalQuadratureChannel(name='BBNAPS2-34', AWG='BBNAPS2', generator='Agilent1', IChannel='ch3', QChannel='ch4', delay=0e-9, ampFactor=1, phaseSkew=0)
    ChannelDict['BBNAPS1-1m1'] = PhysicalMarkerChannel(name='BBNAPS1-1m1', AWG='BBNAPS1')
    ChannelDict['BBNAPS1-2m1'] = PhysicalMarkerChannel(name='BBNAPS1-2m1', AWG='BBNAPS1')
    ChannelDict['BBNAPS1-3m1'] = PhysicalMarkerChannel(name='BBNAPS1-3m1', AWG='BBNAPS1')

    ChannelDict['QPC1-1691'] = Generator(name='QPC1-1691', gateChannel='BBNAPS1-1m1', gateDelay=-50.0e-9, gateBuffer=20e-9, gateMinWidth=100e-9)
    ChannelDict['Agilent1'] = Generator(name='Agilent1', gateChannel='TekAWG1-1m1', gateDelay=-50.0e-9, gateBuffer=20e-9, gateMinWidth=100e-9)
    ChannelDict['Agilent2'] = Generator(name='Agilent2', gateChannel='TekAWG2-3m1', gateDelay=-50.0e-9, gateBuffer=20e-9, gateMinWidth=100e-9)   

    ChannelDict['TekAWG1'] = AWG(name='TekAWG1', model='Tek5000')
    ChannelDict['BBNAPS1'] = AWG(name='BBNAPS1', model='BBNAPS')
    ChannelDict['BBNAPS2'] = AWG(name='BBNAPS2', model='BBNAPS')

    save_channel_info()


    
    