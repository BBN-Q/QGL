'''
Channels is where we store information for mapping virtual (qubit) channel to real channels.

Created on Jan 19, 2012

@author: cryan


'''

import sys
import json
import PulseShapes
import config

from operator import itemgetter
from math import tan,cos,pi
import config

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
        return [[self.ampFactor, self.ampFactor*tan(self.phaseSkew*pi/180)], [0, 1/cos(self.phaseSkew*pi/180)]]
                
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
    if name in ChannelDict and isinstance(ChannelDict[name], Qubit):
        return ChannelDict[name]
    else:
        return Qubit(name, **kwargs)


class Qubit(LogicalChannel):
    '''
    The main class for generating qubit pulses.  Effectively a logical "QuadratureChannel".
    '''
    def __init__(self, name=None, physicalChannel=None, piAmp=1.0, pi2Amp=0.5, shapeFun=PulseShapes.gaussian, pulseLength=20.0e-9, bufferTime=0.0, dragScaling=0, cutoff=2, **kwargs):
        super(Qubit, self).__init__(name=name, physicalChannel=physicalChannel)
        self.shapeFun = shapeFun
        self.pulseLength = pulseLength
        self.bufferTime = bufferTime
        self.piAmp = piAmp
        self.pi2Amp = pi2Amp
        self.dragScaling = dragScaling
        self.cutoff = cutoff

class Generator(Channel):
    '''
    Although not quite a channel, it is tightly linked to channels.
    '''
    def __init__(self, name=None, gateChannel=None, gateBuffer=0.0, gateMinWidth=0.0, gateDelay=0.0):
        super(Generator, self).__init__(name)
        self.gateChannel = gateChannel
        self.gateBuffer = gateBuffer
        self.gateMinWidth = gateMinWidth
        self.gateDelay = gateDelay

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
        ChannelDict.update(json.load(FID, object_hook=json_deserializer))

def json_serializer(obj):
    '''
    Helper function to flatten the channel classes to a dictionary for json serialization.
    We just keep the class name and the properties
    '''
    jsonDict = {'__class__': obj.__class__.__name__}
    #Strip leading underscores off private properties
    newDict = { key.lstrip('_'):value for key,value in obj.__dict__.items()}
    jsonDict.update(newDict)
    #Deal with shape function handles specially
    if 'shapeFun' in jsonDict:
        jsonDict['shapeFun'] = jsonDict['shapeFun'].__name__
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
        if 'shapeFun' in jsonDict:
            jsonDict['shapeFun'] = getattr(PulseShapes, jsonDict['shapeFun'])
        return class_(**jsonDict)

if __name__ == '__main__':
    # create a channel params file
    ChannelDict['q1'] = Qubit(name='q1', piAmp=1.0, pi2Amp=0.5, shapeFun=PulseShapes.drag, pulseLength=40e-9, bufferTime=2e-9, dragScaling=1, physicalChannel='BBNAPS1-12')
    ChannelDict['q2'] = Qubit(name='q2', piAmp=1.0, pi2Amp=0.5, shapeFun=PulseShapes.gaussian, pulseLength=40e-9, bufferTime=2e-9, dragScaling=1, physicalChannel='BBNAPS1-34')
    ChannelDict['M-q1'] = Qubit(name='M-q1', piAmp=1.0, pi2Amp=1.0, shapeFun=PulseShapes.square, pulseLength=200e-9, bufferTime=2e-9, dragScaling=0, physicalChannel='BBNAPS2-34')

    ChannelDict['digitizerTrig'] = LogicalMarkerChannel(name='digitizerTrig', physicalChannel='BBNAPS1-2m1')

    ChannelDict['BBNAPS1-12'] = PhysicalQuadratureChannel(name='BBNAPS1-12', AWG='BBNAPS1', generator='QPC1-1691', IChannel='ch1', QChannel='ch2', delay=0e-9, ampFactor=1, phaseSkew=0)
    ChannelDict['BBNAPS1-34'] = PhysicalQuadratureChannel(name='BBNAPS1-34', AWG='BBNAPS1', generator='QPC1-1691', IChannel='ch3', QChannel='ch4', delay=0e-9, ampFactor=1, phaseSkew=0)
    ChannelDict['BBNAPS2-12'] = PhysicalQuadratureChannel(name='BBNAPS2-12', AWG='BBNAPS2', generator='Agilent1', IChannel='ch1', QChannel='ch2', delay=0e-9, ampFactor=1, phaseSkew=0)
    ChannelDict['BBNAPS2-34'] = PhysicalQuadratureChannel(name='BBNAPS2-34', AWG='BBNAPS2', generator='Agilent1', IChannel='ch3', QChannel='ch4', delay=0e-9, ampFactor=1, phaseSkew=0)
    ChannelDict['BBNAPS1-2m1'] = PhysicalMarkerChannel(name='BBNAPS1-2m1', AWG='BBNAPS1')

    ChannelDict['QPC1-1691'] = Generator(name='QPC1-1691', gateChannel='TekAWG1-ch1m1', gateDelay=-50.0e-9, gateBuffer=20e-9, gateMinWidth=100e-9)
    ChannelDict['Agilent1'] = Generator(name='Agilent1', gateChannel='TekAWG1-ch1m1', gateDelay=-50.0e-9, gateBuffer=20e-9, gateMinWidth=100e-9)
    ChannelDict['Agilent2'] = Generator(name='Agilent2', gateChannel='TekAWG2-ch3m1', gateDelay=-50.0e-9, gateBuffer=20e-9, gateMinWidth=100e-9)   

    ChannelDict['TekAWG1'] = AWG(name='TekAWG1', model='Tek5000')
    ChannelDict['BBNAPS1'] = AWG(name='BBNAPS1', model='BBNAPS')
    ChannelDict['BBNAPS2'] = AWG(name='BBNAPS2', model='BBNAPS')

    save_channel_info()


    
    