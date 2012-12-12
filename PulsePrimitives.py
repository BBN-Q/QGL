# import PulseSequencer
# import Channels
import PulseShapes

from scipy.constants import pi
from PulseSequencer import Pulse

def overrideDefaults(qubit, updateParams):
    '''Helper function to update any parameters passed in and fill in the defaults otherwise.'''
    paramsList = ['pulseType','pulseLength','bufferTime','piAmp','pi2Amp','dragScaling','cutoff']
    #First get the default or updated values
    updateValues = [updateParams[paramName] if paramName in updateParams else getattr(qubit, paramName) for paramName in paramsList]
    #Return a dictionary        
    return {paramName:paramValue for paramName,paramValue in zip(paramsList, updateValues)}

# TODO: update for global caching (no longer chached in channels)
def _cachedPulse(pulseFunc):
    ''' Decorator for caching pulses to keep waveform memory usage down. '''
    def cacheWrap(self):
        if pulseFunc.__name__ not in self.pulseCache:
            self.pulseCache[pulseFunc.__name__] = pulseFunc(self)
        return self.pulseCache[pulseFunc.__name__]
    
    return cacheWrap

def Id(qubit, **kwargs):
    ''' A delay or do-nothing in the form of a pulse i.e. it will take pulseLength+2*bufferTime. '''
    shape = PulseShapes.Delay(**overrideDefaults(qubit, kwargs))
    return Pulse("Id", (qubit), shape, 0.0)

def Xtheta(qubit, amp=0, **kwargs):
    '''  A generic X rotation with a variable amplitude  '''
    params = overrideDefaults(qubit, kwargs)
    shape = getattr(PulseShapes, kwargs['shapeFun'])(amp=amp, phase=0, **params)
    return Pulse("Xtheta", (qubit), shape, 0.0)

def Ytheta(qubit, amp=0, **kwargs):
    ''' A generic Y rotation with a variable amplitude '''
    params = overrideDefaults(qubit, kwargs)
    shape = getattr(PulseShapes, kwargs['shapeFun'])(amp=amp, phase=pi/2, **params)
    return Pulse("Ytheta", (qubit), shape, 0.0)
    
def U90(qubit, phase=0, **kwargs):
    ''' A generic 90 degree rotation with variable phase. '''
    params = overrideDefaults(qubit, kwargs)
    shape = getattr(PulseShapes, kwargs['shapeFun'])(amp=qubit.pi2Amp, phase=phase, **params)
    return Pulse("U90", (qubit), shape, 0.0)

def U180(qubit, phase=0, **kwargs):
    ''' A generic 180 degree rotation with variable phase.  '''
    params = overrideDefaults(qubit, kwargs)
    shape = getattr(PulseShapes, kwargs['shapeFun'])(amp=qubit.piAmp, phase=phase, **params)
    return Pulse("U180", (qubit), shape, 0.0)
    
def Utheta(qubit, amp=0, phase=0, **kwargs):
    '''  A generic rotation with variable amplitude and phase. '''
    params = overrideDefaults(qubit, kwargs)
    shape = getattr(PulseShapes, kwargs['shapeFun'])(amp=amp, phase=phase, **params)
    return Pulse("Utheta", (qubit), shape, 0.0)

#Setup the default 90/180 rotations
# @_cachedPulse
def X(qubit):
    shape = getattr(PulseShapes, qubit.shapeFun)(amp=qubit.piAmp, phase=0, **overrideDefaults(qubit, {}))
    return Pulse("X", (qubit), shape, 0.0)
    
# @_cachedPulse
def X90(qubit):
    shape = getattr(PulseShapes, qubit.shapeFun)(amp=qubit.pi2Amp, phase=0, **overrideDefaults(qubit, {}))
    return Pulse("X90", (qubit), shape, 0.0)

# @_cachedPulse
def Xm(qubit):
    shape = getattr(PulseShapes, qubit.shapeFun)(amp=qubit.piAmp, phase=0.5, **overrideDefaults(qubit, {}))
    return Pulse("Xm", (qubit), shape, 0.0)
    
# @_cachedPulse
def X90m(qubit):
    shape = getattr(PulseShapes, qubit.shapeFun)(amp=qubit.pi2Amp, phase=0.5, **overrideDefaults(qubit, {}))
    return Pulse("X90m", (qubit), shape, 0.0)

# @_cachedPulse
def Y(qubit):
    shape = getattr(PulseShapes, qubit.shapeFun)(amp=qubit.piAmp, phase=0.25, **overrideDefaults(qubit, {}))
    return Pulse("Y", (qubit), shape, 0.0)

# @_cachedPulse
def Y90(qubit):
    shape = getattr(PulseShapes, qubit.shapeFun)(amp=qubit.pi2Amp, phase=0.25, **overrideDefaults(qubit, {}))
    return Pulse("Y90", (qubit), shape, 0.0)

# @_cachedPulse
def Ym(qubit):
    shape = getattr(PulseShapes, qubit.shapeFun)(amp=qubit.piAmp, phase=0.75, **overrideDefaults(qubit, {}))
    return Pulse("Ym", (qubit), shape, 0.0)

# @_cachedPulse
def Y90m(qubit):
    shape = getattr(PulseShapes, qubit.shapeFun)(amp=qubit.pi2Amp, phase=0.75, **overrideDefaults(qubit, {}))
    return Pulse("Y90m", (qubit), shape, 0.0)

## two-qubit primitivies
# @_cachedPulse
def CNOT(source, target):
    # TODO construct the (source, target) channel and pull parameters from there
    # something like: channel = Qubit((source, target))
    shape = getattr(PulseShapes, source.shapeFun)(amp=source.piAmp, phase=0.0, **overrideDefaults(source, {}))
    return Pulse("CNOT", (source, target), shape, 0.0)