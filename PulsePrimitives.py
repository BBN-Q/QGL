import PulseShapes

from scipy.constants import pi
from PulseSequencer import Pulse
from functools import wraps

def overrideDefaults(qubit, updateParams):
    '''Helper function to update any parameters passed in and fill in the defaults otherwise.'''
    paramsList = ['shapeFun','pulseLength','bufferTime','piAmp','pi2Amp','dragScaling','cutoff']
    #First get the default or updated values
    updateValues = [updateParams[paramName] if paramName in updateParams else getattr(qubit, paramName) for paramName in paramsList]
    #Return a dictionary        
    paramDict = {paramName:paramValue for paramName,paramValue in zip(paramsList, updateValues)}
    # pull in the samplingRate from the physicalChannel
    paramDict['samplingRate'] = qubit.physicalChannel.samplingRate
    return paramDict

def _memoize(pulseFunc):
    ''' Decorator for caching pulses to keep waveform memory usage down. '''
    cache = {}
    @wraps(pulseFunc)
    def cacheWrap(*args):
        if args not in cache:
            cache[args] = pulseFunc(args)
        return cache[args]
    return cacheWrap

def Id(qubit, *args, **kwargs):
    '''
    A delay or do-nothing in the form of a pulse i.e. it will take pulseLength+2*bufferTime.
    Accepts the following pulse signatures:
        Id(qubit, [kwargs])
        Id(qubit, delay, [kwargs])
    '''
    params = overrideDefaults(qubit, kwargs)
    if len(args) > 0:
        params['pulseLength'] = args[0]
    shape = PulseShapes.delay(**params)
    return Pulse("Id", qubit, shape, 0, 0.0)

def Xtheta(qubit, amp=0, **kwargs):
    '''  A generic X rotation with a variable amplitude  '''
    params = overrideDefaults(qubit, kwargs)
    shape = params['shapeFun'](amp=amp, **params)
    return Pulse("Xtheta", qubit, shape, 0, 0.0)

def Ytheta(qubit, amp=0, **kwargs):
    ''' A generic Y rotation with a variable amplitude '''
    params = overrideDefaults(qubit, kwargs)
    shape = params['shapeFun'](amp=amp, **params)
    return Pulse("Ytheta", qubit, shape, pi/2, 0.0)
    
def U90(qubit, phase=0, **kwargs):
    ''' A generic 90 degree rotation with variable phase. '''
    params = overrideDefaults(qubit, kwargs)
    shape = params['shapeFun'](amp=qubit.pi2Amp, **params)
    return Pulse("U90", qubit, shape, phase, 0.0)

def U180(qubit, phase=0, **kwargs):
    ''' A generic 180 degree rotation with variable phase.  '''
    params = overrideDefaults(qubit, kwargs)
    shape = params['shapeFun'](amp=qubit.piAmp, **params)
    return Pulse("U180", qubit, shape, phase, 0.0)
    
def Utheta(qubit, amp=0, phase=0, **kwargs):
    '''  A generic rotation with variable amplitude and phase. '''
    params = overrideDefaults(qubit, kwargs)
    shape = params['shapeFun'](amp=amp, **params)
    return Pulse("Utheta", qubit, shape, phase, 0.0)

#Setup the default 90/180 rotations
# @_memoize
def X(qubit):
    shape = qubit.shapeFun(amp=qubit.piAmp, **overrideDefaults(qubit, {}))
    return Pulse("X", qubit, shape, 0, 0.0)
    
# @_memoize
def X90(qubit):
    shape = qubit.shapeFun(amp=qubit.pi2Amp, **overrideDefaults(qubit, {}))
    return Pulse("X90", qubit, shape, 0, 0.0)

# @_memoize
def Xm(qubit):
    shape = qubit.shapeFun(amp=qubit.piAmp, **overrideDefaults(qubit, {}))
    return Pulse("Xm", qubit, shape, pi, 0.0)
    
# @_memoize
def X90m(qubit):
    shape = qubit.shapeFun(amp=qubit.pi2Amp, **overrideDefaults(qubit, {}))
    return Pulse("X90m", qubit, shape, pi, 0.0)

# @_memoize
def Y(qubit):
    shape = qubit.shapeFun(amp=qubit.piAmp, **overrideDefaults(qubit, {}))
    return Pulse("Y", qubit, shape, pi/2, 0.0)

# @_memoize
def Y90(qubit):
    shape = qubit.shapeFun(amp=qubit.pi2Amp, **overrideDefaults(qubit, {}))
    return Pulse("Y90", qubit, shape, pi/2, 0.0)

# @_memoize
def Ym(qubit):
    shape = qubit.shapeFun(amp=qubit.piAmp, **overrideDefaults(qubit, {}))
    return Pulse("Ym", qubit, shape, -pi/4, 0.0)

# @_memoize
def Y90m(qubit):
    shape = qubit.shapeFun(amp=qubit.pi2Amp, **overrideDefaults(qubit, {}))
    return Pulse("Y90m", qubit, shape, -pi/4, 0.0)

## two-qubit primitivies
# @_memoize
def CNOT(source, target):
    # TODO construct the (source, target) channel and pull parameters from there
    # something like: channel = Qubit((source, target))
    shape = source.shapeFun(amp=source.piAmp, **overrideDefaults(source, {}))
    return Pulse("CNOT", (source, target), shape, 0.0, 0.0)