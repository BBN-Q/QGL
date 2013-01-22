import PulseShapes
import Channels

from scipy.constants import pi
from PulseSequencer import Pulse
from functools import wraps

def overrideDefaults(chan, updateParams):
    '''Helper function to update any parameters passed in and fill in the defaults otherwise.'''
    #The default parameter list depends on the channel type so pull out of channel
    #First get the default or updated values
    paramDict = chan.pulseParams.copy()
    paramDict.update(updateParams)
    # pull in the samplingRate from the physicalChannel
    paramDict['samplingRate'] = chan.physicalChannel.samplingRate
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
        params['length'] = args[0]
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
    shape = params['shapeFun'](amp=qubit.pulseParams['pi2Amp'], **params)
    return Pulse("U90", qubit, shape, phase, 0.0)

def U180(qubit, phase=0, **kwargs):
    ''' A generic 180 degree rotation with variable phase.  '''
    params = overrideDefaults(qubit, kwargs)
    shape = params['shapeFun'](amp=qubit.pulseParams['piAmp'], **params)
    return Pulse("U180", qubit, shape, phase, 0.0)
    
def Utheta(qubit, amp=0, phase=0, **kwargs):
    '''  A generic rotation with variable amplitude and phase. '''
    params = overrideDefaults(qubit, kwargs)
    shape = params['shapeFun'](amp=amp, **params)
    return Pulse("Utheta", qubit, shape, phase, 0.0)

#Setup the default 90/180 rotations
# @_memoize
def X(qubit):
    shape = qubit.pulseParams['shapeFun'](amp=qubit.pulseParams['piAmp'], **overrideDefaults(qubit, {}))
    return Pulse("X", qubit, shape, 0, 0.0)
    
# @_memoize
def X90(qubit):
    shape = qubit.pulseParams['shapeFun'](amp=qubit.pulseParams['pi2Amp'], **overrideDefaults(qubit, {}))
    return Pulse("X90", qubit, shape, 0, 0.0)

# @_memoize
def Xm(qubit):
    shape = qubit.pulseParams['shapeFun'](amp=qubit.pulseParams['piAmp'], **overrideDefaults(qubit, {}))
    return Pulse("Xm", qubit, shape, pi, 0.0)
    
# @_memoize
def X90m(qubit):
    shape = qubit.pulseParams['shapeFun'](amp=qubit.pulseParams['pi2Amp'], **overrideDefaults(qubit, {}))
    return Pulse("X90m", qubit, shape, pi, 0.0)

# @_memoize
def Y(qubit):
    shape = qubit.pulseParams['shapeFun'](amp=qubit.pulseParams['piAmp'], **overrideDefaults(qubit, {}))
    return Pulse("Y", qubit, shape, pi/2, 0.0)

# @_memoize
def Y90(qubit):
    shape = qubit.pulseParams['shapeFun'](amp=qubit.pulseParams['pi2Amp'], **overrideDefaults(qubit, {}))
    return Pulse("Y90", qubit, shape, pi/2, 0.0)

# @_memoize
def Ym(qubit):
    shape = qubit.pulseParams['shapeFun'](amp=qubit.pulseParams['piAmp'], **overrideDefaults(qubit, {}))
    return Pulse("Ym", qubit, shape, -pi/4, 0.0)

# @_memoize
def Y90m(qubit):
    shape = qubit.pulseParams['shapeFun'](amp=qubit.pulseParams['pi2Amp'], **overrideDefaults(qubit, {}))
    return Pulse("Y90m", qubit, shape, -pi/4, 0.0)

## two-qubit primitivies
# @_memoize
def CNOT(source, target):
    # construct (source, target) channel and pull parameters from there
    twoQChannel = Channels.QubitFactory(source.name + target.name)
    shape = twoQChannel.shapeFun(amp=twoQChannel.piAmp, **overrideDefaults(twoQChannel, {}))
    return Pulse("CNOT", (source, target), shape, 0.0, 0.0)

## Measurement operators
# @_memoize
def MEAS(qubit, *args, **kwargs):
    '''
    MEAS(q1, ...) constructs a measurement pulse block of a measurment + digitizer trigger.
    Use the single-argument form for an individual readout channel, e.g.
        MEAS(q1)
    Use the multi-argument form for joint readout, e.g.
        MEAS(q1, q2)
    '''
    channelName = "M-" + qubit.name
    for q in args:
        channelName += q.name
    measChannel = Channels.MeasFactory(channelName)
    params = overrideDefaults(measChannel, kwargs)
    # measurement channels should have just an "amp" parameter
    measShape = measChannel.pulseParams['shapeFun'](**params)
    return Pulse("MEAS", measChannel, measShape, 0.0, 0.0) * Pulse("trig", measChannel.trigChan, measChannel.trigChan.pulseParams['shapeFun'](**measChannel.trigChan.pulseParams), 0.0, 0.0)
