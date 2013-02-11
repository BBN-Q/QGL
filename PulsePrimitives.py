import PulseShapes
import Channels

from math import pi, sin, cos, acos, sqrt
import numpy as np
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
    return Pulse("Ym", qubit, shape, -pi/2, 0.0)

# @_memoize
def Y90m(qubit):
    shape = qubit.pulseParams['shapeFun'](amp=qubit.pulseParams['pi2Amp'], **overrideDefaults(qubit, {}))
    return Pulse("Y90m", qubit, shape, -pi/2, 0.0)

def arb_axis_drag(qubit, nutFreq, rotAngle=0, polarAngle=0, aziAngle=0, **kwargs):
    """
    Single qubit arbitrary axis pulse implemented with phase ramping and frame change.
    For now we assume gaussian shape. 

    Parameters
    ----------
    qubit : logical channel
    rotAngle : effective rotation rotAngle (radians)
    polarAngle : polar angle of rotation axis (radians)
    aziAngle : azimuthal (radians)
    """
    params = overrideDefaults(qubit, kwargs)

    if params['length'] > 0:
        #Start from a gaussian shaped pulse
        gaussPulse = PulseShapes.gaussian(amp=1, **params)

        #To calculate the phase ramping we'll need the sampling rate
        sampRate = qubit.physicalChannel.samplingRate

        #Scale to achieve to the desired rotation
        calScale = (rotAngle/2/pi)*sampRate/sum(gaussPulse)

        #Calculate the phase ramp steps to achieve the desired Z component to the rotation axis
        phaseSteps = -2*pi*cos(polarAngle)*calScale*gaussPulse/sampRate

        #Calculate Z DRAG correction to phase steps
        #beta is a conversion between XY drag scaling and Z drag scaling
        beta = params['dragScaling']/sampRate
        instantaneousDetuning = beta*(2*pi*calScale*sin(polarAngle)*gaussPulse)**2
        phaseSteps = phaseSteps + instantaneousDetuning*(1.0/sampRate)
        #center phase ramp around the middle of the pulse time stes
        phaseRamp = np.cumsum(phaseSteps) - phaseSteps/2

        frameChange = sum(phaseSteps);

        shape = (1.0/nutFreq)*sin(polarAngle)*calScale*np.exp(1j*aziAngle)*gaussPulse*np.exp(1j*phaseRamp)

    elif abs(polarAngle) < 1e-10:
        #Otherwise assuem we have a zero-length Z rotation 
        frameChange = -rotAngle;
        shape = np.array([], dtype=np.complex128)
    else:
        raise ValueError('Non-zero transverse rotation with zero-length pulse.')
    
    return Pulse("ArbAxis", qubit, shape, 0.0, frameChange)

def AC(qubit, cliffNum):
    """

    The set of 24 Atomic Clifford single qubit pulses.

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel) 
    cliffNum : the zero-indexed Clifford number 

    Returns
    -------
    pulse object
    """

    #Figure out the approximate nutation frequency calibration from the X180 and the samplingRate
    Xp = X(qubit)
    xpulse = Xp.shape
    nutFreq = 0.5/(sum(xpulse)/qubit.physicalChannel.samplingRate);


    #Now a big else if chain for to get the specific Clifford
    if cliffNum == 0:
        #Identity gate
        return Id(qubit, length=0)
    elif cliffNum == 1:
        #X90
        return X90(qubit)
    elif cliffNum == 2:
        #X180
        return X(qubit)
    elif cliffNum == 3:
        #X90m
        return X90m(qubit)
    elif cliffNum == 4:
        #Y90
        return Y90(qubit)
    elif cliffNum == 5:
        #Y180
        return Y(qubit)
    elif cliffNum == 6:
        #Y90m
        return Y90m(qubit)
    elif cliffNum == 7:
        #Z90
        return arb_axis_drag(qubit, nutFreq, rotAngle=pi/2, length=0)
    elif cliffNum == 8:
        #Z180
        return arb_axis_drag(qubit, nutFreq, rotAngle=pi, length=0)
    elif cliffNum == 9:
        #Z90m
        return arb_axis_drag(qubit, nutFreq, rotAngle=-pi/2, length=0)
    elif cliffNum == 10:
        #X+Y 180
        return U180(qubit, phase=pi/4)
    elif cliffNum == 11:
        #X-Y 180
        return U180(qubit, phase=-pi/4)
    elif cliffNum == 12:
        #X+Z 180(Hadamard)
        return arb_axis_drag(qubit, nutFreq, rotAngle=pi, polarAngle=pi/4)
    elif cliffNum == 13:
        #X-Z 180
        return arb_axis_drag(qubit, nutFreq, rotAngle=pi, polarAngle=pi/4, aziAngle=pi)
    elif cliffNum == 14:
        #Y+Z 180
        return arb_axis_drag(qubit, nutFreq, rotAngle=pi, polarAngle=pi/4, aziAngle=pi/2)
    elif cliffNum == 15:
        #Y-Z 180
        return arb_axis_drag(qubit, nutFreq, rotAngle=pi, polarAngle=pi/4, aziAngle=-pi/2)
    elif cliffNum == 16:
        #X+Y+Z 120
        return arb_axis_drag(qubit, nutFreq, rotAngle=2*pi/3, polarAngle=acos(1/sqrt(3)), aziAngle=pi/4)
    elif cliffNum == 17:
        #X+Y+Z -120 (equivalent to -X-Y-Z 120)
        return arb_axis_drag(qubit, nutFreq, rotAngle=2*pi/3, polarAngle=pi-acos(1/sqrt(3)), aziAngle=5*pi/4)
    elif cliffNum == 18:
        #X-Y+Z 120
        return arb_axis_drag(qubit, nutFreq, rotAngle=2*pi/3, polarAngle=acos(1/sqrt(3)), aziAngle=-pi/4)
    elif cliffNum == 19:
        #X-Y+Z 120 (equivalent to -X+Y-Z 120)
        return arb_axis_drag(qubit, nutFreq, rotAngle=2*pi/3, polarAngle=pi-acos(1/sqrt(3)), aziAngle=3*pi/4)
    elif cliffNum == 20:
        #X+Y-Z 120
        return arb_axis_drag(qubit, nutFreq, rotAngle=2*pi/3, polarAngle=pi-acos(1/sqrt(3)), aziAngle=pi/4)
    elif cliffNum == 21:
        #X+Y-Z -120 (equivalent to -X-Y+Z 120)
        return arb_axis_drag(qubit, nutFreq, rotAngle=2*pi/3, polarAngle=acos(1/sqrt(3)), aziAngle=5*pi/4)
    elif cliffNum == 22:
        #-X+Y+Z 120 
        return arb_axis_drag(qubit, nutFreq, rotAngle=2*pi/3, polarAngle=acos(1/sqrt(3)), aziAngle=3*pi/4)
    elif cliffNum == 23:
        #-X+Y+Z -120 (equivalent to X-Y-Z 120)
        return arb_axis_drag(qubit, nutFreq, rotAngle=2*pi/3, polarAngle=pi-acos(1/sqrt(3)), aziAngle=-pi/4)
    else:
        raise ValueError('Clifford number must be between 0 and 23')


## two-qubit primitivies
# @_memoize
def CNOT(source, target):
    # construct (source, target) channel and pull parameters from there
    twoQChannel = Channels.QubitFactory(source.name + target.name)
    shape = twoQChannel.pulseParams['shapeFun'](amp=twoQChannel.pulseParams['piAmp'], **overrideDefaults(twoQChannel, {}))
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
