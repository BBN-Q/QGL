'''
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
import PulseShapes
import Channels
import operator

from math import pi, sin, cos, acos, sqrt
import numpy as np
from PulseSequencer import Pulse, TAPulse
from functools import wraps

def overrideDefaults(chan, updateParams):
    '''Helper function to update any parameters passed in and fill in the defaults otherwise.'''
    # The default parameter list depends on the channel type so pull out of channel
    # Then update passed values
    paramDict = chan.pulseParams.copy()
    paramDict.update(updateParams)
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

def Id(channel, *args, **kwargs):
    '''
    A delay or no-op in the form of a pulse.
    Accepts the following pulse signatures:
        Id(channel, [kwargs])
        Id(channel, delay, [kwargs])
    '''
    params = overrideDefaults(channel, kwargs)
    if len(args) > 0 and isinstance(args[0], (int,float)):
        params['length'] = args[0]

    return TAPulse("Id", channel, params['length'], 0)

# the most generic pulse is Utheta
def Utheta(qubit, amp=0, phase=0, label='Utheta', **kwargs):
    '''  A generic rotation with variable amplitude and phase. '''
    params = overrideDefaults(qubit, kwargs)
    return Pulse(label, qubit, params, amp, phase, 0.0)

# generic pulses around X, Y, and Z axes
def Xtheta(qubit, amp=0, label='Xtheta', **kwargs):
    '''  A generic X rotation with a variable amplitude  '''
    return Utheta(qubit, amp, 0, label=label, **kwargs)

def Ytheta(qubit, amp=0, label='Ytheta', **kwargs):
    ''' A generic Y rotation with a variable amplitude '''
    return Utheta(qubit, amp, pi/2, label=label, **kwargs)

def Ztheta(qubit, angle=0, label='Ztheta', **kwargs):
    # special cased because it can be done with a frame update
    return TAPulse(label, qubit, length=0, amp=0, phase=0, frameChange=-angle)

#Setup the default 90/180 rotations
# @_memoize
def X(qubit, **kwargs):
    return Xtheta(qubit, qubit.pulseParams['piAmp'], label="X", **kwargs)

# @_memoize
def X90(qubit, **kwargs):
    return Xtheta(qubit, qubit.pulseParams['pi2Amp'], label="X90", **kwargs)

# @_memoize
def Xm(qubit, **kwargs):
    return Xtheta(qubit, -qubit.pulseParams['piAmp'], label="Xm", **kwargs)

# @_memoize
def X90m(qubit, **kwargs):
    return Xtheta(qubit, -qubit.pulseParams['pi2Amp'], label="X90m", **kwargs)

# @_memoize
def Y(qubit, **kwargs):
    return Ytheta(qubit, qubit.pulseParams['piAmp'], label="Y", **kwargs)

# @_memoize
def Y90(qubit, **kwargs):
    return Ytheta(qubit, qubit.pulseParams['pi2Amp'], label="Y90", **kwargs)

# @_memoize
def Ym(qubit, **kwargs):
    return Ytheta(qubit, -qubit.pulseParams['piAmp'], label="Ym", **kwargs)

# @_memoize
def Y90m(qubit, **kwargs):
    return Ytheta(qubit, -qubit.pulseParams['pi2Amp'], label="Y90m", **kwargs)

# @_memoize
def Z(qubit, **kwargs):
    return Ztheta(qubit, pi, label="Z", **kwargs)

# @_memoize
def Z90(qubit, **kwargs):
    return Ztheta(qubit, pi/2, label="Z90", **kwargs)

# @_memoize
def Z90m(qubit, **kwargs):
    return Ztheta(qubit, -pi/2, label="Z90m", **kwargs)

# 90/180 degree rotations with control over the rotation axis
def U90(qubit, phase=0, **kwargs):
    ''' A generic 90 degree rotation with variable phase. '''
    return Utheta(qubit, qubit.pulseParams['pi2Amp'], phase, label="U90", **kwargs)

def U(qubit, phase=0, **kwargs):
    ''' A generic 180 degree rotation with variable phase.  '''
    return Utheta(qubit, qubit.pulseParams['piAmp'], phase, label="U", **kwargs)

def arb_axis_drag(qubit, nutFreq, rotAngle=0, polarAngle=0, aziAngle=0, **kwargs):
    """
    Single qubit arbitrary axis pulse implemented with phase ramping and frame change.
    For now we assume gaussian shape.

    Parameters
    ----------
    qubit : logical channel
    nutFreq: effective nutation frequency per unit of drive amplitude (Hz)
    rotAngle : effective rotation rotAngle (radians)
    polarAngle : polar angle of rotation axis (radians)
    aziAngle : azimuthal (radians)
    """
    params = overrideDefaults(qubit, kwargs)

    # TODO: figure out way to reduce code duplication between this and the pulse shape
    if params['length'] > 0:
        #To calculate the phase ramping we'll need the sampling rate
        sampRate = qubit.physChan.samplingRate

        #Start from a gaussian shaped pulse
        gaussPulse = PulseShapes.gaussian(amp=1, samplingRate=sampRate, **params).real

        #Scale to achieve to the desired rotation
        calScale = (rotAngle/2/pi)*sampRate/sum(gaussPulse)

        #Calculate the phase ramp steps to achieve the desired Z component to the rotation axis
        phaseSteps = -2*pi*cos(polarAngle)*calScale*gaussPulse/sampRate

        #Calculate Z DRAG correction to phase steps
        #beta is a conversion between XY drag scaling and Z drag scaling
        beta = params['dragScaling']/sampRate
        instantaneousDetuning = beta*(2*pi*calScale*sin(polarAngle)*gaussPulse)**2
        phaseSteps = phaseSteps + instantaneousDetuning*(1.0/sampRate)

        frameChange = sum(phaseSteps);

    elif abs(polarAngle) < 1e-10:
        #Otherwise assume we have a zero-length Z rotation
        frameChange = -rotAngle;
    else:
        raise ValueError('Non-zero transverse rotation with zero-length pulse.')

    params['nutFreq'] = nutFreq
    params['rotAngle'] = rotAngle
    params['polarAngle'] = polarAngle
    params['shapeFun'] = PulseShapes.arb_axis_drag
    return Pulse("ArbAxis", qubit, params, 1.0, aziAngle, frameChange)

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
    xpulse = Xp.amp * Xp.shape
    nutFreq = 0.5/(sum(xpulse)/qubit.physChan.samplingRate);


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
        return Z90(qubit)
    elif cliffNum == 8:
        #Z180
        return Z(qubit)
    elif cliffNum == 9:
        #Z90m
        return Z90m(qubit)
    elif cliffNum == 10:
        #X+Y 180
        return U(qubit, phase=pi/4)
    elif cliffNum == 11:
        #X-Y 180
        return U(qubit, phase=-pi/4)
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
def CNOT(source, target, **kwargs):
    # construct (source, target) channel and pull parameters from there
    channel = Channels.QubitFactory(source.label + target.label)
    params = overrideDefaults(channel, kwargs)
    return Pulse("CNOT", (source, target), params, channel.pulseParams['piAmp'], 0.0, 0.0)

def flat_top_gaussian(chan, riseFall, length, amp, phase=0):
    """
    A square pulse with rising and falling gaussian shape
    """
    return Utheta(chan, length=riseFall, amp=amp, phase=phase, shapeFun=PulseShapes.gaussOn) + \
           Utheta(chan, length=length, amp=amp, phase=phase, shapeFun=PulseShapes.square) + \
           Utheta(chan, length=riseFall, amp=amp, phase=phase, shapeFun=PulseShapes.gaussOff)

def echoCR(controlQ, targetQ, amp=1, phase=0, length=200e-9, riseFall=20e-9, lastPi=True):
    """
    An echoed CR pulse.  Used for calibration of CR gate
    """
    CRchan = Channels.EdgeFactory(controlQ, targetQ)
    if not CRchan.isforward(controlQ, targetQ):
        raise ValueError('Could not find an edge with control qubit {0}'.format(controlQ))
    seq = [flat_top_gaussian(CRchan, amp=amp, riseFall=riseFall, length=length, phase=phase),
           X(controlQ),
           flat_top_gaussian(CRchan, amp=amp, riseFall=riseFall, length=length, phase=phase+np.pi)]
    if lastPi:
        seq += [X(controlQ)]
    return seq

def ZX90_CR(controlQ, targetQ, **kwargs):
    """
    A calibrated CR ZX90 pulse.  Uses 'amp' for the pulse amplitude, 'phase' for its phase (in deg).
    """
    CRchan = Channels.EdgeFactory(controlQ, targetQ)
    params = overrideDefaults(CRchan, kwargs)
    return echoCR(controlQ, targetQ, amp = params['amp'], phase = params['phase']/180*np.pi, length=params['length'], riseFall=params['riseFall'])

def CNOT_CR(controlQ, targetQ, **kwargs):
    edge = Channels.EdgeFactory(controlQ, targetQ)

    if edge.isforward(controlQ, targetQ):
        # control and target for CNOT and CR match
        return ZX90_CR(controlQ, targetQ, **kwargs) + [Z90m(controlQ) * X90m(targetQ)]
    else:
        # control and target for CNOT and CR are inverted
        return [Y90(controlQ) * Y90(targetQ),
                X(controlQ) * X(targetQ)] + \
                ZX90_CR(targetQ, controlQ, **kwargs) + \
               [Z90(targetQ),
                X90(controlQ) * Y90(targetQ),
                Y90m(controlQ) * X(targetQ)]

## Measurement operators
# @_memoize
def MEAS(*args, **kwargs):
    '''
    MEAS(q1, ...) constructs a measurement pulse block of a measurment
    Use the single-argument form for an individual readout channel, e.g.
        MEAS(q1)
    Use tuple-form for joint readout, e.g.
        MEAS((q1, q2))
    Use multi-argument form for joint simultaneous readout.
    '''
    def create_meas_pulse(qubit):
        if isinstance(qubit, Channels.Qubit):
            #Deal with single qubit readout channel
            channelName = "M-" + qubit.label
        elif isinstance(qubit, tuple):
            #Deal with joint readout channel
            channelName = "M-"
            for q in qubit:
                channelName += q.label
        measChan = Channels.MeasFactory(channelName)
        params = overrideDefaults(measChan, kwargs)
        params['frequency'] = measChan.autodyneFreq
        params['baseShape'] = params.pop('shapeFun')
        params['shapeFun'] = PulseShapes.autodyne
        amp = params.pop('amp')
        return Pulse("MEAS", measChan, params, amp, 0.0, 0.0)

    return reduce(operator.mul, [create_meas_pulse(qubit) for qubit in args])

# Gating/blanking pulse primitives
def BLANK(chan, length):
    return TAPulse("BLANK", chan.gateChan, length, 1, 0, 0)
