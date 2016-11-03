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
from . import PulseShapes
from . import Channels
from . import ChannelLibrary
from . import config
import operator

from math import pi, sin, cos, acos, sqrt
import numpy as np
from .PulseSequencer import Pulse, TAPulse, align
from functools import wraps, reduce


def overrideDefaults(chan, updateParams):
    '''Helper function to update any parameters passed in and fill in the defaults otherwise.'''
    # The default parameter list depends on the channel type so pull out of channel
    # Then update passed values
    paramDict = chan.pulseParams.copy()
    paramDict.update(updateParams)
    return paramDict


def _memoize(pulseFunc):
    """ Decorator for caching pulses to keep waveform memory usage down. """
    # namespacce the cache so we can access (and reset) from elsewhere
    _memoize.cache = {}

    @wraps(pulseFunc)
    def cacheWrap(*args, **kwargs):
        if kwargs:
            return pulseFunc(*args, **kwargs)
        key = (pulseFunc, args)
        if key not in _memoize.cache:
            _memoize.cache[key] = pulseFunc(*args)
        return _memoize.cache[key]

    return cacheWrap

def clear_pulse_cache():
    _memoize.cache = {}

@_memoize
def Id(channel, *args, **kwargs):
    '''
    A delay or no-op in the form of a pulse.
    Accepts the following pulse signatures:
        Id(channel, [kwargs])
        Id(channel, delay, [kwargs])
    '''
    params = overrideDefaults(channel, kwargs)
    if len(args) > 0 and isinstance(args[0], (int, float)):
        params['length'] = args[0]

    return TAPulse("Id",
                   channel,
                   params['length'],
                   0,
                   ignoredStrParams=['amp', 'phase', 'frameChange'])


# the most generic pulse is Utheta
def Utheta(qubit,
           amp=0,
           phase=0,
           label='Utheta',
           ignoredStrParams=[],
           **kwargs):
    '''  A generic rotation with variable amplitude and phase. '''
    params = overrideDefaults(qubit, kwargs)
        #amp and phase are now pulse parameters rather than shape parameters
    if "amp" in params:
      del params["amp"]
    if "phase" in params:
      del params["phase"]
    return Pulse(label, qubit, params, amp, phase, 0.0, ignoredStrParams)


# generic pulses around X, Y, and Z axes
def Xtheta(qubit, amp=0, label='Xtheta', ignoredStrParams=[], **kwargs):
    '''  A generic X rotation with a variable amplitude  '''
    ignoredStrParams += ['phase', 'frameChange']
    return Utheta(qubit, amp, 0, label, ignoredStrParams, **kwargs)


def Ytheta(qubit, amp=0, label='Ytheta', ignoredStrParams=[], **kwargs):
    ''' A generic Y rotation with a variable amplitude '''
    ignoredStrParams += ['phase', 'frameChange']
    return Utheta(qubit, amp, pi / 2, label, ignoredStrParams, **kwargs)


def Ztheta(qubit,
           angle=0,
           label='Ztheta',
           ignoredStrParams=['amp', 'phase', 'length'],
           **kwargs):
    # special cased because it can be done with a frame update
    return TAPulse(label,
                   qubit,
                   length=0,
                   amp=0,
                   phase=0,
                   frameChange=-angle,
                   ignoredStrParams=ignoredStrParams)


#Setup the default 90/180 rotations
@_memoize
def X90(qubit, **kwargs):
    return Xtheta(qubit,
                  qubit.pulseParams['pi2Amp'],
                  label="X90",
                  ignoredStrParams=['amp'],
                  **kwargs)

@_memoize
def X90m(qubit, **kwargs):
    return Xtheta(qubit,
                  -qubit.pulseParams['pi2Amp'],
                  label="X90m",
                  ignoredStrParams=['amp'],
                  **kwargs)

@_memoize
def Y90(qubit, **kwargs):
    return Ytheta(qubit,
                  qubit.pulseParams['pi2Amp'],
                  label="Y90",
                  ignoredStrParams=['amp'],
                  **kwargs)

@_memoize
def Y90m(qubit, **kwargs):
    return Ytheta(qubit,
                  -qubit.pulseParams['pi2Amp'],
                  label="Y90m",
                  ignoredStrParams=['amp'],
                  **kwargs)

#90 degree rotation with control over the rotation axis
@_memoize
def U90(qubit, phase=0, **kwargs):
    """ A generic 90 degree rotation with variable phase. """
    if "label" not in kwargs:
        kwargs["label"] = "U90"
    return Utheta(qubit,
        qubit.pulseParams['pi2Amp'],
        phase,
        ignoredStrParams=['amp'],
        **kwargs)

if config.pulse_primitives_lib == 'standard':
    # pi rotations formed by different choice of pulse amplitude
    @_memoize
    def X(qubit, **kwargs):
        return Xtheta(qubit,
                      qubit.pulseParams['piAmp'],
                      label="X",
                      ignoredStrParams=['amp'],
                      **kwargs)

    @_memoize
    def Xm(qubit, **kwargs):
        return Xtheta(qubit,
                      -qubit.pulseParams['piAmp'],
                      label="Xm",
                      ignoredStrParams=['amp'],
                      **kwargs)

    @_memoize
    def Y(qubit, **kwargs):
        return Ytheta(qubit,
                      qubit.pulseParams['piAmp'],
                      label="Y",
                      ignoredStrParams=['amp'],
                      **kwargs)

    @_memoize
    def Ym(qubit, **kwargs):
        return Ytheta(qubit,
                      -qubit.pulseParams['piAmp'],
                      label="Ym",
                      ignoredStrParams=['amp'],
                      **kwargs)

    @_memoize
    def U(qubit, phase=0, **kwargs):
        """ A generic 180 degree rotation with variable phase.  """
        if "label" not in kwargs:
            kwargs["label"] = "U"
        return Utheta(qubit,
            qubit.pulseParams['piAmp'],
            phase,
            ignoredStrParams=['amp'],
            **kwargs)

elif config.pulse_primitives_lib == 'all90':
    # pi rotations formed by two pi/2 rotations
    @_memoize
    def X(qubit, **kwargs):
        return X90(qubit, **kwargs) + X90(qubit, **kwargs)

    @_memoize
    def Xm(qubit, **kwargs):
        return X90m(qubit, **kwargs) + X90m(qubit, **kwargs)

    @_memoize
    def Y(qubit, **kwargs):
        return Y90(qubit, **kwargs) + Y90(qubit, **kwargs)

    @_memoize
    def Ym(qubit, **kwargs):
        return Y90m(qubit, **kwargs) + Y90m(qubit, **kwargs)

    @_memoize
    def U(qubit, phase=0, **kwargs):
        """ A generic 180 degree rotation with variable phase.  """
        return U90(qubit, phase, *kwargs) + U90(qubit, phase, *kwargs)

else:
    raise Exception("Invalid pulse library")

@_memoize
def Z(qubit, **kwargs):
    return Ztheta(qubit, pi, label="Z", **kwargs)


@_memoize
def Z90(qubit, **kwargs):
    return Ztheta(qubit, pi / 2, label="Z90", **kwargs)


@_memoize
def Z90m(qubit, **kwargs):
    return Ztheta(qubit, -pi / 2, label="Z90m", **kwargs)


def arb_axis_drag(qubit,
                  nutFreq,
                  rotAngle=0,
                  polarAngle=0,
                  aziAngle=0,
                  **kwargs):
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
        gaussPulse = PulseShapes.gaussian(amp=1,
                                          samplingRate=sampRate,
                                          **params).real

        #Scale to achieve to the desired rotation
        calScale = (rotAngle / 2 / pi) * sampRate / sum(gaussPulse)

        #Calculate the phase ramp steps to achieve the desired Z component to the rotation axis
        phaseSteps = -2 * pi * cos(
            polarAngle) * calScale * gaussPulse / sampRate

        #Calculate Z DRAG correction to phase steps
        #beta is a conversion between XY drag scaling and Z drag scaling
        beta = params['dragScaling'] / sampRate
        instantaneousDetuning = beta * (2 * pi * calScale * sin(polarAngle) *
                                        gaussPulse)**2
        phaseSteps = phaseSteps + instantaneousDetuning * (1.0 / sampRate)

        frameChange = sum(phaseSteps)

    elif abs(polarAngle) < 1e-10:
        #Otherwise assume we have a zero-length Z rotation
        frameChange = -rotAngle
    else:
        raise ValueError(
            'Non-zero transverse rotation with zero-length pulse.')

    params['nutFreq'] = nutFreq
    params['rotAngle'] = rotAngle
    params['polarAngle'] = polarAngle
    params['shapeFun'] = PulseShapes.arb_axis_drag
    return Pulse(kwargs["label"] if "label" in kwargs else "ArbAxis", qubit,
                 params, 1.0, aziAngle, frameChange)


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
    nutFreq = 0.5 / (sum(xpulse) / qubit.physChan.samplingRate)

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
        return U(qubit, phase=pi / 4, label="AC_10")
    elif cliffNum == 11:
        #X-Y 180
        return U(qubit, phase=-pi / 4, label="AC_11")
    elif cliffNum == 12:
        #X+Z 180(Hadamard)
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=pi,
                             polarAngle=pi / 4,
                             label="AC_12")
    elif cliffNum == 13:
        #X-Z 180
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=pi,
                             polarAngle=pi / 4,
                             aziAngle=pi,
                             label="AC_13")
    elif cliffNum == 14:
        #Y+Z 180
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=pi,
                             polarAngle=pi / 4,
                             aziAngle=pi / 2,
                             label="AC_14")
    elif cliffNum == 15:
        #Y-Z 180
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=pi,
                             polarAngle=pi / 4,
                             aziAngle=-pi / 2,
                             label="AC_15")
    elif cliffNum == 16:
        #X+Y+Z 120
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=2 * pi / 3,
                             polarAngle=acos(1 / sqrt(3)),
                             aziAngle=pi / 4,
                             label="AC_16")
    elif cliffNum == 17:
        #X+Y+Z -120 (equivalent to -X-Y-Z 120)
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=2 * pi / 3,
                             polarAngle=pi - acos(1 / sqrt(3)),
                             aziAngle=5 * pi / 4,
                             label="AC_17")
    elif cliffNum == 18:
        #X-Y+Z 120
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=2 * pi / 3,
                             polarAngle=acos(1 / sqrt(3)),
                             aziAngle=-pi / 4,
                             label="AC_18")
    elif cliffNum == 19:
        #X-Y+Z 120 (equivalent to -X+Y-Z -120)
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=2 * pi / 3,
                             polarAngle=pi - acos(1 / sqrt(3)),
                             aziAngle=3 * pi / 4,
                             label="AC_19")
    elif cliffNum == 20:
        #X+Y-Z 120
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=2 * pi / 3,
                             polarAngle=pi - acos(1 / sqrt(3)),
                             aziAngle=pi / 4,
                             label="AC_20")
    elif cliffNum == 21:
        #X+Y-Z -120 (equivalent to -X-Y+Z 120)
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=2 * pi / 3,
                             polarAngle=acos(1 / sqrt(3)),
                             aziAngle=5 * pi / 4,
                             label="AC_21")
    elif cliffNum == 22:
        #-X+Y+Z 120
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=2 * pi / 3,
                             polarAngle=acos(1 / sqrt(3)),
                             aziAngle=3 * pi / 4,
                             label="AC_22")
    elif cliffNum == 23:
        #-X+Y+Z -120 (equivalent to X-Y-Z 120)
        return arb_axis_drag(qubit,
                             nutFreq,
                             rotAngle=2 * pi / 3,
                             polarAngle=pi - acos(1 / sqrt(3)),
                             aziAngle=-pi / 4,
                             label="AC_23")
    else:
        raise ValueError('Clifford number must be between 0 and 23')


def DiAC(qubit, cliffNum, Xonly = False):
    """

    The set of 24 Diatomic Clifford single qubit pulses. Use X90(m), Z, Z90
    if Xonly; use also Y90(m) otherwise

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    cliffNum : the zero-indexed Clifford number

    Returns
    -------
    pulse object
    """
    #Now a big else if chain for to get the specific Clifford
    if cliffNum == 0:
        #Identity gate
        return Id(qubit, length=0)
    elif cliffNum == 1:
        #X90
        return X90(qubit)
    elif cliffNum == 2:
        #X180
        return X90(qubit)+X90(qubit)
    elif cliffNum == 3:
        #X90m
        return X90m(qubit)
    elif cliffNum == 4:
        #Y90
        return X90(qubit) + Z90m(qubit) + X90(qubit) + Z(qubit) if Xonly else Y90(qubit)
    elif cliffNum == 5:
        #Y180
        return X90(qubit) + X90(qubit) + Z(qubit) if Xonly else Y90(qubit) + Y90(qubit)
    elif cliffNum == 6:
        #Y90m
        return X90(qubit) + Z90(qubit) + X90(qubit) + Z(qubit) if Xonly else Y90m(qubit)
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
        return X90(qubit) + X90(qubit) + Z90(qubit) if Xonly else Y90(qubit) + Y90(qubit) + Z90m(qubit)
    elif cliffNum == 11:
        #X-Y 180
        return X90(qubit) + X90(qubit) + Z90m(qubit) if Xonly else Y90(qubit) + Y90(qubit) + Z90(qubit)
    elif cliffNum == 12:
        #X+Z 180(Hadamard)
        return Z(qubit) + X90(qubit) + Z90m(qubit) + X90(qubit) + Z(qubit) if Xonly else Z(qubit) + Y90(qubit)
    elif cliffNum == 13:
        #X-Z 180
        return Z(qubit) + X90(qubit) + Z90(qubit) + X90(qubit) + Z(qubit) if Xonly else Z(qubit) + Y90m(qubit)
    elif cliffNum == 14:
        #Y+Z 180
        return Z90(qubit) + X90(qubit) + Z90m(qubit) + X90(qubit) + Z90m(qubit) if Xonly else Z90(qubit) + Y90(qubit) + Z90(qubit)
    elif cliffNum == 15:
        #Y-Z 180
        return Z90(qubit) + X90(qubit) + Z90(qubit) + X90(qubit) + Z90m(qubit) if Xonly else Z90(qubit) + Y90m(qubit) + Z90(qubit)
    elif cliffNum == 16:
            #X+Y+Z -120 (equivalent to -X-Y-Z 120)
        return Z90(qubit) + X90(qubit) + Z90m(qubit) + X90(qubit) + Z(qubit) if Xonly else Z90(qubit) + Y90(qubit)
    elif cliffNum == 17:
        #X+Y+Z 120
        return Z(qubit) + X90(qubit) + Z90m(qubit) + X90(qubit) + Z90m(qubit) if Xonly else Z(qubit) + Y90(qubit) + Z90(qubit)
    elif cliffNum == 18:
            #X-Y+Z 120 (equivalent to -X+Y-Z 120)
            return X90(qubit) + Z90(qubit) + X90(qubit) + Z90m(qubit) if Xonly else Y90m(qubit) + Z90(qubit)
    elif cliffNum == 19:
        #X-Y+Z -120
        return Z90m(qubit) + X90(qubit) + Z90m(qubit) + X90(qubit) + Z(qubit) if Xonly else Z90m(qubit) + Y90(qubit)
    elif cliffNum == 20:
            #X+Y-Z -120 (equivalent to -X-Y+Z 120)
            return Z(qubit) + X90(qubit) + Z90(qubit) + X90(qubit) + Z90m(qubit) if Xonly else Z(qubit) + Y90m(qubit) + Z90(qubit)
    elif cliffNum == 21:
        #X+Y-Z 120
        return Z90(qubit) + X90(qubit) + Z90(qubit) + X90(qubit) + Z(qubit) if Xonly else Z90(qubit) + Y90m(qubit)
    elif cliffNum == 22:
        #-X+Y+Z -120 (equivalent to X-Y-Z 120)
        return X90(qubit) + Z90m(qubit) + X90(qubit) + Z90m(qubit) if Xonly else Y90(qubit) + Z90(qubit)
    elif cliffNum == 23:
        #-X+Y+Z 120
        return Z90m(qubit) + X90(qubit) + Z90(qubit) + X90(qubit) + Z(qubit) if Xonly else Z90m(qubit) + Y90m(qubit)
    else:
        raise ValueError('Clifford number must be between 0 and 23')

## two-qubit primitivies
def CNOT(source, target, **kwargs):
    # construct (source, target) channel and pull parameters from there
    channel = ChannelLibrary.EdgeFactory(source, target)
    channel.pulseParams['piAmp'] = channel.pulseParams['amp']
    p = X(channel, **kwargs)
    return p._replace(label="CNOT")


def flat_top_gaussian(chan,
                      riseFall,
                      length,
                      amp,
                      phase=0,
                      label="flat_top_gaussian"):
    """
    A constant pulse with rising and falling gaussian shape
    """
    return Utheta(chan, length=riseFall, amp=amp, phase=phase, shapeFun=PulseShapes.gaussOn, label=label+"_rise") + \
           Utheta(chan, length=length, amp=amp, phase=phase, shapeFun=PulseShapes.constant, label=label+"_top") + \
           Utheta(chan, length=riseFall, amp=amp, phase=phase, shapeFun=PulseShapes.gaussOff, label=label+"_fall")


def echoCR(controlQ,
           targetQ,
           amp=1,
           phase=0,
           length=200e-9,
           riseFall=20e-9,
           lastPi=True,
           echoQs=[]):
    """
    An echoed CR pulse.  Used for calibration of CR gate
    """
    CRchan = ChannelLibrary.EdgeFactory(controlQ, targetQ)
    if not CRchan.isforward(controlQ, targetQ):
        raise ValueError(
            'Could not find an edge with control qubit {0}'.format(controlQ))
    seq = [flat_top_gaussian(CRchan,
                             amp=amp,
                             riseFall=riseFall,
                             length=length,
                             phase=phase,
                             label="echoCR_first_half"),
                             reduce(operator.mul,
                                    [X(q) for q in echoQs + [controlQ]]),
           flat_top_gaussian(CRchan,
                             amp=amp,
                             riseFall=riseFall,
                             length=length,
                             phase=phase + np.pi,
                             label="echoCR_second_half")]
    if lastPi:
        seq += [reduce(operator.mul,
               [X(q) for q in echoQs + [controlQ]])]
    return seq


def ZX90_CR(controlQ, targetQ, echoQs=[], **kwargs):
    """
    A calibrated CR ZX90 pulse.  Uses 'amp' for the pulse amplitude, 'phase' for its phase (in deg).
    echoQs: list of qubits to decouple during the gate
    """
    CRchan = ChannelLibrary.EdgeFactory(controlQ, targetQ)
    params = overrideDefaults(CRchan, kwargs)
    return echoCR(controlQ,
                  targetQ,
                  echoQs = echoQs,
                  amp=params['amp'],
                  phase=params['phase'],
                  length=params['length'],
                  riseFall=params['riseFall'])


def CNOT_CR(controlQ, targetQ, **kwargs):
    edge = ChannelLibrary.EdgeFactory(controlQ, targetQ)

    if edge.isforward(controlQ, targetQ):
        # control and target for CNOT and CR match
        return ZX90_CR(controlQ, targetQ, **
                       kwargs) + [Z90m(controlQ) * X90m(targetQ)]
    else:
        # control and target for CNOT and CR are inverted
        return [Y90(controlQ) * Y90(targetQ),
                X(controlQ) * X(targetQ)] + \
                ZX90_CR(targetQ, controlQ, **kwargs) + \
               [Z90(targetQ),
                X90(controlQ) * Y90(targetQ),
                Y90m(controlQ) * X(targetQ)]


## Measurement operators
@_memoize
def MEAS(qubit, **kwargs):
    '''
    MEAS(q1) measures a qubit. Applies to the pulse with the label M-q1
    '''
    channelName = "M-" + qubit.label
    measChan = ChannelLibrary.MeasFactory(channelName)
    params = overrideDefaults(measChan, kwargs)
    if measChan.measType == 'autodyne':
        params['frequency'] = measChan.autodyneFreq
        params['baseShape'] = params.pop('shapeFun')
        params['shapeFun'] = PulseShapes.autodyne
    amp = params.pop('amp')
    ignoredStrParams = ['phase', 'frameChange']
    if 'amp' not in kwargs:
        ignoredStrParams.append('amp')
    return Pulse("MEAS", measChan, params, amp, 0.0, 0.0, ignoredStrParams)


#MEAS and ring-down time on one qubit, echo on every other
def MeasEcho(qM, qD, delay, piShift=None, phase=0, **kwargs):
    '''
    qM : qubit to be measured (single qubit)
    qD : qubits to be echoed (single qubit or tuple)
    delay : idle time after M-qM
    PiShift: relative shift of the echo pulse from the center of the pulse block (in s, to the right if >0)
    phase : rotation axis of the echo pulse
    '''
    channelName = "M-" + qM.label
    measChan = ChannelLibrary.MeasFactory(channelName)
    params = overrideDefaults(measChan, kwargs)

    if not isinstance(qD, tuple):
        qD = (qD, )
    measChan = ChannelLibrary.MeasFactory('M-%s' % qM.label)
    if piShift:
        if piShift > 0:
            measEcho = align(
                (MEAS(qM, **params) + TAPulse('Id', measChan, delay, 0)) *
                reduce(operator.mul,
                       [Id(q, piShift) + U(q, phase=phase) for q in qD]))
        elif piShift < 0:
            measEcho = align(
                (MEAS(qM, **params) + TAPulse('Id', measChan, delay, 0)) *
                reduce(operator.mul,
                       [U(q, phase=phase) + Id(q, -piShift) for q in qD]))
    else:
        measEcho = align((MEAS(qM, **params) + TAPulse('Id', measChan, delay, 0)) *
                         reduce(operator.mul, [U(q, phase=phase) for q in qD]))
    measEcho.label = 'MEAS'  #to generate the digitizer trigger
    return measEcho


# Gating/blanking pulse primitives
def BLANK(chan, length):
    return TAPulse("BLANK", chan.gateChan, length, 1, 0, 0)
