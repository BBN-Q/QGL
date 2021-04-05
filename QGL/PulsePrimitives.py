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
from . import ChannelLibraries
from . import config
import operator

from math import pi, sin, cos, acos, sqrt
import numpy as np
from .PulseSequencer import Pulse, TAPulse, CompoundGate, align
from functools import wraps, reduce


def overrideDefaults(chan, updateParams):
    '''Helper function to update any parameters passed in and fill in the defaults otherwise.'''
    # The default parameter list depends on the channel type so pull out of channel
    # Then update passed values
    paramDict = chan.pulse_params.copy()
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
           angle=0,
           phase=0,
           label='Utheta',
           ignoredStrParams=[],
           **kwargs):
    '''
    A generic rotation with variable angle and phase.
    This primitive needs to convert the requested rotation angle into a Pulse
    amplitude. If the user would like to specify a specific amplitude, he/she
    can specify it directly via the 'amp' keyword argument.
    '''
    params = overrideDefaults(qubit, kwargs)
    # amp and phase are now pulse parameters rather than shape parameters
    if "amp" in params:
        del params["amp"]
    if "phase" in params:
        del params["phase"]
    if "frequency" in kwargs:
        frequency = kwargs["frequency"]
    else:
        frequency = None
    # allow override of angle -> amplitude lookup if the user provides an "amp"
    # keyword argument
    if "amp" in kwargs:
        amp = kwargs["amp"]
    else:
        # construct an angle -> amplitude lookup table
        # TODO should this live in the Channel object instead?
        angle2amp = {
            pi    :  qubit.pulse_params['piAmp'],
            -pi   : -qubit.pulse_params['piAmp'],
            pi/2  :  qubit.pulse_params['pi2Amp'],
            -pi/2 : -qubit.pulse_params['pi2Amp'],
        }
        if angle in angle2amp:
            amp = angle2amp[angle]
        else:
            # linearly scale based upon the 'pi/2' amplitude
            amp  = (angle / (pi/2)) * qubit.pulse_params['pi2Amp']
    return Pulse(label, qubit, params, amp, phase, 0.0, ignoredStrParams, frequency=frequency)


# generic pulses around X, Y, and Z axes
def Xtheta(qubit, angle=0, label='Xtheta', ignoredStrParams=None, **kwargs):
    '''  A generic X rotation with a variable rotation angle  '''
    if ignoredStrParams is None:
        ignoredStrParams = ['phase', 'frameChange']
    else:
        ignoredStrParams += ['phase', 'frameChange']
    return Utheta(qubit, angle, 0, label, ignoredStrParams, **kwargs)


def Ytheta(qubit, angle=0, label='Ytheta', ignoredStrParams=None, **kwargs):
    ''' A generic Y rotation with a variable rotation angle '''
    if ignoredStrParams is None:
        ignoredStrParams = ['phase', 'frameChange']
    else:
        ignoredStrParams += ['phase', 'frameChange']
    return Utheta(qubit, angle, pi/2, label, ignoredStrParams, **kwargs)


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
                  pi/2,
                  label="X90",
                  ignoredStrParams=['amp'],
                  **kwargs)

@_memoize
def X90m(qubit, **kwargs):
    return Xtheta(qubit,
                  -pi/2,
                  label="X90m",
                  ignoredStrParams=['amp'],
                  **kwargs)

@_memoize
def Y90(qubit, **kwargs):
    return Ytheta(qubit,
                  pi/2,
                  label="Y90",
                  ignoredStrParams=['amp'],
                  **kwargs)

@_memoize
def Y90m(qubit, **kwargs):
    return Ytheta(qubit,
                  -pi/2,
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
                  pi/2,
                  phase,
                  ignoredStrParams=['amp'],
                  **kwargs)

@_memoize
def X(qubit, **kwargs):
    if config.pulse_primitives_lib == 'standard':
        return Xtheta(qubit,
                      pi,
                      label="X",
                      ignoredStrParams=['amp'],
                      **kwargs)
    elif config.pulse_primitives_lib == 'all90':
        return X90(qubit, **kwargs) + X90(qubit, **kwargs)
    else:
        raise Exception("Invalid pulse_primitives_lib. Must be standard or all90")

@_memoize
def Xm(qubit, **kwargs):
    if config.pulse_primitives_lib == 'standard':
        return Xtheta(qubit,
                      -pi,
                      label="Xm",
                      ignoredStrParams=['amp'],
                      **kwargs)
    elif config.pulse_primitives_lib == 'all90':
        return X90m(qubit, **kwargs) + X90m(qubit, **kwargs)
    else:
        raise Exception("Invalid pulse_primitives_lib. Must be standard or all90")

@_memoize
def Y(qubit, **kwargs):
    if config.pulse_primitives_lib == 'standard':
        return Ytheta(qubit,
                      pi,
                      label="Y",
                      ignoredStrParams=['amp'],
                      **kwargs)
    elif config.pulse_primitives_lib == 'all90':
        return Y90(qubit, **kwargs) + Y90(qubit, **kwargs)
    else:
        raise Exception("Invalid pulse_primitives_lib. Must be standard or all90")

@_memoize
def Ym(qubit, **kwargs):
    if config.pulse_primitives_lib == 'standard':
        return Ytheta(qubit,
                      -pi,
                      label="Ym",
                      ignoredStrParams=['amp'],
                      **kwargs)
    elif config.pulse_primitives_lib == 'all90':
        return Y90m(qubit, **kwargs) + Y90m(qubit, **kwargs)
    else:
        raise Exception("Invalid pulse_primitives_lib. Must be standard or all90")

@_memoize
def U(qubit, phase=0, **kwargs):
    """ A generic 180 degree rotation with variable phase.  """
    if "label" not in kwargs:
        kwargs["label"] = "U"
    if config.pulse_primitives_lib == 'standard':
        return Utheta(qubit,
                      pi,
                      phase,
                      ignoredStrParams=['amp'],
                      **kwargs)
    elif config.pulse_primitives_lib == 'all90':
        return U90(qubit, phase, *kwargs) + U90(qubit, phase, *kwargs)
    else:
        raise Exception("Invalid pulse_primitives_lib. Must be standard or all90")

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
        sampRate = qubit.phys_chan.sampling_rate

        #Start from a gaussian shaped pulse
        gaussPulse = PulseShapes.gaussian(amp=1,
                                          sampling_rate=sampRate,
                                          **params).real

        #Scale to achieve to the desired rotation
        calScale = (rotAngle / 2 / pi) * sampRate / sum(gaussPulse)

        #Calculate the phase ramp steps to achieve the desired Z component to the rotation axis
        phaseSteps = -2*pi * cos(polarAngle) * calScale * gaussPulse / sampRate

        #Calculate Z DRAG correction to phase steps
        #beta is a conversion between XY drag scaling and Z drag scaling
        beta = params['drag_scaling'] / sampRate
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
    params['shape_fun'] = PulseShapes.arb_axis_drag
    return Pulse(kwargs["label"] if "label" in kwargs else "ArbAxis", qubit,
                 params, 1.0, aziAngle, frameChange)


def DiatomicPulse(qubit, a, b, c):
  return (Ztheta(qubit, angle=c) + X90(qubit) + 
          Ztheta(qubit, angle=b) + X90(qubit) + 
          Ztheta(qubit, angle=a))

def ZYZPulse(qubit, a, b, c):
  Ypulse = Id(qubit) if np.isclose(b, 0.0) else Ytheta(qubit, angle=b)
  return Ztheta(qubit, angle=c)+Ypulse+Ztheta(qubit, angle=a)

def XYXPulse(qubit, α, β, γ):
  p1 =  Id(qubit) if np.isclose(γ, 0.0) else Xtheta(qubit, angle=γ)
  p2 =  Id(qubit) if np.isclose(β, 0.0) else Ytheta(qubit, angle=β)
  p3 =  Id(qubit) if np.isclose(α, 0.0) else Xtheta(qubit, angle=α)
  return p1+p2+p3

## two-qubit primitivies

# helper used by echoCR
def flat_top_gaussian(chan,
                      riseFall,
                      length,
                      amp,
                      phase=0,
                      label="flat_top_gaussian"):
    """
    A constant pulse with rising and falling gaussian shape
    """
    # Note: use_cos is a hack that should definitely be fixed more elegantly later
    use_cos = False
    if riseFall == 0:
        p = Utheta(chan, length=length, amp=amp, phase=phase, shape_fun=PulseShapes.constant, label=label+"_top")
    else:
        if use_cos:
            #print('Using cosine ZX90')
            p =  Utheta(chan, length=riseFall, amp=amp, phase=phase, shape_fun=PulseShapes.cosOn, label=label+"_rise") + \
                 Utheta(chan, length=length, amp=amp, phase=phase, shape_fun=PulseShapes.constant, label=label+"_top") + \
                 Utheta(chan, length=riseFall, amp=amp, phase=phase, shape_fun=PulseShapes.cosOff, label=label+"_fall")
        else:
            p =  Utheta(chan, length=riseFall, amp=amp, phase=phase, shape_fun=PulseShapes.gaussOn, label=label+"_rise") + \
                 Utheta(chan, length=length, amp=amp, phase=phase, shape_fun=PulseShapes.constant, label=label+"_top") + \
                 Utheta(chan, length=riseFall, amp=amp, phase=phase, shape_fun=PulseShapes.gaussOff, label=label+"_fall")
    return p._replace(label=label)


def echoCR(controlQ,
           targetQ,
           amp=1,
           phase=0,
           length=200e-9,
           riseFall=20e-9,
           lastPi=True, canc_amp = 0.0, canc_phase=np.pi/2):
    """
    An echoed CR pulse.  Used for calibration of CR gate
    """
    CRchan = ChannelLibraries.EdgeFactory(controlQ, targetQ)
    if not CRchan.isforward(controlQ, targetQ):
        raise ValueError(
            'Could not find an edge with control qubit {0}'.format(controlQ))

    seq = [flat_top_gaussian(CRchan,
                             amp=amp,
                             riseFall=riseFall,
                             length=length,
                             phase=phase,
                             label="echoCR_first_half"),
           X(controlQ),
           flat_top_gaussian(CRchan,
                             amp=amp,
                             riseFall=riseFall,
                             length=length,
                             phase=phase + np.pi,
                             label="echoCR_second_half")]
    if canc_amp != 0:
        seq[0]*=flat_top_gaussian(targetQ,
                            amp=canc_amp,
                            riseFall=riseFall,
                            length=length,
                            phase=canc_phase,
                            label="cancCR_first_half")
        seq[2]*=flat_top_gaussian(targetQ,
                            amp=canc_amp,
                            riseFall=riseFall,
                            length=length,
                            phase=canc_phase + np.pi,
                            label="cancCR_second_half")
    if lastPi:
        seq += [X(controlQ)]
    return CompoundGate(seq)


def ZX90_CR(controlQ, targetQ, **kwargs):
    """
    A calibrated CR ZX90 pulse.  Uses 'amp' for the pulse amplitude, 'phase' for its phase (in deg).
    """
    CRchan = ChannelLibraries.EdgeFactory(controlQ, targetQ)
    params = overrideDefaults(CRchan, kwargs)
    canc_amp = params['canc_amp'] if 'canc_amp' in params else 0
    canc_phase = params['canc_phase'] if 'canc_phase' in params else 0
    return echoCR(controlQ,
                  targetQ,
                  amp=params['amp'],
                  phase=params['phase'],
                  length=params['length'],
                  riseFall=params['riseFall'],
                  canc_phase=canc_phase,
                  canc_amp=canc_amp)


def CNOT_CR(controlQ, targetQ, **kwargs):
    edge = ChannelLibraries.EdgeFactory(controlQ, targetQ)

    if edge.isforward(controlQ, targetQ):
        # control and target for CNOT and CR match
        seq = ZX90_CR(controlQ, targetQ, **kwargs).seq + \
            [Z90m(controlQ) * X90m(targetQ)]
    else:
        # control and target for CNOT and CR are inverted
        seq = [Y90(controlQ) * Y90(targetQ),
                X(controlQ) * X(targetQ)] + \
                ZX90_CR(targetQ, controlQ, **kwargs).seq + \
               [Z90(targetQ),
                X90(controlQ) * Y90(targetQ),
                Y90m(controlQ) * X(targetQ)]
    return CompoundGate(seq)

def CNOT_simple(source, target, **kwargs):
    # construct (source, target) channel and pull parameters from there
    channel = ChannelLibraries.EdgeFactory(source, target)
    channel.pulse_params['piAmp'] = channel.pulse_params['amp']
    # add "pi2Amp" too so that Utheta can construct its angle2amp lookup table
    channel.pulse_params['pi2Amp'] = channel.pulse_params['amp'] / 2
    p = X(channel, **kwargs)
    return p._replace(label="CNOT")

@_memoize
def CNOT(source, target, **kwargs):
    channel = ChannelLibraries.EdgeFactory(source, target)
    if hasattr(channel, 'cnot_impl') and channel.cnot_impl:
        cnot_impl_name = channel.cnot_impl
    else:
        cnot_impl_name = config.cnot_implementation
    cnot_impl = globals()[cnot_impl_name]
    return cnot_impl(source, target, **kwargs)

# The worker method for MEAS and MEASA
def _MEAS(qubit, **kwargs):

    '''
    _MEAS(q1) implements both MEAS and MEASA, but because of
    the way memoize works, we want to distinguish them and
    memoize them separately.  (this may change if the way
    memoization works is changed)

    TODO: this is annoying because measuring the same qubit
    but storing the result to two different addresses creates
    two distinct pulses, unless we also memoize the waveforms
    themselves.
    '''

    channelName = "M-" + qubit.label
    measChan = ChannelLibraries.MeasFactory(channelName)
    params = overrideDefaults(measChan, kwargs)
    if measChan.meas_type == 'autodyne':
        params['frequency'] = measChan.autodyne_freq
        params['baseShape'] = params.pop('shape_fun')
        params['shape_fun'] = PulseShapes.autodyne
    amp = params.pop('amp')
    ignoredStrParams = ['phase', 'frameChange']
    if 'amp' not in kwargs:
        ignoredStrParams.append('amp')
    meas_label = "MEAS_no_trig" if 'dig_trig' in kwargs and not kwargs['dig_trig'] else "MEAS"
    if 'maddr' in kwargs:
        maddr = kwargs['maddr']
    else:
        maddr = (-1, 0)

    return Pulse(meas_label, measChan, params,
            amp, 0.0, 0.0, ignoredStrParams, maddr=maddr)

## Measurement operators
@_memoize
def MEAS(qubit, **kwargs):
    '''
    MEAS(q1) measures a qubit (applying the pulse with the label M-q1)

    This is the 'traditional' measurement; it does not require
    a measurement address, and will gripe if one is provided in
    the kwargs

    Note that in the future there may be hardware that requires a
    measurement address (even if it's just a placeholder)
    '''

    if 'maddr' in kwargs:
        raise ValueError('MEAS must not have a maddr kwarg')

    return _MEAS(qubit, **kwargs)

@_memoize
def MEASA(qubit, maddr, **kwargs):
    '''
    MEASA(q1) measures a qubit (applying the pulse with the label M-q1)
    and stores the result in the given address in "measurement memory".

    Note that a failure will occur if the hardware does not support
    measurement memory and a MEASA is requested.

    There may be special values for maddr that change the behavior
    of this operation, but they have not been specified yet.
    '''

    return _MEAS(qubit, maddr=maddr, **kwargs)

#MEAS and ring-down time on one qubit, echo on every other
def MeasEcho(qM, qD, delay, piShift=None, phase=0):
    '''
    qM : qubit to be measured (single qubit)
    qD : qubits to be echoed (single qubit or tuple)
    delay : idle time after M-qM
    PiShift: relative shift of the echo pulse from the center of the pulse block (in s, to the right if >0)
    phase : rotation axis of the echo pulse
    '''
    if not isinstance(qD, tuple):
        qD = (qD, )
    measChan = ChannelLibraries.MeasFactory('M-%s' % qM.label)
    if piShift:
        if piShift > 0:
            measEcho = (MEAS(qM) + TAPulse('Id', measChan, delay, 0)) * reduce(operator.mul, [Id(q, piShift) + U(q, phase=phase) for q in qD])
        elif piShift < 0:
            measEcho = (MEAS(qM) + TAPulse('Id', measChan, delay, 0)) * reduce(operator.mul, [U(q, phase=phase) + Id(q, -piShift) for q in qD])
    else:
        measEcho = (MEAS(qM) + TAPulse('Id', measChan, delay, 0)) * reduce(operator.mul, [U(q, phase=phase) for q in qD])
    measEcho.label = 'MEAS'  #to generate the digitizer trigger
    return measEcho

# Gating/blanking pulse primitives
def BLANK(chan, length):
    return TAPulse("BLANK", chan.gate_chan, length, 1, 0, 0)

def TRIG(marker_chan, length):
    '''TRIG(marker_chan, length) generates a trigger output of amplitude 1 on
    a LogicalMarkerChannel.
    '''
    if not isinstance(marker_chan, Channels.LogicalMarkerChannel):
        raise ValueError("TRIG pulses can only be generated on LogicalMarkerChannels.")
    return TAPulse("TRIG", marker_chan, length, 1.0, 0., 0.)
