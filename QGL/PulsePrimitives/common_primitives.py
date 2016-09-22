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
from .. import PulseShapes
from .. import Channels
from .. import ChannelLibrary
import operator

from math import pi, sin, cos, acos, sqrt
import numpy as np
from ..PulseSequencer import Pulse, TAPulse, align
from functools import wraps, reduce
from .helpers import overrideDefaults, _memoize, clear_pulse_cache

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

@_memoize
def Z(qubit, **kwargs):
    return Ztheta(qubit, pi, label="Z", **kwargs)


@_memoize
def Z90(qubit, **kwargs):
    return Ztheta(qubit, pi / 2, label="Z90", **kwargs)


@_memoize
def Z90m(qubit, **kwargs):
    return Ztheta(qubit, -pi / 2, label="Z90m", **kwargs)


# 90/180 degree rotations with control over the rotation axis
def U90(qubit, phase=0, **kwargs):
    """ A generic 90 degree rotation with variable phase. """
    if "label" not in kwargs:
        kwargs["label"] = "U90"
    return Utheta(qubit,
                  qubit.pulseParams['pi2Amp'],
                  phase,
                  ignoredStrParams=['amp'],
                  **kwargs)


def U(qubit, phase=0, **kwargs):
    """ A generic 180 degree rotation with variable phase.  """
    if "label" not in kwargs:
        kwargs["label"] = "U"
    return Utheta(qubit,
                  qubit.pulseParams['piAmp'],
                  phase,
                  ignoredStrParams=['amp'],
                  **kwargs)


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

# Gating/blanking pulse primitives
def BLANK(chan, length):
    return TAPulse("BLANK", chan.gateChan, length, 1, 0, 0)
