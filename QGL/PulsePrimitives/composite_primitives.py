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
from .common_primitives import Id, Z, X90, X90m, Y90, Y90m, Z90, Z90m, Utheta, Xtheta, Ytheta, Ztheta, U90, U, arb_axis_drag
from .. import config
if config.pulse_primitives_lib == 'standard':
    from .standard_primitives import X, Xm, Y, Ym, AC
elif config.pulse_primitives_lib == 'all90':
    from .all90_primitives import X, Xm, Y, Ym, AC
else:
    raise Exception("Invalid pulse library")

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
           lastPi=True):
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
                             label="echoCR_first_half"), X(controlQ),
           flat_top_gaussian(CRchan,
                             amp=amp,
                             riseFall=riseFall,
                             length=length,
                             phase=phase + np.pi,
                             label="echoCR_second_half")]
    if lastPi:
        seq += [X(controlQ)]
    return seq


def ZX90_CR(controlQ, targetQ, **kwargs):
    """
    A calibrated CR ZX90 pulse.  Uses 'amp' for the pulse amplitude, 'phase' for its phase (in deg).
    """
    CRchan = ChannelLibrary.EdgeFactory(controlQ, targetQ)
    params = overrideDefaults(CRchan, kwargs)
    return echoCR(controlQ,
                  targetQ,
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
    measChan = ChannelLibrary.MeasFactory('M-%s' % qM.label)
    if piShift:
        if piShift > 0:
            measEcho = align(
                (MEAS(qM) + TAPulse('Id', measChan, delay, 0)) *
                reduce(operator.mul,
                       [Id(q, piShift) + U(q, phase=phase) for q in qD]))
        elif piShift < 0:
            measEcho = align(
                (MEAS(qM) + TAPulse('Id', measChan, delay, 0)) *
                reduce(operator.mul,
                       [U(q, phase=phase) + Id(q, -piShift) for q in qD]))
    else:
        measEcho = align((MEAS(qM) + TAPulse('Id', measChan, delay, 0)) *
                         reduce(operator.mul, [U(q, phase=phase) for q in qD]))
    measEcho.label = 'MEAS'  #to generate the digitizer trigger
    return measEcho
