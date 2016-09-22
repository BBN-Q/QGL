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

from .helpers import overrideDefaults, _memoize, clear_pulse_cache
from .common_primitives import *
#Setup the default 90/180 rotations
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
