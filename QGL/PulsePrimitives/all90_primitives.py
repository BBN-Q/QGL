'''
Copyright 2016 Raytheon BBN Technologies

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

def AC(qubit, cliffNum, Xonly = False):
    """

    The set of 24 Diatomic Clifford single qubit pulses.

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
        return X90(qubit) + X90(qubit)
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
