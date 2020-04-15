"""
Manipulating Cliffords.  Mainly for RB purposes.
"""
import numpy as np
from scipy.linalg import expm
from numpy import pi
from itertools import product
from random import choice
import operator
from functools import reduce

from .tools.clifford_tools import *
from .tools.euler_angles import xyx_angles
from .PulsePrimitives import *

###
### Single Qubit Cliffords
###

def generator_pulse(G):
    """
	A function that returns the pulse corresponding to a generator
	Randomly chooses between the -p and -m versions for the 180's
	"""
    generatorPulses = {0: (Id, ),
                       1: (X90, ),
                       3: (X90m, ),
                       4: (Y90, ),
                       6: (Y90m, ),
                       2: (X, Xm),
                       5: (Y, Ym)}
    return choice(generatorPulses[G])

def StdClifford(qubit, cliffNum):
    """

    The set of 24 clifford pulses from the Pauli group generators.
    Note that this will randomly select between +/- X and +- Y pulses.

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    cliffNum : the zero-indexed Clifford number

    Returns
    -------
    pulse object
    """
    genSeq = generatorSeqs[choice(C1Seqs[cliffNum])]
    return reduce(operator.add, [generator_pulse(g)(qubit) for g in genSeq])

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

    #Figure out the approximate nutation frequency calibration from the X180 and the sampling_rate
    Xp = X(qubit)
    xpulse = Xp.amp * Xp.shape
    nutFreq = 0.5 / (sum(xpulse) / qubit.phys_chan.sampling_rate)

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

def get_DiAC_phases(cliffNum):
    """

    Returns the phases (in multiples of pi) of the three Z gates dressing the two X90
    pulses comprising the DiAC pulse correspoding to cliffNum
    e.g., get_DiAC_phases(1) returns a=0, b=1, c=1, in
    Ztheta(a) + X90 + Ztheta(b) + X90 + Ztheta(c) = Id
    """
    DiAC_table = [
    [0, 1, 1],
    [0.5, -0.5, 0.5],
    [0, 0, 0],
    [0.5, 0.5, 0.5],
    [0, -0.5, 1],
    [0, 0, 1],
    [0, 0.5, 1],
    [0, 1, -0.5],
    [0, 1, 0],
    [0, 1, 0.5],
    [0, 0, 0.5],
    [0, 0, -0.5],
    [1, -0.5, 1],
    [1, 0.5, 1],
    [0.5, -0.5, -0.5],
    [0.5, 0.5, -0.5],
    [0.5, -0.5, 1],
    [1, -0.5, -0.5],
    [0, 0.5, -0.5],
    [-0.5, -0.5, 1],
    [1, 0.5, -0.5],
    [0.5, 0.5, 1],
    [0, -0.5, -0.5],
    [-0.5, 0.5, 1]]
    return DiAC_table[cliffNum]

def DiAC(qubit, cliffNum, compiled = True):
    """

    The set of 24 Diatomic Clifford single qubit pulses. Each pulse is decomposed
    as Ztheta(a) + X90 + Ztheta(b) + X90 + Ztheta(c) if compiled = False,
    uses also Y90, Y90m and shorter sequences if compiled = True

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    cliffNum : the zero-indexed Clifford number

    Returns
    -------
    pulse object
    """
    #Now a big else if chain for to get the specific Clifford
    if not compiled:
        DiAC_phases = get_DiAC_phases(cliffNum)
        return Ztheta(qubit, angle = DiAC_phases[0]*np.pi) + X90(qubit) + Ztheta(qubit, angle = DiAC_phases[1]*np.pi) + \
        X90(qubit) + Ztheta(qubit, angle = DiAC_phases[2]*np.pi)
    else:
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
            return Y90(qubit)
        elif cliffNum == 5:
            #Y180
            return Y90(qubit) + Y90(qubit)
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
            return Y90(qubit) + Y90(qubit) + Z90m(qubit)
        elif cliffNum == 11:
            #X-Y 180
            return Y90(qubit) + Y90(qubit) + Z90(qubit)
        elif cliffNum == 12:
            #X+Z 180(Hadamard)
            return Z(qubit) + Y90(qubit)
        elif cliffNum == 13:
            #X-Z 180
            return Z(qubit) + Y90m(qubit)
        elif cliffNum == 14:
            #Y+Z 180
            return Z90(qubit) + Y90(qubit) + Z90(qubit)
        elif cliffNum == 15:
            #Y-Z 180
            return Z90(qubit) + Y90m(qubit) + Z90(qubit)
        elif cliffNum == 16:
                #X+Y+Z -120 (equivalent to -X-Y-Z 120)
            return Z90(qubit) + Y90(qubit)
        elif cliffNum == 17:
            #X+Y+Z 120
            return Z(qubit) + Y90(qubit) + Z90(qubit)
        elif cliffNum == 18:
            #X-Y+Z 120 (equivalent to -X+Y-Z 120)
            return Y90m(qubit) + Z90(qubit)
        elif cliffNum == 19:
            #X-Y+Z -120
            return Z90m(qubit) + Y90(qubit)
        elif cliffNum == 20:
            #X+Y-Z -120 (equivalent to -X-Y+Z 120)
            return Z(qubit) + Y90m(qubit) + Z90(qubit)
        elif cliffNum == 21:
            #X+Y-Z 120
            return Z90(qubit) + Y90m(qubit)
        elif cliffNum == 22:
            #-X+Y+Z -120 (equivalent to X-Y-Z 120)
            return Y90(qubit) + Z90(qubit)
        elif cliffNum == 23:
            #-X+Y+Z 120
            return Z90m(qubit) + Y90m(qubit)
        else:
            raise ValueError('Clifford number must be between 0 and 23')

def XYXClifford(qubit, cliff_num):
    """
    The set of 24 Diatomic Clifford single qubit pulses. Each pulse is decomposed
    as Rx(α)Ry(β)Rx(γ).

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    cliffNum : the zero-indexed Clifford number

    Returns
    -------
    pulse object
    """
    α, β, γ = xyx_angles(C1[cliff_num])

    p1 =  Id(qubit) if np.isclose(γ, 0.0) else Xtheta(qubit, angle=γ)
    p2 =  Id(qubit) if np.isclose(β, 0.0) else Ytheta(qubit, angle=β)
    p3 =  Id(qubit) if np.isclose(α, 0.0) else Xtheta(qubit, angle=α)

    return p1 + p2 + p3

###
### Two qubit Cliffords
###

clifford_map = {}
clifford_map['STD'] = StdClifford
clifford_map['DIAC'] = DiAC 
clifford_map['AC'] = AC 
clifford_map['XYX'] = XYXClifford

def Cx2(c1, c2, q1, q2, kind='std'):
    """
    Helper function to create pulse block for a pair of single-qubit Cliffords
    """
    
    clifford_fun = clifford_map[kind.upper()]
    seq1 = clifford_fun(q1, c1)
    seq2 = clifford_fun(q2, c2)

    #Create the pulse block
    return seq1 * seq2

def entangling_seq(gate, q1, q2):
    """
    Helper function to create the entangling gate sequence
    """
    if gate == "CNOT":
        return ZX90_CR(q2, q1)
    elif gate == "iSWAP":
        return [ZX90_CR(q2, q1) , Y90m(q1) * Y90m(q2), ZX90_CR(q2, q1)]
    elif gate == "SWAP":
        return [ZX90_CR(q2, q1), Y90m(q1) * Y90m(q2), ZX90_CR(
            q2, q1), (X90(q1) + Y90m(q1)) * X90(q2), ZX90_CR(q2, q1)]

def TwoQubitClifford(q1, q2, cliffNum, kind='std'):

    if kind.upper() not in clifford_map.keys():
        raise ValueError(f"Unknown clifford type: must be one of {clifford.map.keys()}.")

    c = C2Seqs[cliffNum]
    seq = [Cx2(c[0][0], c[0][1], q1, q2, kind=kind)]
    if c[1]:
        seq += entangling_seq(c[1], q1, q2)
    if c[2]:
        seq += [Cx2(c[2][0], c[2][1], q1, q2, kind=kind)]
    return seq