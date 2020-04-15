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
from .PulsePrimitives import *

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

def Cx2(c1, c2, q1, q2):
    """
	Helper function to create pulse block for a pair of single-qubit Cliffords
	"""
    #Create list of pulse objects on the qubits
    seq1 = clifford_seq(c1, q1)
    seq2 = clifford_seq(c2, q2)

    #Create the pulse block
    return reduce(operator.add, seq1) * reduce(operator.add, seq2)


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


def entangling_mat(gate):
    """
	Helper function to create the entangling gate matrix
	"""
    echoCR = expm(1j * pi / 4 * np.kron(pX, pZ))
    if gate == "CNOT":
        return echoCR
    elif gate == "iSWAP":
        return reduce(lambda x, y: np.dot(y, x),
                      [echoCR, np.kron(C1[6], C1[6]), echoCR])
    elif gate == "SWAP":
        return reduce(lambda x, y: np.dot(y, x),
                      [echoCR, np.kron(C1[6], C1[6]), echoCR, np.kron(
                          np.dot(C1[6], C1[1]), C1[1]), echoCR])
    else:
        raise ValueError("Entangling gate must be one of: CNOT, iSWAP, SWAP.")


def clifford_seq(c, q1, q2=None):
    """
	Return a sequence of pulses that implements a clifford C
	"""
    #If qubit2 not defined assume 1 qubit
    if not q2:
        genSeq = generatorSeqs[choice(C1Seqs[c])]
        return [generator_pulse(g)(q1) for g in genSeq]
    else:
        #Look up the sequence for the integer
        c = C2Seqs[c]
        seq = [Cx2(c[0][0], c[0][1], q1, q2)]
        if c[1]:
            seq += entangling_seq(c[1], q1, q2)
        if c[2]:
            seq += [Cx2(c[2][0], c[2][1], q1, q2)]
        return seq
