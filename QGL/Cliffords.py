"""
Manipulating Cliffords.  Mainly for RB purposes.
"""
import numpy as np
from scipy.linalg import expm
from numpy import pi
from itertools import product
from random import choice
import operator

from PulsePrimitives import *

#Single qubit paulis
pX = np.array([[0, 1],[1, 0]], dtype=np.complex128)
pY = np.array([[0, -1j],[1j, 0]], dtype=np.complex128)
pZ = np.array([[1, 0],[0, -1]], dtype=np.complex128)
pI = np.eye(2, dtype=np.complex128)

def pauli_mats(n):
	"""
	Return a list of n-qubit Paulis as numpy array.
	"""
	assert n > 0, "You need at least 1 qubit!"
	if n == 1:
		return [pI, pX, pY, pZ]
	else:
		paulis = pauli_mats(n-1)
		return [np.kron(p1, p2) for p1,p2 in product([pI, pX, pY, pZ], paulis)]

#Basis single-qubit Cliffords with an arbitrary enumeration
C1 = {}
C1[0] = pI
C1[1] = expm(-1j*(pi/4)*pX)
C1[2] = expm(-2j*(pi/4)*pX)
C1[3] = expm(-3j*(pi/4)*pX)
C1[4] = expm(-1j*(pi/4)*pY)
C1[5] = expm(-2j*(pi/4)*pY)
C1[6] = expm(-3j*(pi/4)*pY)
C1[7] = expm(-1j*(pi/4)*pZ)
C1[8] = expm(-2j*(pi/4)*pZ)
C1[9] = expm(-3j*(pi/4)*pZ)
C1[10] = expm(-1j*(pi/2)*(1/np.sqrt(2))*(pX+pY))
C1[11] = expm(-1j*(pi/2)*(1/np.sqrt(2))*(pX-pY))
C1[12] = expm(-1j*(pi/2)*(1/np.sqrt(2))*(pX+pZ))
C1[13] = expm(-1j*(pi/2)*(1/np.sqrt(2))*(pX-pZ))
C1[14] = expm(-1j*(pi/2)*(1/np.sqrt(2))*(pY+pZ))
C1[15] = expm(-1j*(pi/2)*(1/np.sqrt(2))*(pY-pZ))
C1[16] = expm(-1j*(pi/3)*(1/np.sqrt(3))*(pX+pY+pZ))
C1[17] = expm(-2j*(pi/3)*(1/np.sqrt(3))*(pX+pY+pZ))
C1[18] = expm(-1j*(pi/3)*(1/np.sqrt(3))*(pX-pY+pZ))
C1[19] = expm(-2j*(pi/3)*(1/np.sqrt(3))*(pX-pY+pZ))
C1[20] = expm(-1j*(pi/3)*(1/np.sqrt(3))*(pX+pY-pZ))
C1[21] = expm(-2j*(pi/3)*(1/np.sqrt(3))*(pX+pY-pZ))
C1[22] = expm(-1j*(pi/3)*(1/np.sqrt(3))*(-pX+pY+pZ))
C1[23] = expm(-2j*(pi/3)*(1/np.sqrt(3))*(-pX+pY+pZ))

#A little memoize decorator
def memoize(function):
	cache = {}
	def decorated(*args):
		if args not in cache:
			cache[args] = function(*args)
		return cache[args]
	return decorated

@memoize
def clifford_multiply(c1, c2):
    '''
    Multiplication table for single qubit cliffords.  Note this assumes c1 is applied first.
    i.e.  clifford_multiply(c1, c2) calculates c2*c1
    '''
    tmpMult = np.dot(C1[c2],C1[c1])
    checkArray = np.array([np.abs(np.trace(np.dot(tmpMult.transpose().conj(),C1[x]))) for x in range(24)])
    return checkArray.argmax()

#We can usually (without atomic Cliffords) only apply a subset of the single-qubit Cliffords
#i.e. the pulses that we can apply: Id, X90, X90m, Y90, Y90m, X, Y
generatorPulses = [0, 1, 3, 4, 6, 2, 5]

def generator_pulse(G):
	"""
	A function that returns the pulse corresponding to a generator 
	Randomly chooses between the -p and -m versions for the 180's
	"""
	generatorPulses = {0:(Id,), 1:(X90,), 3:(X90m,), 4:(Y90,), 6:(Y90m,), 2:(X,Xm), 5:(Y,Ym)}
	return choice(generatorPulses[G])

#Get all combinations of generator sequences up to length three
generatorSeqs = [x for x in product(generatorPulses,repeat=1)] + \
                [x for x in product(generatorPulses,repeat=2)] + \
				[x for x in product(generatorPulses,repeat=3)]

#Find the effective unitary for each generator sequence
reducedSeqs = np.array([reduce(clifford_multiply,x) for x in generatorSeqs ])

#Pick first generator sequence (and thus shortest) that gives each Clifford and then
#also add all those that have the same length

#First for each of the 24 single-qubit Cliffords find which sequences create them
allC1Seqs = [np.nonzero(reducedSeqs==x)[0] for x in range(24)]
#And the length of the first one for all 24
minSeqLengths = [len(generatorSeqs[seqs[0]]) for seqs in allC1Seqs]
#Now pull out all those that are the same length as the first one
C1Seqs = []
for minLength,seqs in zip(minSeqLengths, allC1Seqs):
	C1Seqs.append(filter(lambda seqNum: len(generatorSeqs[seqNum])==minLength, seqs))


C2Seqs = []

"""
The IBM paper has the Sgroup (rotation n*(pi/3) rotations about the X+Y+Z axis)
Sgroup = [C[0], C[16], C[17]]

The two qubit Cliffords can be written down as the product of 
1. A choice of one of 24^2 C \otimes C single-qubit Cliffords
2. Optionally an entangling gate from CNOT, iSWAP and SWAP
3. Optional one of 9 S \otimes S gate

Therefore, we'll enumerate the two-qubit Clifford as a three tuple ((c1,c2), Entangling, (s1,s2))

"""

#1. All pairs of single-qubit Cliffords
for c1, c2 in product(range(24), repeat=2):
	C2Seqs.append( ((c1,c2), None, None) )

#2. The CNOT-like class, replacing the CNOT with a echoCR
#TODO: sort out whether we need to explicitly encorporate the single qubit rotations into the trailing S gates
# The leading single-qubit Cliffords are fully sampled so they should be fine

for (c1,c2), (s1,s2) in product(product(range(24), repeat=2), product([0, 16, 17], repeat=2)):
	C2Seqs.append( ((c1,c2), "CNOT", (s1,s2)) )

#3. iSWAP like class - replacing iSWAP with (echoCR - (Y90m*Y90m) - echoCR)
for (c1,c2), (s1,s2) in product(product(range(24), repeat=2), product([0, 16, 17], repeat=2)):
	C2Seqs.append( ((c1,c2), "iSWAP", (s1,s2)) )

#4. SWAP like class
for c1, c2 in product(range(24), repeat=2):
	C2Seqs.append( ((c1,c2), "SWAP", None) )

def Cx2(c1,c2, q1, q2):
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
		return ZX90_CR(q2, q1) + [Y90m(q1)*Y90m(q2)] +  ZX90_CR(q2, q1)
	elif gate == "SWAP":
		return ZX90_CR(q2, q1) + [Y90m(q1)*Y90m(q2)] + ZX90_CR(q2, q1) + [(X90(q1)+Y90m(q1))*X90(q2)] + ZX90_CR(q2, q1)

def entangling_mat(gate):
	"""
	Helper function to create the entangling gate matrix
	"""
	echoCR = expm(1j*pi/4*np.kron(pX, pZ))
	if gate == "CNOT":
		return echoCR
	elif gate == "iSWAP":
		return reduce(lambda x,y: np.dot(y,x), [echoCR, np.kron(C1[6], C1[6]), echoCR])
	elif gate == "SWAP":
		return reduce(lambda x,y: np.dot(y,x), [echoCR, np.kron(C1[6], C1[6]), echoCR, np.kron(np.dot(C1[6], C1[1]), C1[1]), echoCR])
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

@memoize
def clifford_mat(c, numQubits):
	"""
	Return the matrix unitary the implements the qubit clifford C
	"""
	assert numQubits <=2, "Oops! I only handle one or two qubits"
	if numQubits == 1:
		return C1[c]
	else:
		c = C2Seqs[c]
		mat = np.kron(clifford_mat(c[0][0], 1), clifford_mat(c[0][1], 1))
		if c[1]:
			mat = np.dot(entangling_mat(c[1]), mat)
		if c[2]:
			mat = np.dot(np.kron(clifford_mat(c[2][0], 1), clifford_mat(c[2][1], 1)), mat)
		return mat

def inverse_clifford(cMat):
	dim = cMat.shape[0]
	if dim == 2:
		for ct in range(24):
			if np.isclose(np.abs(np.dot(cMat, clifford_mat(ct, 1)).trace()), dim):
				return ct
	elif dim == 4:
		for ct in range(len(C2Seqs)):
			if np.isclose(np.abs(np.dot(cMat, clifford_mat(ct, 2)).trace()), dim):
				return ct
	else:
		raise Exception("Expected 2 or 4 qubit dimensional matrix.")


	#If we got here something is wrong
	raise Exception("Couldn't find inverse clifford")








