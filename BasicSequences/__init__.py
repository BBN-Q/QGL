from RabiAmp import RabiAmp
from Ramsey import Ramsey
from FlipFlop import FlipFlop
from SPAM import SPAM
from RB import SingleQubitRB, SingleQubitRB_AC, SingleQubitRBT




from itertools import product
import operator
from ..PulsePrimitives import Id, X

def create_cal_seqs(qubits, numCals):
	"""
	Helper function to create a set of calibration sequences.
	"""
	calSet = [Id, X]
	calSeqs = [reduce(operator.mul, [p(q) for p,q in zip(pulseSet, qubits)]) for pulseSet in product(calSet, repeat=len(qubits))]
	return reduce(operator.add, [[[seq]]*numCals for seq in calSeqs])
