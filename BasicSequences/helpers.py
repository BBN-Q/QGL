from itertools import product
import operator
from ..PulsePrimitives import Id, X, MEAS

def create_cal_seqs(qubits, numRepeats, measChans=None):
	"""
	Helper function to create a set of calibration sequences.

	Parameters
	----------
	qubits : logical channels, e.g. (q1,) or (q1,q2) (tuple) 
	numRepeats = number of times to repeat calibration sequences (int)
	"""
	if measChans is None:
		measChans = qubits

	calSet = [Id, X]
	calSeqs = [reduce(operator.mul, [p(q) for p,q in zip(pulseSet, qubits)]) for pulseSet in product(calSet, repeat=len(qubits))]

	#Add on the measurement operator and repeat
	return reduce(operator.add, [[[seq]+[MEAS(*measChans)]]*numRepeats for seq in calSeqs])
