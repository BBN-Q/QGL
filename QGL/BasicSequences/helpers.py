from itertools import product
import operator
from ..PulsePrimitives import Id, X, MEAS
from ..ControlFlow import qwait
from functools import reduce


def create_cal_seqs(qubits, numRepeats, measChans=None, waitcmp=False):
    """
	Helper function to create a set of calibration sequences.

	Parameters
	----------
	qubits : logical channels, e.g. (q1,) or (q1,q2) (tuple)
	numRepeats = number of times to repeat calibration sequences (int)
	waitcmp = True if the sequence contains branching
	"""
    if measChans is None:
        measChans = qubits

    calSet = [Id, X]
    #Make all combination for qubit calibration states for n qubits and repeat
    calSeqs = [reduce(operator.mul, [p(q) for p, q in zip(pulseSet, qubits)])
               for pulseSet in product(calSet, repeat=len(qubits))
               for _ in range(numRepeats)]

    #Add on the measurement operator.
    measBlock = reduce(operator.mul, [MEAS(q) for q in qubits])
    return [[seq, measBlock, qwait('CMP')] if waitcmp else [seq, measBlock]
            for seq in calSeqs]

def time_descriptor(times, desired_units="us"):
    if desired_units == "s":
        scale = 1
    elif desired_units == "ms":
        scale = 1e3
    elif desired_units == "us" or desired_units == "μs":
        scale = 1e6
    elif desired_units == "ns":
        scale = 1e9
    axis_descriptor = {
        'name': 'time',
        'unit': desired_units,
        'points': list(scale * times)
    }
    return axis_descriptor
