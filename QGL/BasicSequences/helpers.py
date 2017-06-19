# coding=utf-8

from itertools import product
import operator
from ..PulsePrimitives import Id, X, MEAS
from ..ControlFlow import qwait
from functools import reduce

def create_cal_seqs(qubits, numRepeats, measChans=None, waitcmp=False, delay=None):
    """
	Helper function to create a set of calibration sequences.

	Parameters
	----------
	qubits : logical channels, e.g. (q1,) or (q1,q2) (tuple)
	numRepeats = number of times to repeat calibration sequences (int)
	waitcmp = True if the sequence contains branching
    delay: optional time between state preparation and measurement (s)
	"""
    if measChans is None:
        measChans = qubits

    calSet = [Id, X]
    #Make all combination for qubit calibration states for n qubits and repeat
    cal_seqs = [reduce(operator.mul, [p(q) for p, q in zip(pulseSet, qubits)])
               for pulseSet in product(calSet, repeat=len(qubits))
               for _ in range(numRepeats)]

    #Add on the measurement operator.
    measBlock = reduce(operator.mul, [MEAS(q) for q in qubits])
    #Add optional delay
    full_cal_seqs = [[seq, Id(qubits[0], delay), measBlock] if delay else [seq, measBlock] for seq in cal_seqs]
    if waitcmp:
        [cal_seq.append(qwait(kind='CMP')) for cal_seq in full_cal_seqs]
    return full_cal_seqs

def cal_descriptor(qubits, numRepeats):
    states = ['0', '1']
    # generate state set in same order as we do above in create_cal_seqs()
    state_set = [reduce(operator.add, s) for s in product(states, repeat=len(qubits))]
    descriptor = {
        'name': 'calibration',
        'unit': 'state',
        'partition': 2,
        'points': []
    }
    for state in state_set:
        descriptor['points'] += [state] * numRepeats
    return descriptor

def delay_descriptor(delays, desired_units="us"):
    if desired_units == "s":
        scale = 1
    elif desired_units == "ms":
        scale = 1e3
    elif desired_units == "us" or desired_units == u"μs":
        scale = 1e6
    elif desired_units == "ns":
        scale = 1e9
    axis_descriptor = {
        'name': 'delay',
        'unit': desired_units,
        'points': list(scale * delays),
        'partition': 1
    }
    return axis_descriptor
