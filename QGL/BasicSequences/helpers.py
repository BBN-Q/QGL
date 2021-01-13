# coding=utf-8

from functools import reduce
from itertools import product
import numpy as np
import operator
from typing import Iterable, Union, List, Mapping, Any, Dict

from ..PulsePrimitives import Id, X, MEAS
from ..ControlFlow import qwait
import QGL.Channels as Channels
import QGL.PulsePrimitives as PulseSequencer

def create_cal_seqs(qubits: Channels.LogicalChannel, 
                    numRepeats: int, 
                    measChans: Iterable[Channels.LogicalChannel] = None, 
                    waitcmp: bool = False, 
                    delay: Union[int, float] = None) -> List[List[PulseSequencer.Pulse]]:
    """
    Helper function to create a set of calibration sequences.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
  	numRepeats : int
        Number of times to repeat calibration sequences
  	waitcmp : boolean, optional
        True if the sequence contains branching (sequencing based on values in
        the register). Default = False.
    delay : int/float, optional
        Time between state preparation and measurement (seconds)

    Returns
    -------
    seq : QGL sequence interable
        A list of QGL sequence containing calibration pulses and measurements

    Examples
    --------
    >>> create_cal_seqs((q2,q3), numRepeats=2)
    [[Id(q2)⊗ Id(q3),
      MEAS(M-q2, shape_fun=<autodyne>)⊗ MEAS(M-q3, shape_fun=<autodyne>)],
     [Id(q2)⊗ Id(q3),
      MEAS(M-q2, shape_fun=<autodyne>)⊗ MEAS(M-q3, shape_fun=<autodyne>)],
     [Id(q2)⊗ X(q3),
      MEAS(M-q2, shape_fun=<autodyne>)⊗ MEAS(M-q3, shape_fun=<autodyne>)],
     [Id(q2)⊗ X(q3),
      MEAS(M-q2, shape_fun=<autodyne>)⊗ MEAS(M-q3, shape_fun=<autodyne>)],
     [X(q2)⊗ Id(q3),
      MEAS(M-q2, shape_fun=<autodyne>)⊗ MEAS(M-q3, shape_fun=<autodyne>)],
     [X(q2)⊗ Id(q3),
      MEAS(M-q2, shape_fun=<autodyne>)⊗ MEAS(M-q3, shape_fun=<autodyne>)],
     [X(q2)⊗ X(q3),
      MEAS(M-q2, shape_fun=<autodyne>)⊗ MEAS(M-q3, shape_fun=<autodyne>)],
     [X(q2)⊗ X(q3),
      MEAS(M-q2, shape_fun=<autodyne>)⊗ MEAS(M-q3, shape_fun=<autodyne>)]]
	"""
    if measChans is None:
        measChans = qubits

    calSet = [Id, X]
    #Make all combination for qubit calibration states for n qubits and repeat
    cal_seqs = [reduce(operator.mul, [p(q) for p, q in zip(pulseSet, qubits)])
               for pulseSet in product(calSet, repeat=len(qubits))
               for _ in range(numRepeats)]

    #Add on the measurement operator.
    measBlock = reduce(operator.mul, [MEAS(q) for q in measChans])
    #Add optional delay
    full_cal_seqs = [[seq, Id(qubits[0], delay), measBlock] if delay else [seq, measBlock] for seq in cal_seqs]
    if waitcmp:
        [cal_seq.append(qwait(kind='CMP')) for cal_seq in full_cal_seqs]
    return full_cal_seqs

def cal_descriptor(qubits: Iterable[Channels.LogicalChannel], 
                   numRepeats: int, 
                   partition: int = 2, 
                   states = ['0', '1']) -> Dict[str, Any]:
    # generate state set in same order as we do above in create_cal_seqs()
    state_set = [reduce(operator.add, s) for s in product(states, repeat=len(qubits))]
    descriptor = {
        'name': 'calibration',
        'unit': 'state',
        'partition': partition,
        'points': []
    }
    for state in state_set:
        descriptor['points'] += [state] * numRepeats
    return descriptor

def delay_descriptor(delays: np.ndarray, 
                     desired_units: str = "us") -> Dict[str, Any]:
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
        # Make sure delays is a numpy array so can multiply it by a float safely
        'points': list(scale * np.array(delays)),
        'partition': 1
    }
    return axis_descriptor
