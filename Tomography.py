'''
Helper functions for adding tomography routines.
'''
from itertools import product
import operator
from PulsePrimitives import *

def create_tomo_blocks(qubits, numPulses, alignment='parallel'):
	'''
	Helper function to create the tomography pulse block in either parallel or serial.
	'''
	#Tomography pulse sets
	if numPulses == 4:
		tomoSet = [Id, X90, Y90, X]
	elif numPulses == 6:
		tomoSet = [Id, X90, X90m, Y90, Y90m, X]
	else:
		raise ValueError("Only able to handle numPulses=4 or 6")

	#Create all combinations of pulses for the number of qubits
	return [reduce(operator.mul, [p(q) for p,q in zip(pulseSet, qubits)]) for pulseSet in product(tomoSet, repeat=len(qubits))]

def state_tomo(seq, qubits=None, numPulses=4):
	'''
	Apply state tomography readout pulses and measurement.
	Expects a single entry list sequence
	'''
	return [seq + [tomoBlock,  MEAS(*qubits)]
				 for tomoBlock in create_tomo_blocks(qubits, numPulses)]

def process_tomo(seq, qubits=None, numPulses=4):
	'''
	Apply process tomography state prep. readout pulses and measurement.
	Expects a single entry list sequence
	'''
	return [[prepBlock] + seq + [readoutBlock,  MEAS(*qubits)]
				 for prepBlock, readoutBlock in product(create_tomo_blocks(qubits, numPulses), repeat=len(qubits)) ]

