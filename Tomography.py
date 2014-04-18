'''
Helper functions for adding tomography routines.

Copyright 2013 Raytheon BBN Technologies

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

def state_tomo(seq, qubits, numPulses=4, measChans=None):
	'''
	Apply state tomography readout pulses and measurement.

	Parameters
	-----------
	seq : a single entry list sequence to perform tomography on
	qubits : which qubits to act on 
	numPulses : number of readout pulses
	measChans : tuple of measurement channels to readout (defaults to individual qubit channels)
	'''
	if measChans is None:
		measChans = qubits
	
	return [seq + [tomoBlock,  MEAS(*measChans)]
				 for tomoBlock in create_tomo_blocks(qubits, numPulses)]

def process_tomo(seq, qubits, numPulses=4, measChans=None):
	'''
	Apply process tomography state prep. and readout pulses and measurement.

	Parameters
	-----------
	seq : a single entry list sequence to perform tomography on
	qubits : which qubits to act on 
	numPulses : number of prep/readout pulses
	measChans : tuple of measurement channels to readout (defaults to individual qubit channels)
	'''
	if measChans is None:
		measChans = qubits
	
	return [[prepBlock] + seq + [readoutBlock,  MEAS(*measChans)]
				 for prepBlock, readoutBlock in product(create_tomo_blocks(qubits, numPulses), repeat=2)]

