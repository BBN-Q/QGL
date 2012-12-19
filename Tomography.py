'''
Helper functions for adding tomography routines.
'''

from PulsePrimitives import *

def state_tomo(seq, qubits=None, numPulses=4):

	#Tomography pulses
	if numPulses == 4:
		tomoSet = [Id, X90, Y90, X]
	elif numPulses == 6:
		tomoSet = [Id, X90, X90m, Y90, Y90m, X]
	else:
		raise NameError("Only able to handle numPulses=4 or 6")

	q1, q2 = qubits
	return [seq + [p1(q1)*p2(q2), MEAS(q1,q2)] for p1 in tomoSet for p2 in tomoSet]

