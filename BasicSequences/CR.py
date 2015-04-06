from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
import QGL.PulseShapes
from helpers import create_cal_seqs

def PiRabi(controlQ, targetQ, CRchan, lengths, riseFall=40e-9, amp=1, phase=0, calRepeats=2, showPlot=False):
	"""
	Variable length CX experiment.

	Parameters
	----------
	controlQ : logical channel for the control qubit (LogicalChannel)
	CRchan: logical channel for the cross-resonance pulse (LogicalChannel) 
	lengths : pulse lengths of the CR pulse to sweep over (iterable)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""

	seqs = [[Id(controlQ)] + flat_top_gaussian(CRchan, riseFall, amp=amp, phase=phase, length=l) \
	+ [MEAS(targetQ)*MEAS(controlQ)] for l in lengths]+[[X(controlQ)] + flat_top_gaussian(CRchan, riseFall, amp=amp, phase=phase, length=l)\
	+ [X(controlQ), MEAS(targetQ)*MEAS(controlQ)] for l in lengths] + create_cal_seqs([targetQ,controlQ], calRepeats, measChans=(targetQ,controlQ))

	fileNames = compile_to_hardware(seqs, 'PiRabi/PiRabi')
	print(fileNames)

	if showPlot:
		plot_pulse_files(fileNames)


def EchoCRLen(controlQ, targetQ, CRchan, lengths, riseFall=40e-9, amp=1, phase=0, calRepeats=2, showPlot=False):
	"""
	Variable length CX experiment, with echo pulse sandwiched between two CR opposite-phase pulses.

	Parameters
	----------
	controlQ : logical channel for the control qubit (LogicalChannel)
	CRchan: logical channel for the cross-resonance pulse (LogicalChannel) 
	lengths : pulse lengths of the CR pulse to sweep over (iterable)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
	seqs = [[Id(controlQ)] + echoCR(controlQ, CRchan, length=l, phase=phase, riseFall=riseFall) + [Id(controlQ), MEAS(targetQ)*MEAS(controlQ)]\
	 for l in lengths]+ [[X(controlQ)] + echoCR(controlQ, CRchan, length=l, phase= phase, riseFall=riseFall) + [X(controlQ), MEAS(targetQ)*MEAS(controlQ)]\
	  for l in lengths] + create_cal_seqs((targetQ,controlQ), calRepeats, measChans=(targetQ,controlQ))

	fileNames = compile_to_hardware(seqs, 'EchoCR/EchoCR')
	print(fileNames)

	if showPlot:
		plot_pulse_files(fileNames)


def EchoCRPhase(controlQ, targetQ, CRchan, phases, riseFall=40e-9, amp=1, length=100e-9, calRepeats=2, showPlot=False):
	"""
	Variable phase CX experiment, with echo pulse sandwiched between two CR opposite-phase pulses.

	Parameters
	----------
	controlQ : logical channel for the control qubit (LogicalChannel)
	CRchan: logical channel for the cross-resonance pulse (LogicalChannel) 
	phases : pulse phases of the CR pulse to sweep over (iterable)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
	seqs = [[Id(controlQ)] + echoCR(controlQ, CRchan, length=length, phase=ph, riseFall=riseFall) + [X90(targetQ)*Id(controlQ), MEAS(targetQ)*MEAS(controlQ)] \
	for ph in phases]+[[X(controlQ)] + echoCR(controlQ, CRchan, length=length, phase= ph, riseFall = riseFall) + [X90(targetQ)*X(controlQ), MEAS(targetQ)*MEAS(controlQ)]\
	 for ph in phases]+create_cal_seqs((targetQ,controlQ), calRepeats, measChans=(targetQ,controlQ))

	fileNames = compile_to_hardware(seqs, 'EchoCR/EchoCR')
	print(fileNames)

	if showPlot:
		plot_pulse_files(fileNames)




