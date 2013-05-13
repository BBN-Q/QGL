from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files

def HahnEcho(qubit, pulseSpacings, showPlot=False):
	"""
	A single pulse Hahn echo. 

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel) 
	pulseSpacings : pulse spacings to sweep over; the t in 90-t-180-t-180 (iterable)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
	seqs = [ [X90(qubit), Id(qubit, t), X(qubit), Id(qubit,t), X90(qubit) ,MEAS(qubit)] 
					 for t in pulseSpacings]

	fileNames = compile_to_hardware(seqs, 'Echo/Echo')
	print(fileNames)

	if showPlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin

