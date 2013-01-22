'''
Variable pulse spacing Ramsey with optional TPPI.
'''
from QGL import *
from scipy.constants import pi

def Ramsey(qubit, pulseSpacings, TPPIFreq=0, showPlot=False):
	phases = 2*pi*TPPIFreq*pulseSpacings
	seqs = [[X90(qubit), Id(qubit, d), U90(qubit, phase=phase), MEAS(qubit)] 
				for d,phase in zip(pulseSpacings, phases)]

	fileNames = compile_to_hardware(seqs, 'Ramsey')

	if showPlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin

