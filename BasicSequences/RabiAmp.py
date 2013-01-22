'''
Variable amplitude rabi nutation experiment.
'''
from QGL import *

def RabiAmp(qubit, amps, makePlot=False, phase=0):
	seqs = [[Utheta(qubit, amp=amp, phase=phase), MEAS(qubit)] for amp in amps]
	
	fileNames = compile_to_hardware(seqs, 'Rabi')

	if makePlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin


