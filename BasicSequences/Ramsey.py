from QGL import *
from scipy.constants import pi

def Ramsey(qubit, pulseSpacings, TPPIFreq=0, showPlot=False):
	"""

	Variable pulse spacing Ramsey (pi/2 - tau - pi/2) with optional TPPI.
	
	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel) 
	pulseSpacings : pulse spacings (iterable; seconds)
	TPPIFreq : frequency for TPPI phase updates of second Ramsey pulse (Hz)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""	

	#Create the phases for the TPPI
	phases = 2*pi*TPPIFreq*pulseSpacings

	#Create the basic Ramsey sequence
	seqs = [[X90(qubit), Id(qubit, d), U90(qubit, phase=phase), MEAS(qubit)] 
				for d,phase in zip(pulseSpacings, phases)]

	#Tack on the calibration scalings
	seqs += [[Id(qubit), MEAS(qubit)], [Id(qubit), MEAS(qubit)], [X(qubit), MEAS(qubit)], [X(qubit), MEAS(qubit)]]
	print('Number of sequences: {0}'.format(len(seqs)))

	fileNames = compile_to_hardware(seqs, 'Ramsey/Ramsey')
	print(fileNames)

	if showPlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin

