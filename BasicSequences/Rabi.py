from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
import QGL.PulseShapes
from helpers import create_cal_seqs


def RabiAmp(qubit, amps, phase=0, showPlot=False):
	"""

	Variable amplitude Rabi nutation experiment.

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel)
	amps : pulse amplitudes to sweep over (iterable)
	phase : phase of the pulse (radians)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
	seqs = [[Utheta(qubit, amp=amp, phase=phase), MEAS(qubit)] for amp in amps]

	fileNames = compile_to_hardware(seqs, 'Rabi/Rabi')
	print(fileNames)

	if showPlot:
		plot_pulse_files(fileNames)

def RabiWidth(qubit, widths, amp=1, phase=0, shapeFun=QGL.PulseShapes.tanh, showPlot=False):
	"""

	Variable pulse width Rabi nutation experiment.

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel)
	widths : pulse widths to sweep over (iterable)
	phase : phase of the pulse (radians, default = 0)
	shapeFun : shape of pulse (function, default = PulseShapes.tanh)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
	seqs = [[Utheta(qubit, length=l, amp=amp, phase=phase, shapeFun=shapeFun), MEAS(qubit)] for l in widths]

	fileNames = compile_to_hardware(seqs, 'Rabi/Rabi')
	print(fileNames)

	if showPlot:
		plot_pulse_files(fileNames)


def RabiAmp_TwoQubits(qubit1, qubit2, amps, amps2, phase=0, showPlot=False, meas=[1,1],docals=False):
	"""

	Variable amplitude Rabi nutation experiment for up to two qubits, with measurement on both. Need to be extended
	to arbitrary number of qubits

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel)
	amps : pulse amplitudes to sweep over for qubit 1(iterable)
	amps2: pulse amplitudes to sweep over for qubit 2(iterable, same index)
	phase : phase of the pulses (radians)
	showPlot : whether to plot (boolean)
	meas: list of 1/0 for measurement on/off

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
	meas1 = MEAS(qubit1) if meas[0]==1 else Id(qubit1)
  	meas2 = MEAS(qubit2) if meas[1]==1 else Id(qubit2)
	seqs = [[Utheta(qubit1, amp=amp, phase=phase)*Utheta(qubit2, amp=amp2,phase=phase), meas1*meas2] for (amp,amp2) in zip(amps,amps2)]

	if docals:
		seqs += create_cal_seqs((qubit1,qubit2), 2, measChans=(qubit1,qubit2))

	fileNames = compile_to_hardware(seqs, 'Rabi/Rabi')
	print(fileNames)

	if showPlot:
		plot_pulse_files(fileNames)

def RabiAmpPi(qubit, mqubit, amps, phase=0, showPlot=False):
	"""

	Variable amplitude Rabi nutation experiment.

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel)
	amps : pulse amplitudes to sweep over (iterable)
	phase : phase of the pulse (radians)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
	seqs = [[X(mqubit),Utheta(qubit, amp=amp, phase=phase), X(mqubit), MEAS(mqubit)] for amp in amps]

	fileNames = compile_to_hardware(seqs, 'Rabi/Rabi')
	print(fileNames)

	if showPlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin

def SingleShot(qubit, showPlot = False):
	"""
	2-segment sequence with qubit prepared in |0> and |1>, useful for single-shot fidelity measurements and kernel calibration
    """
	seqs = [[Id(qubit), MEAS(qubit)], [X(qubit), MEAS(qubit)]]
	filenames = compile_to_hardware(seqs, 'SingleShot/SingleShot')
	print(filenames)

	if showPlot:
		plot_pulse_files(filenames)

def PulsedSpec(qubit, specOn = True, showPlot = False):
	"""
	Measurement preceded by a qubit pulse if specOn = True
    """
	qPulse = X(qubit) if specOn else Id(qubit)
	seqs = [[qPulse, MEAS(qubit)]]
	filenames = compile_to_hardware(seqs, 'Spec/Spec')
	print(filenames)

	if showPlot:
		plot_pulse_files(filenames)

def Swap(qubit, mqubit, delays, showPlot=False):
	"""

	Variable amplitude Rabi nutation experiment.

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel)
	amps : pulse amplitudes to sweep over (iterable)
	phase : phase of the pulse (radians)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
	seqs = [[X(qubit), X(mqubit), Id(mqubit, d), MEAS(mqubit)*MEAS(qubit)] for d in delays] + create_cal_seqs((mqubit,qubit), 2, measChans=(mqubit,qubit))

	fileNames = compile_to_hardware(seqs, 'Rabi/Rabi')
	print(fileNames)

	if showPlot:
		plotWin = plot_pulse_files(fileNames)
		return plotWin
