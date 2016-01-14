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


def RabiAmp_NQubits(qubits, amps, phase=0, showPlot=False, measChans=None, docals=False, calRepeats=2):
	"""

	Variable amplitude Rabi nutation experiment for an arbitrary number of qubits simultaneously

	Parameters
	----------
	qubits : tuple of logical channels to implement sequence (LogicalChannel)
	amps : pulse amplitudes to sweep over for all qubits (iterable)
	phase : phase of the pulses (radians)
	showPlot : whether to plot (boolean)
    measChans : tuble of qubits to be measured (LogicalChannel)
    docals, calRepeats: enable calibration sequences, repeated calRepeats times

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
	if measChans is None:
		measChans = qubits

	seqs = [[reduce(operator.mul, [Utheta(q, amp=amp, phase=phase) for q in qubits]),MEAS(*measChans)] for amp in amps]

	if docals:
		seqs += create_cal_seqs(qubits, calRepeats, measChans=measChans)

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
