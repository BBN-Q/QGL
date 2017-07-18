from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
import QGL.PulseShapes
from .helpers import create_cal_seqs, delay_descriptor, cal_descriptor
from functools import reduce


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

    axis_descriptor = [{
        'name': 'amplitude',
        'unit': None,
        'points': list(amps),
        'partition': 1
    }]

    metafile = compile_to_hardware(seqs, 'Rabi/Rabi', axis_descriptor=axis_descriptor)

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def RabiWidth(qubit,
              widths,
              amp=1,
              phase=0,
              shape_fun=QGL.PulseShapes.tanh,
              showPlot=False):
    """

	Variable pulse width Rabi nutation experiment.

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel)
	widths : pulse widths to sweep over (iterable)
	phase : phase of the pulse (radians, default = 0)
	shape_fun : shape of pulse (function, default = PulseShapes.tanh)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
    seqs = [[Utheta(qubit,
                    length=l,
                    amp=amp,
                    phase=phase,
                    shape_fun=shape_fun), MEAS(qubit)] for l in widths]

    metafile = compile_to_hardware(seqs, 'Rabi/Rabi',
        axis_descriptor=[delay_descriptor(widths)])

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def RabiAmp_NQubits(qubits,
                    amps,
                    phase=0,
                    showPlot=False,
                    measChans=None,
                    docals=False,
                    calRepeats=2):
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

    measBlock = reduce(operator.mul, [MEAS(q) for q in qubits])
    seqs = [[reduce(operator.mul,
                    [Utheta(q, amp=amp, phase=phase) for q in qubits]),
             measBlock] for amp in amps]

    if docals:
        seqs += create_cal_seqs(qubits, calRepeats, measChans=measChans)

    axis_descriptor = [
        {
            'name': 'amplitude',
            'unit': None,
            'points': list(amps),
            'partition': 1
        },
        cal_descriptor(qubits, calRepeats)
    ]

    metafile = compile_to_hardware(seqs, 'Rabi/Rabi',
        axis_descriptor=axis_descriptor)

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


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
    seqs = [[X(mqubit), Utheta(qubit, amp=amp, phase=phase), X(mqubit),
             MEAS(mqubit)] for amp in amps]

    axis_descriptor = [{
        'name': 'amplitude',
        'unit': None,
        'points': list(amps),
        'partition': 1
    }]

    metafile = compile_to_hardware(seqs, 'Rabi/Rabi', axis_descriptor=axis_descriptor)

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def SingleShot(qubit, showPlot=False):
    """
	2-segment sequence with qubit prepared in |0> and |1>, useful for single-shot fidelity measurements and kernel calibration
    """
    seqs = [[Id(qubit), MEAS(qubit)], [X(qubit), MEAS(qubit)]]

    axis_descriptor = {
        'name': 'state',
        'unit': 'state',
        'points': ["0", "1"],
        'partition': 1
    }

    metafile = compile_to_hardware(seqs, 'SingleShot/SingleShot')

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def PulsedSpec(qubit, specOn=True, showPlot=False):
    """
	Measurement preceded by a qubit pulse if specOn = True
    """
    qPulse = X(qubit) if specOn else Id(qubit)
    seqs = [[qPulse, MEAS(qubit)]]
    metafile = compile_to_hardware(seqs, 'Spec/Spec')

    if showPlot:
        plot_pulse_files(metafile)
    return metafile
