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
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    amps : int/float iterable
        Array-like iterable of amplitudes to sweep over. [-1, 1]
    phase : float, optional
        Phase of the Rabi pulse (radians). Default = 0
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = RabiAmp(q1, np.linspace(-1.0, 1.0, 101));
    Compiled 101 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""
    seqs = [[Utheta(qubit, amp=amp, phase=phase), MEAS(qubit)] for amp in amps]

    axis_descriptor = [{
        'name': 'amplitude',
        'unit': None,
        'points': list(amps),
        'partition': 1
    }]

    metafile = compile_to_hardware(seqs,
                                   'Rabi/Rabi',
                                   axis_descriptor=axis_descriptor)

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
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    widths : int/float iterable
        Lengths of the Rabi pulse to sweep through (seconds). 4 ns minimum.
    phase : float, optional
        Phase of the Rabi pulse (radians). Default = 0.
    shape_fun : function, optional
        Pulse shape to use for the RabiWidth experiments.  Some useful ones are
        avaliable in QGL.PulseShapes but you can write your own.  See those in
        QGL.PulseShapes for examples.
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = Rabiwidth(q1, np.linspace(20.0e-9, 2020.0e-9, 101));
    Compiled 101 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""
    seqs = [[Utheta(qubit,
                    length=l,
                    amp=amp,
                    phase=phase,
                    shape_fun=shape_fun), MEAS(qubit)] for l in widths]

    metafile = compile_to_hardware(seqs,
                                   'Rabi/Rabi',
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
    Variable amplitude Rabi nutation experiment for an arbitrary number of
    qubits simultaneously

    Parameters
    ----------
    qubits : Channels.LogicalChannel tupple
        A hashable (immutable) tupple of qubits for the Rabi experiment
    amps : int/float iterable
        Array-like iterable of amplitudes to sweep over. [-1, 1]
    phase : float, optional
        Phase of the Rabi pulse (radians). Default = 0
    showPlot : boolean, optional
        Whether to plot
    measChans : Channels.LogicalChannel tupple, optional
        A hashable (immutable) tupple of qubits to measured.
    docals : boolean, optional
        Whether to append calibration pulses to the end of the sequence
    calRepeats : int, optional
        How many times to repeat calibration scalings (default 2)

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = RabiAmp_NQubits((q1,q2), np.linspace(-1.0, 1.0, 101));
    Compiled 101 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""
    if measChans is None:
        measChans = qubits

    measBlock = reduce(operator.mul, [MEAS(q) for q in measChans])
    if shape_fun:
        seqs = [[reduce(operator.mul,
                        [Utheta(q, amp=amp, phase=phase, shape_fun) for q in qubits]),
                 measBlock] for amp in amps]
    else:
    seqs = [[reduce(operator.mul,
                    [Utheta(q, amp=amp, phase=phase) for q in qubits]),
             measBlock] for amp in amps]

    if docals:
        seqs += create_cal_seqs(qubits, calRepeats, measChans=measChans)
    else:
        calRepeats = 0

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
	Variable amplitude Rabi nutation experiment with the state of a second,
    spectator qubit flipped for the duration of the Rabi pulse.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel for Rabi pulse
    mqubit : Channels.LogicalChannel
        Logical channel to invert for the Rabi pulse duration
    amps : int/float iterable
        Array-like iterable of amplitudes to sweep over. [-1, 1]
    phase : float, optional
        Phase of the Rabi pulse (radians). Default = 0
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = RabiAmpPi(q1, q2, np.linspace(-1.0, 1.0, 101));
    Compiled 101 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""
    seqs = [[X(mqubit), Utheta(qubit, amp=amp, phase=phase), X(mqubit),
             MEAS(mqubit)] for amp in amps]

    axis_descriptor = [{
        'name': 'amplitude',
        'unit': None,
        'points': list(amps),
        'partition': 1
    }]

    metafile = compile_to_hardware(seqs,
                                   'Rabi/Rabi',
                                   axis_descriptor=axis_descriptor)

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def SingleShot(qubit, showPlot=False):
    """
    2-segment sequence with qubit prepared in |0> and |1>, useful for
    single-shot fidelity measurements and kernel calibration.  It produces a
    simple sequence of [[Id(qubit), MEAS(qubit)], [X(qubit), MEAS(qubit)]]

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to apply single-shot sequence
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = SingleShot(q1);
    Compiled 1 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
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

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to apply single-shot sequence
    specon : boolean, optional
        Toggles the spectroscopy pulse on and off.  Default = True
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = PulsedSpec(q1);
    Compiled 1 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """
    qPulse = X(qubit) if specOn else Id(qubit)
    seqs = [[qPulse, MEAS(qubit)]]
    metafile = compile_to_hardware(seqs, 'Spec/Spec')

    if showPlot:
        plot_pulse_files(metafile)
    return metafile
