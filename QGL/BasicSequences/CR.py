from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..ChannelLibraries import EdgeFactory
from ..PulseSequencePlotter import plot_pulse_files
from .helpers import create_cal_seqs, delay_descriptor, cal_descriptor
import numpy as np
from itertools import product

def PiRabi(controlQ,
           targetQ,
           lengths,
           riseFall=40e-9,
           amp=1,
           phase=0,
           calRepeats=2,
           showPlot=False):
    """
	Variable length CX experiment.

    Parameters
    ----------
    controlQ : Channels.LogicalChannel
        Logical channel for the control qubit
    targetQ: Channels.LogicalChannel
        Logical channel for the target qubit
    lengths : int/float iterable
        Pulse lengths of the CR pulse to sweep over (seconds).  4 ns minimum.
    riseFall : float, optional
        Rise/fall time of the CR pulse (seconds)
    amp : float, optional
        Amplitude of the CR pulse. Valid range: [0.0, 1.0].
    phase : float, optional
        Phase of the CR pulse (radians)
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files.

    Examples
    --------
    >>> mf = PiRabi(q1, q2, np.linspace(20.0e-9, 200.02e-6, 101));
    Compiled 210 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""

    CRchan = EdgeFactory(controlQ, targetQ)
    seqs = [[Id(controlQ),
             flat_top_gaussian(CRchan, riseFall, amp=amp, phase=phase, length=l),
             MEAS(targetQ)*MEAS(controlQ)] for l in lengths] + \
           [[X(controlQ),
             flat_top_gaussian(CRchan, riseFall, amp=amp, phase=phase, length=l),
             X(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for l in lengths] + \
              create_cal_seqs([targetQ,controlQ], calRepeats, measChans=(targetQ,controlQ))

    metafile = compile_to_hardware(seqs, 'PiRabi/PiRabi',
        axis_descriptor=[
            delay_descriptor(np.concatenate((lengths, lengths))),
            cal_descriptor((controlQ, targetQ), calRepeats)
        ])

    if showPlot:
        plot_pulse_files(metafile)

    return metafile


def EchoCRLen(controlQ,
              targetQ,
              lengths,
              riseFall=40e-9,
              amp=1,
              phase=0,
              calRepeats=2,
              showPlot=False, canc_amp=0, canc_phase=np.pi/2):
    """
	Variable length CX experiment, with echo pulse sandwiched between two CR
    opposite-phase pulses.  This is primarily used as a subroutine
    in calibration.

	Parameters
	----------
	controlQ : Channels.LogicalChannel
        Logical channel for the control qubit
	targetQ : Channels.LogicalChannel
        Logical channel for the target qubit
	lengths : int/float iterable
        Pulse lengths of the CR pulse to sweep over (seconds)
	riseFall : float, optional
        Rise/fall time of the CR pulse (seconds)
	amp : float, optional
        Amplitude of the CR pulse. Valid range: [0.0, 1.0].
	phase : float, optional
        Phase of the CR pulse (radian)
	calRepeats : int, optional
        Number of calibrations repeats for each 2-qubit state basis state
	showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = EchoCRLen(q1, q2, np.linspace(20.0e-9, 200.02e-6, 101));
    Compiled 210 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""
    seqs = [[Id(controlQ),
             echoCR(controlQ, targetQ, length=l, phase=phase, amp=amp, riseFall=riseFall, canc_amp=canc_amp, canc_phase=canc_phase),
             Id(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for l in lengths] + \
           [[X(controlQ),
             echoCR(controlQ, targetQ, length=l, phase= phase, amp=amp, riseFall=riseFall, canc_amp=canc_amp, canc_phase=canc_phase),
             X(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for l in lengths] + \
           create_cal_seqs((controlQ,targetQ), calRepeats, measChans=(targetQ,controlQ))

    metafile = compile_to_hardware(seqs, 'EchoCR/EchoCR',
        axis_descriptor=[
            delay_descriptor(np.concatenate((lengths, lengths))),
            cal_descriptor((controlQ, targetQ), calRepeats)
        ])

    if showPlot:
        plot_pulse_files(metafile)

    return metafile


def EchoCRPhase(controlQ,
                targetQ,
                phases,
                riseFall=40e-9,
                amp=1,
                length=100e-9,
                calRepeats=2,
                showPlot=False,
                canc_amp=0,
                canc_phase=np.pi/2):
    """
    Variable phase CX experiment, with echo pulse sandwiched between two CR
    opposite-phase pulses.  This is primarily used as a subroutine
    in calibration.

	Parameters
	----------
	controlQ : Channels.LogicalChannel
        Logical channel for the control qubit
	targetQ : Channels.LogicalChannel
        Logical channel for the target qubit
	phases : float iterable
        Pulse phases of the CR pulse to sweep over (radians)
	riseFall : float, optional
        Rise/fall time of the CR pulse (seconds)
	amp : float, optional
        Amplitude of the CR pulse. Valid range: [0.0, 1.0].
	length : float, optional
        Duration of each of the two flat parts of the CR pulse (seconds)
	calRepeats : int, optional
        Number of calibrations repeats for each 2-qubit state basis state
	showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = EchoCRPhase(q1, q2, np.linspace(0.0, np.pi, 51));
    Compiled 110 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""
    seqs = [[Id(controlQ),
             echoCR(controlQ, targetQ, length=length, phase=ph, amp=amp, riseFall=riseFall, canc_amp=canc_amp, canc_phase=canc_phase),
             X90(targetQ)*Id(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for ph in phases] + \
           [[X(controlQ),
             echoCR(controlQ, targetQ, length=length, phase= ph, amp=amp, riseFall = riseFall, canc_amp=canc_amp, canc_phase=canc_phase),
             X90(targetQ)*X(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for ph in phases] + \
             create_cal_seqs((controlQ, targetQ), calRepeats, measChans=(targetQ,controlQ))

    axis_descriptor = [
        {
            'name': 'phase',
            'unit': 'radians',
            'points': list(phases)+list(phases),
            'partition': 1
        },
        cal_descriptor((controlQ, targetQ), calRepeats)
    ]

    metafile = compile_to_hardware(seqs, 'EchoCR/EchoCR',
        axis_descriptor=axis_descriptor)

    if showPlot:
        plot_pulse_files(metafile)

    return metafile


def EchoCRAmp(controlQ,
              targetQ,
              amps,
              riseFall=40e-9,
              length=50e-9,
              phase=0,
              calRepeats=2,
              showPlot=False):
    """
    Variable amplitude CX experiment, with echo pulse sandwiched between two
    CR opposite-phase pulses.

	Parameters
	----------
    controlQ : Channels.LogicalChannel
        Logical channel for the control qubit
    targetQ : Channels.LogicalChannel
        Logical channel for the target qubit
    amps : float iterable
        Pulse amplitudes of the CR pulse to sweep over. Valid range: [0.0, 1.0]
    riseFall : float, optional
        Rise/fall time of the CR pulse (seconds)
    length : float, optional
        Duration of each of the two flat parts of the CR pulse (seconds)
    phase : float, optional
        Phase of the CR pulse (radians)
    calRepeats : int, optional
        Number of calibrations repeats for each 2-qubit state basis state
    showPlot : whether to plot (boolean)

    Returns
    -------
    metafile : path to a json metafile with details about the sequences and paths to compiled machine files

    Examples
    --------
    >>> mf = EchoCRAmp(q1, q2, np.linspace(0.7, 0.9, 101));
    Compiled 105 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""
    seqs = [[Id(controlQ),
             echoCR(controlQ, targetQ, length=length, phase=phase, riseFall=riseFall,amp=a),
             Id(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for a in amps] + \
           [[X(controlQ),
             echoCR(controlQ, targetQ, length=length, phase= phase, riseFall=riseFall,amp=a),
             X(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for a in amps] + \
           create_cal_seqs((controlQ ,targetQ), calRepeats, measChans=(targetQ,controlQ))

    axis_descriptor = [
        {
            'name': 'amplitude',
            'unit': None,
            'points': list(amps)+list(amps),
            'partition': 1
        },
        cal_descriptor((controlQ, targetQ), calRepeats)
    ]

    metafile = compile_to_hardware(seqs, 'EchoCR/EchoCR',
        axis_descriptor=axis_descriptor)
    if showPlot:
        plot_pulse_files(metafile)

    return metafile

def CRtomo_seq(controlQ,
               targetQ,
               lengths,
               phase,
               amp=0.8,
               riseFall=20e-9,
               calRepeats=2):
    """
    Variable length CX experiment, for Hamiltonian tomography.

    Parameters
    ----------
    controlQ : Channels.LogicalChannel
        Logical channel for the control qubit
    targetQ : Channels.LogicalChannel
        Logical channel for the target qubit
    lengths : int/float iterable
        Pulse lengths of the CR pulse to sweep over (seconds)
    phase : float
        Phase of the CR pulse (radians)
    amps : float, optional
        Pulse amplitude of the CR pulse. Valid range: [0.0, 1.0]
    riseFall : float, optional
        Rise/fall time of the CR pulse (seconds)
    calRepeats : int, optional
        Number of calibrations repeats for each 2-qubit state basis state

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = CRtomo_seq(q2, q3, np.linspace(20.0e-9, 200.02e-6, 101), \
            phase=0.0);
    Compiled 610 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """
    CRchan = ChannelLibraries.EdgeFactory(controlQ, targetQ)
    tomo_pulses = [Y90m, X90, Id]
    seqs = [[Id(controlQ),
         flat_top_gaussian(CRchan, amp=amp, riseFall=riseFall, length=l, phase=phase, label="CR"),
         Id(controlQ)*tomo_pulse(targetQ),
         MEAS(targetQ)] for l,tomo_pulse in product(lengths, tomo_pulses)] + \
       [[X(controlQ),
         flat_top_gaussian(CRchan, amp=amp, riseFall=riseFall, length=l, phase=phase, label="CR"),
         X(controlQ)*tomo_pulse(targetQ),
         MEAS(targetQ)] for l,tomo_pulse in product(lengths, tomo_pulses)] + \
       create_cal_seqs((targetQ,), 2,)
    metafile = compile_to_hardware(seqs, 'CR/CR',
        axis_descriptor=[
            delay_descriptor(np.concatenate((np.repeat(lengths,3), np.repeat(lengths,3)))),
            cal_descriptor((targetQ,), 2)
        ])
    return metafile
