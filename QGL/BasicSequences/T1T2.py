from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from scipy.constants import pi
from .helpers import create_cal_seqs, delay_descriptor, cal_descriptor


def InversionRecovery(qubit,
                      delays,
                      showPlot=False,
                      calRepeats=2,
                      suffix=False):
    """
    Inversion recovery experiment to measure qubit T1

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    delays : int/float iterable
        Delays after inversion before measurement (seconds). 4 ns minimum.
    showPlot : boolean, optional
        Whether to plot
    calRepeats : int, optional
        How many repetitions of calibration pulses

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = InversionRecovery(q1, np.linspace(20.0e-9, 200.02e-6, 101));
    Compiled 105 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    #Create the basic sequences
    seqs = [[X(qubit), Id(qubit, d), MEAS(qubit)] for d in delays]

    #Tack on the calibration scalings
    seqs += create_cal_seqs((qubit, ), calRepeats)

    metafile = compile_to_hardware(seqs,
        'T1' + ('_' + qubit.label) * suffix + '/T1' +
        ('_' + qubit.label) * suffix,
        axis_descriptor=[
            delay_descriptor(delays),
            cal_descriptor((qubit,), calRepeats)
        ])

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def Ramsey(qubit,
           pulseSpacings,
           TPPIFreq=0,
           showPlot=False,
           calRepeats=2,
           suffix=False):
    """
    Variable pulse spacing Ramsey (pi/2 - tau - pi/2) with optional TPPI.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    pulseSpacings : int/float iterable
        Pulse spacings (seconds). 4 ns minimum.
    TPPIFreq : float, optional
        Frequency for TPPI phase updates of second Ramsey pulse (Hz)
    showPlot : boolean, optional
        Whether to plot
    calRepeats : int, optional
        How many repetitions of calibration pulses

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = Ramsey(q1, np.linspace(20.0e-9, 200.02e-6, 101), TPPIFreq=1.0e6);
    Compiled 105 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    #Create the phases for the TPPI
    phases = 2 * pi * TPPIFreq * pulseSpacings

    #Create the basic Ramsey sequence
    seqs = [[X90(qubit), Id(qubit, d), U90(qubit, phase=phase), MEAS(qubit)]
            for d, phase in zip(pulseSpacings, phases)]

    #Tack on the calibration scalings
    seqs += create_cal_seqs((qubit, ), calRepeats)

    metafile = compile_to_hardware(seqs,
        'Ramsey' + ('_' + qubit.label) * suffix + '/Ramsey' +
        ('_' + qubit.label) * suffix,
        axis_descriptor=[
            delay_descriptor(pulseSpacings),
            cal_descriptor((qubit,), calRepeats)
        ])

    if showPlot:
        plot_pulse_files(metafile)
    return metafile
