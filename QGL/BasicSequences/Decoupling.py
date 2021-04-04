from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from .helpers import create_cal_seqs, delay_descriptor, cal_descriptor

from typing import Iterable, Union

def HahnEcho(qubit: Channels.LogicalChannel,
             pulseSpacings: Iterable[Union[int, float]],
             periods: int = 0, 
             calRepeats: int = 2, 
             showPlot: bool = False) -> str:
    """
	A single pulse Hahn echo with variable phase of second pi/2 pulse.

	Parameters
	----------
	qubit :  Channels.LogicalChannel
        Logical channel to implement sequence
	pulseSpacings : int iterable
        Pulse spacings to sweep over; the t in 90-t-180-t-180 (iterable)
	periods : int, optional
        Number of artificial oscillations
	calRepeats : int, optional
        How many times to repeat calibration scalings (default 2)
	showPlot : boolean, optional
        Whether to plot

	Returns
	-------
	metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = HahnEcho(q1, np.linspace(20.0e-9, 200.02e-6, 101), periods=6);
    Compiled 105 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""
    seqs = []
    for k in range(len(pulseSpacings)):
        seqs.append([X90(qubit), Id(qubit, pulseSpacings[k]), \
        Y(qubit), Id(qubit,pulseSpacings[k]), \
        U90(qubit,phase=2*pi*periods/len(pulseSpacings)*k), MEAS(qubit)])

#Tack on the calibration scalings
    seqs += create_cal_seqs((qubit, ), calRepeats)

    metafile = compile_to_hardware(seqs, 'Echo/Echo',
        axis_descriptor=[
            delay_descriptor(2 * pulseSpacings),
            cal_descriptor((qubit,), calRepeats)
        ])

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def CPMG(qubit: Channels.LogicalChannel,
         numPulses: Iterable[int], 
         pulseSpacing: float, 
         calRepeats: int = 2, 
         showPlot: bool = False):
    """
	CPMG pulse train with fixed pulse spacing. Note this pulse spacing is centre to centre,
	i.e. it accounts for the pulse width

	Parameters
	----------
	qubit : Channels.LogicalChannel
        Logical channel to implement sequence
	numPulses : int iterable
        Number of 180 pulses; should be even
	pulseSpacing : float
        Spacing between the 180's (seconds)
	calRepeats : int, optional
        How many times to repeat calibration scalings (default 2)
	showPlot : bool, optional
        Whether to plot (boolean)

	Returns
	-------
	metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = CPMG(q2, [2**n for n in range(8)], 20.0e-6);
    Compiled 12 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""
    #First setup the t-180-t block
    CPMGBlock = [Id(qubit, (pulseSpacing - qubit.pulse_params['length']) / 2),
                 Y(qubit),
                 Id(qubit, (pulseSpacing - qubit.pulse_params['length']) / 2)]

    seqs = [[X90(qubit)] + CPMGBlock * rep + [X90(qubit), MEAS(qubit)]
            for rep in numPulses]

    #Tack on the calibration scalings
    seqs += create_cal_seqs((qubit, ), calRepeats)

    metafile = compile_to_hardware(seqs, 'CPMG/CPMG',
        axis_descriptor=[
            # NOTE: numPulses is often not a numpy array, so cannot multiply
            # by a float. But thankfully, np.array(np.array) = np.array so
            # this is always a good move here.
            delay_descriptor(pulseSpacing * np.array(numPulses)),
            cal_descriptor((qubit,), calRepeats)
        ])

    if showPlot:
        plot_pulse_files(metafile)
    return metafile

def JAZZ(qubit: Channels.LogicalChannel,
         spec_qubit: Channels.LogicalChannel,
         pulseSpacings: Iterable[Union[int, float]],
         periods: int = 0, 
         calRepeats: int = 2, 
         showPlot: bool = False) -> str:
    """
    Implement a JAZZ sequence on a pair of coupled qubits. A Ramsey experiment
    is run on the first qubit and pi pulses are toggled on the second.

    Parameters
    ----------
    qubit :  Channels.LogicalChannel
        Logical channel to measure
    spec_qubit :  Channels.LogicalChannel
        Logical channel to echo
    pulseSpacings : int iterable
        Pulse spacings to sweep over; the t in 90-t-180-t-180 (iterable)
    periods : int, optional
        Number of artificial oscillations
    calRepeats : int, optional
        How many times to repeat calibration scalings (default 2)
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = JAZZ(q1, q2, np.linspace(20.0e-9, 200.02e-6, 101), periods=6);
    Compiled 105 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """
    seqs = []
    for k in range(len(pulseSpacings)):
        seqs.append([X90(qubit), Id(qubit, pulseSpacings[k]), \
        Y(qubit)*X(spec_qubit), Id(qubit,pulseSpacings[k]), \
        U90(qubit,phase=2*pi*periods/len(pulseSpacings)*k), MEAS(qubit)])

    for k in range(len(pulseSpacings)):
        seqs.append([X90(qubit)*X(spec_qubit), Id(qubit, pulseSpacings[k]), \
        Y(qubit)*X(spec_qubit), Id(qubit,pulseSpacings[k]), \
        U90(qubit,phase=2*pi*periods/len(pulseSpacings)*k), MEAS(qubit)])

    # Tack on the calibration scalings
    seqs += create_cal_seqs((qubit, ), calRepeats)

    metafile = compile_to_hardware(seqs, 'JAZZ/JAZZ',
        axis_descriptor=[
            delay_descriptor(np.hstack(np.squeeze([2 * pulseSpacings] *2))),
            cal_descriptor((qubit,), calRepeats)
        ])

    if showPlot:
        plot_pulse_files(metafile)
    return metafile
