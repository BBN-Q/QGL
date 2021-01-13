from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files

from itertools import chain
from numpy import pi
from typing import Iterable, Union

def SPAM(qubit: Channels.LogicalChannel,
         angleSweep: Iterable[Union[int, float]], 
         maxSpamBlocks: int = 10, 
         showPlot: bool = False) -> str:
    """
    X-Y sequence (X-Y-X-Y)**n to determine quadrature angles or mixer
    correction.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    angleSweep : int/float iterable
        Array-like iterable of angles to sweep over (radians).
    maxSpamBlocks : int, optional
        Number of (X-Y-X-Y) sequences to include
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = SPAM(q1, np.linspace(-1.0, 1.0, 11));
    Compiled 122 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    def spam_seqs(angle):
        """
        Helper function to create a list of sequences increasing SPAM blocks
        with a given angle.
        """
        SPAMBlock = [X(qubit), U(qubit, phase=pi / 2 + angle), X(qubit),
                     U(qubit, phase=pi / 2 + angle)]
        return [[Y90(qubit)] + SPAMBlock * rep + [X90(qubit)]
                for rep in range(maxSpamBlocks)]

    #Insert an identity at the start of every set to mark them off
    seqs = list(chain.from_iterable([[[Id(qubit)]] + spam_seqs(angle)
                                     for angle in angleSweep]))

    #Add a final pi for reference
    seqs.append([X(qubit)])

    #Add the measurment block to every sequence
    measBlock = MEAS(qubit)
    for seq in seqs:
        seq.append(measBlock)

    metafile = compile_to_hardware(seqs, 'SPAM/SPAM')

    if showPlot:
        plot_pulse_files(metafile)
    return metafile
