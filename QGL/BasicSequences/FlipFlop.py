from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from itertools import chain

from typing import Iterable

def FlipFlop(qubit: Channels.LogicalChannel,
             dragParamSweep: Iterable[float], 
             maxNumFFs: int = 10, 
             showPlot: bool = False) -> str:
    """
    Flip-flop sequence (X90-X90m)**n to determine off-resonance or DRAG
    parameter optimization.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    dragParamSweep : float iterable
        Array-like drag parameter values to sweep over
    maxNumFFs : int, optional
        Maximum number of flip-flop pairs (default 10)
    showPlot : bool, optional
        Whether to plot (boolean)

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> mf = FlipFlop(q2, np.linspace(0.1,0.9))
    Compiled 551 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
	"""

    def flipflop_seqs(drag_scaling):
        """ Helper function to create a list of sequences with a specified drag parameter. """
        qubit.pulse_params['drag_scaling'] = drag_scaling
        return [[X90(qubit)] + [X90(qubit), X90m(qubit)] * rep + [Y90(qubit)]
                for rep in range(maxNumFFs)]

    #Insert an identity at the start of every set to mark them off
    originalScaling = qubit.pulse_params['drag_scaling']
    seqs = list(chain.from_iterable([[[Id(qubit)]] + flipflop_seqs(dragParam)
                                     for dragParam in dragParamSweep]))
    qubit.pulse_params['drag_scaling'] = originalScaling

    #Add a final pi for reference
    seqs.append([X(qubit)])

    #Add the measurment block to every sequence
    measBlock = MEAS(qubit)
    for seq in seqs:
        seq.append(measBlock)

    metafile = compile_to_hardware(seqs, 'FlipFlop/FlipFlop')

    if showPlot:
        plot_pulse_files(metafile)

    return metafile
