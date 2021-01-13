"""
Sequences for optimizing gating timing.
"""
from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from typing import Iterable, Union


def sweep_gateDelay(qubit: Channels.LogicalChannel, 
                    sweepPts: Iterable[Union[int,float]]) -> str:
    """
    Sweep the gate delay associated with a qubit channel using a simple Id, Id,
    X90, X90 seqeuence.

    Parameters
    ---------
    qubit : Channels.LogicalChannel
        Qubit channel for which to create sequences
    sweepPts : int/float iterable
        Iterable to sweep the gate delay over (seconds)

    Returns
    -------
    void : string
        This functions produces a set of files enumerating the sweepPts
        given in the parameters and returns nothing.  This function is currently
        not used and will be depricated in the future.

    Examples
    --------
    >>> sweepgateDelay(q1, np.linspace(20.0e-9, 220.0e-9, 101))
    """

    generator = qubit.phys_chan.generator
    oldDelay = generator.gateDelay

    for ct, delay in enumerate(sweepPts):
        seqs = [[Id(qubit, length=120e-9), Id(qubit), MEAS(qubit)],
                [Id(qubit, length=120e-9), MEAS(qubit)],
                [Id(qubit, length=120e-9), X90(qubit), MEAS(qubit)],
                [Id(qubit, length=120e-9), X90(qubit), MEAS(qubit)]]

        generator.gateDelay = delay

        metafile = compile_to_hardware(seqs, 
                                       'BlankingSweeps/GateDelay', 
                                       suffix='_{}'.format(ct + 1))

    generator.gateDelay = oldDelay

    return metafile
