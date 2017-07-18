"""
Sequences for optimizing gating timing.
"""
from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware


def sweep_gateDelay(qubit, sweepPts):
    """
    Sweep the gate delay associated with a qubit channel using a simple Id, Id, X90, X90
    seqeuence.
    
    Parameters
    ---------
    qubit : logical qubit to create sequences for
    sweepPts : iterable to sweep the gate delay over.
    """

    generator = qubit.phys_chan.generator
    oldDelay = generator.gateDelay

    for ct, delay in enumerate(sweepPts):
        seqs = [[Id(qubit, length=120e-9), Id(qubit), MEAS(qubit)],
                [Id(qubit, length=120e-9), MEAS(qubit)],
                [Id(qubit, length=120e-9), X90(qubit), MEAS(qubit)],
                [Id(qubit, length=120e-9), X90(qubit), MEAS(qubit)]]

        generator.gateDelay = delay

        compile_to_hardware(seqs,
                            'BlankingSweeps/GateDelay',
                            suffix='_{}'.format(ct + 1))

    generator.gateDelay = oldDelay
