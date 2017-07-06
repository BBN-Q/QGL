from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from .helpers import create_cal_seqs
from itertools import product
import operator
from ..ControlFlow import *
from functools import reduce


@qfunction
def qreset(qubits, signVec, measDelay, buf, reg_size=None, TDM_map=None):
    """
    for each qubit, build the set of feedback actions to perform when
    receiving a zero or one in the comparison register
    reg_size, TDM_map: optional arguments to reset a subset of the qubit register (see Reset)
    """
    if not reg_size:
        reg_size = len(qubits)
        TDM_map = np.arange(reg_size,0,-1)

    FbGates = []
    for ct, q in enumerate(qubits):
        if signVec[ct] == 0:
            FbGates.append([gate(q) for gate in [Id, X]])
        else:  # inverted logic
            FbGates.append([gate(q) for gate in [X, Id]])
    FbSeq = [reduce(operator.mul, x) for x in product(*FbGates)]

    # load register
    seq = [Id(qubits[0], measDelay), qwait(kind='CMP'), Id(qubits[0], buf)]
    # create a branch for each possible comparison value
    for ct in range(2**reg_size):
        # duplicate branches for the irrelevant results if reg_size > len(TDM_map)
        meas_result = [(ct & TDM_bit)>0 for TDM_bit in 2**(np.array(TDM_map)-1)]
        branch_idx = sum([t*2**(len(qubits)-ind-1) for ind,t in enumerate((meas_result))])
        seq += qif(ct, [FbSeq[branch_idx]])
        
    return seq


def Reset(qubits,
          measDelay=1e-6,
          signVec=None,
          doubleRound=True,
          buf=20e-9,
          showPlot=False,
          measChans=None,
          docals=True,
          calRepeats=2,
          reg_size=None,
          TDM_map=None):
    """

    Preparation, simultanoeus reset, and measurement of an arbitrary number of qubits

    Parameters
    ----------
    qubits : tuple of logical channels to implement sequence (LogicalChannel)
    measDelay : delay between end of measuerement and LOADCMP
    signVec : conditions for feedback. Tuple of 0 (flip if signal is above threshold) and 1 (flip if below) for each qubit. Default = 0 for all qubits
    doubleRound : if true, double round of feedback
    showPlot : whether to plot (boolean)
    measChans : tuble of qubits to be measured (LogicalChannel)
    docals, calRepeats: enable calibration sequences, repeated calRepeats times
    reg_size: total number of qubits, including those that are not reset. Default set to len(qubits)
    TDM_map: map each qubit to a TDM digital input. Default: np.array(qN, qN-1, ..., q1) from MSB to LSB.
    
    Returns
    -------
    plotHandle : handle to plot window to prevent destruction
    """
    if measChans is None:
        measChans = qubits

    if signVec == None:
        signVec = (0, ) * len(qubits)

    seqs = [prep + [qreset(qubits, signVec, measDelay, buf, reg_size=reg_size, TDM_map=TDM_map)]
            for prep in create_cal_seqs(qubits, 1)]
    measBlock = reduce(operator.mul, [MEAS(q) for q in qubits])
    if doubleRound:
        for seq in seqs:
            seq += [measBlock]
            seq.append(qreset(qubits, signVec, measDelay, buf, reg_size=reg_size, TDM_map=TDM_map))

    # add final measurement
    for seq in seqs:
        seq += [measBlock, Id(qubits[0], measDelay), qwait(kind='CMP')]

    if docals:
        seqs += create_cal_seqs(qubits,
                                calRepeats,
                                measChans=measChans,
                                waitcmp=True)

    metafile = compile_to_hardware(seqs, 'Reset/Reset')

    if showPlot:
        plot_pulse_files(metafile)
    return metafile
