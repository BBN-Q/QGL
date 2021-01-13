from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from .helpers import create_cal_seqs
from itertools import product
import operator
from ..ControlFlow import *
from ..TdmInstructions import *
from functools import reduce
from typing import Iterable, Union, Tuple

@qfunction
def qreset(qubits: Channels.LogicalChannel, 
           signVec: Tuple[bool], 
           measDelay: Union[int,float], 
           buf: Union[int,float], 
           reg_size: int = None, 
           TDM_map: Iterable[Union[int,bool]] = None) -> list:
    """
    For each qubit, build the set of feedback actions to perform when
    receiving a zero or one in the comparison register

    Parameters
    ----------
    qubits : Channels.LogicalChannel tuple
        A hashable (immutable) tuple of qubits to reset
    signVec : boolean tuple
        A hashable (immutable) tuple of binary values from the compairison
        register indicating the measured state of each qubit in the register
        before reset.
    measDelay : int/float
        Delay after measurement before performing the LOADCMP comparison with
        value in the register (seconds)
    buf : int/float
        Wait time between (seconds)
    reg_size : int, optional
        Size of the register in number of qubits, including those not reset.
        Default value is set to len(qubits).
    TDM_map : bit mask, optional
        Map each qubit to a TDM digital input.  If True, arguments reset a
        subset of the qubit register (see Reset).
        Default: np.array(qN, qN-1, ..., q1) from MSB to LSB.

    Returns
    -------
    seq : QGL.ControlFlow.Call
        QGL sequence with the qreset calls

    Examples
    --------
    >>> qreset((q1, q2), (0,1), 2e-6, 2e-6);
    CALL(H:)
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
        # duplicate branches for the irrelevant results
        # if reg_size > len(TDM_map)
        meas_result = [(ct & TDM_bit)>0 for TDM_bit in 2**(np.array(TDM_map)-1)]
        branch_idx = sum([t*2**(len(qubits)-ind-1)
                          for ind,t in enumerate((meas_result))])
        seq += qif(ct, [FbSeq[branch_idx]])

    return seq


def Reset(qubits: Iterable[Channels.LogicalChannel],
          measDelay: Union[int,float]=1e-6,
          signVec: Tuple[bool] = None,
          doubleRound: bool = True,
          buf: Union[int,float] = 20e-9,
          showPlot: bool = False,
          measChans: Channels.LogicalChannel = None,
          add_cals: bool = True,
          calRepeats: int = 2,
          reg_size: int = None,
          TDM_map: Iterable[Union[int,bool]]=None) -> str:
    """
    Preparation, simultanoeus reset, and measurement of an arbitrary number
    of qubits

    Parameters
    ----------
    qubits : Channels.LogicalChannel tuple
        A hashable (immutable) tuple of qubits to reset
    measDelay : int/float, optional
        Delay after measurement before performing the LOADCMP compairison with
        value in the register (seconds)
    signVec : boolean tuple, optional
        conditions for feedback. Tuple of 0 (flip if signal is above threshold) and 1 (flip if below) for each qubit. Default = 0 for all qubits
    doubleRound : boolean, optional
        If true, do two rounds of feedback
    showPlot : boolean, optional
        Whether to plot
    measChans : LogicalChannel tuple, optional
        A hashable (immutable) tuple of qubits to measured.
    add_cals : boolean, optional
        Whether to append calibration pulses to the end of the sequence
    calRepeats : int, optional
        How many times to repeat calibration scalings (default 2)
    reg_size : int, optional
        Size of the register in number of qubits, including those not reset.
        Default value is set to len(qubits).
    TDM_map : bit mask, optional
        Map each qubit to a TDM digital input.  If True, arguments reset a
        subset of the qubit register (see Reset).
        Default: np.array(qN, qN-1, ..., q1) from MSB to LSB.

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths to
        compiled machine files

    Examples
    --------
    >>> Reset((q1, q2));
    Compiled 12 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """
    if measChans is None:
        measChans = qubits

    if signVec == None:
        signVec = (0, ) * len(qubits)

    seqs = [prep + [qreset(qubits,
                           signVec,
                           measDelay,
                           buf,
                           reg_size=reg_size,
                           TDM_map=TDM_map)]
                           for prep in create_cal_seqs(qubits, 1)]
    measBlock = reduce(operator.mul, [MEAS(q) for q in qubits])
    if doubleRound:
        for seq in seqs:
            seq += [measBlock]
            seq.append(qreset(qubits,
                              signVec,
                              measDelay,
                              buf,
                              reg_size=reg_size,
                              TDM_map=TDM_map))

    # add final measurement
    for seq in seqs:
        seq += [measBlock, Id(qubits[0], measDelay), qwait(kind='CMP')]

    if add_cals:
        seqs += create_cal_seqs(qubits,
                                calRepeats,
                                measChans=measChans,
                                waitcmp=True)

    metafile = compile_to_hardware(seqs, 'Reset/Reset')

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


# do not make it a subroutine for now
def BitFlip3(data_qs: Iterable[Channels.LogicalChannel],
             ancilla_qs: Iterable[Channels.LogicalChannel],
             theta: Union[int,float] = None,
             phi: Union[int,float] = None,
             nrounds: int = 1,
             meas_delay: Union[int,float] = 1e-6,
             add_cals: bool = False,
             calRepeats: int = 2) -> str:
    """
    Encoding on 3-qubit bit-flip code, followed by n rounds of syndrome
    detection, and final correction using the n results.

    Parameters
    ----------
    data_qs : Channels.LogicalChannel tuple
        A hashable (immutable) tuple of qubits of the 3 code qubits
    ancilla_qs : Channels.LogicalChannel tuple
        A hashable (immutable) tuple of qubits of the 2 syndrome qubits
    theta : int/float, optional
        Longitudinal rotation angle for the encoded state (radians).
        Default = None.
    phi : int/float, optional
        Azimuthal rotation angle for the encoded state (radians).
        Default = None.
    nrounds: int, optional
        Number of consecutive measurements
    measDelay : int/float, optional
        Delay between syndrome check rounds (seconds)
    add_cals : boolean, optional
        Whether to append calibration pulses to the end of the sequence
    calRepeats : int, optional
        How many times to repeat calibration scalings (default 2)

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths to
        compiled machine files

    Examples
    --------
    >>> mf = BitFlip3((q1, q2, q3), (q4, q5));
    Compiled 12 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """
    if len(data_qs) != 3 or len(ancilla_qs) != 2:
        raise Exception("Wrong number of qubits")
    seqs =  [
    DecodeSetRounds(1,0,nrounds),
    Invalidate(10, 2*nrounds),
    Invalidate(11, 0x1)]

    # encode single-qubit state into 3 qubits
    if theta and phi:
        seqs+=[Utheta(data_qs[1], theta, phi),
               CNOT(data_qs[1], data_qs[0]),
               CNOT(data_qs[1], data_qs[2])]

    # multiple rounds of syndrome measurements
    for n in range(nrounds):
        seqs+= [CNOT(data_qs[0],ancilla_qs[0])*CNOT(data_qs[1],ancilla_qs[1])],
        seqs+= [CNOT(data_qs[1], ancilla_qs[0])*CNOT(data_qs[2],ancilla_qs[1])],
        seqs+= [MEASA(ancilla_qs[0], maddr=(10, 2*n))*
                MEASA(ancilla_qs[1], maddr=(10, 2*n+1)),
                Id(ancilla_qs[0], meas_delay),
                    MEAS(data_qs[0], amp=0)*
                    MEAS(data_qs[1], amp=0)*
                    MEAS(data_qs[2], amp=0)]
                    # virtual msmt's just to keep the number of segments
                    # uniform across digitizer channels
    seqs+=Decode(10, 11, 2*nrounds)
    seqs+=qwait("RAM",11)
    seqs+=[MEAS(data_qs[0])*
           MEAS(data_qs[1])*
           MEAS(data_qs[2])*
           MEAS(ancilla_qs[0], amp=0)*
           MEAS(ancilla_qs[1], amp=0)]
           # virtual msmt's

    # apply corrective pulses depending on the decoder result
    FbGates = []
    for q in data_qs:
        FbGates.append([gate(q) for gate in [Id, X]])
    FbSeq = [reduce(operator.mul, x) for x in product(*FbGates)]
    for k in range(8):
        seqs += qif(k, [FbSeq[k]])
    if add_cals:
        seqs += create_cal_seqs(qubits,
        calRepeats)
    metafile = compile_to_hardware(seqs, 'BitFlip/BitFlip', tdm_seq=True)
    return metafile

def MajorityVoteN(qubits: Iterable[Channels.LogicalChannel],
                  nrounds: int,
                  prep: Iterable[bool] = [],
                  meas_delay: float = 1e-6,
                  add_cals: bool = False,
                  calRepeats: int = 2) -> str:
    """
    Majority vote across multiple measurement results (same or different qubits)

    Parameters
    ----------
    qubits : Channels.LogicalChannel tuple
        A hashable (immutable) tuple of qubits for majority vote
    nrounds: int
        Number of consecutive measurements
    prep : boolean iterable, optional
        Array of binary values mapping X(q) pulses to the list of qubits
        proivided. Ex: (q1,q2), prep=(1,0) -> would apply a pi pulse to q1
        before the majority vote measurement. Default = []
    measDelay : int/float, optional
        Delay between syndrome check rounds (seconds)
    add_cals : boolean, optional
        Whether to append calibration pulses to the end of the sequence
    calRepeats : int, optional
        How many times to repeat calibration scalings (default 2)

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths to
        compiled machine files

    Examples
    --------
    >>> mf = MajorityVoteN((q1, q2, q3), 10);
    Compiled 1 sequences.
    o INVALIDATE(channel=None, addr=0x1, mask=0x0)
    o WRITEADDR(channel=None, addr=0x1, value=0xfffff)
    MAJORITYMASK(in_addr=1, out_addr=0)
    o INVALIDATE(channel=None, addr=0xa, mask=0xfffff)
    o INVALIDATE(channel=None, addr=0xb, mask=0x1)
    MAJORITY(in_addr=a, out_addr=b)
    >>> mf
    '/path/to/exp/exp-meta.json'
    """
    nqubits = len(qubits)
    seqs = [MajorityMask(1, 0, nrounds*nqubits),
           Invalidate(10, nrounds*nqubits),
           Invalidate(11, 1)]
    if prep:
       seqs += [reduce(operator.mul,
                    [X(q) for n,q in enumerate(qubits) if prep[n]])]
    for n in range(nrounds):
       seqs += [reduce(operator.mul,
               [MEASA(q, (10, nqubits*n+m)) for m,q in enumerate(qubits)]),
               Id(qubits[0],meas_delay)]
    seqs+=MajorityVote(10,11, nrounds*nqubits)
    seqs+=qwait("RAM", 11)
    seqs+=[Id(qubits[0],100e-9)]
    seqs+=qif(1,[X(qubits[0])]) # placeholder for any conditional operation
    seqs=[seqs]
    if add_cals:
        seqs += create_cal_seqs(qubits,
        calRepeats)
    metafile = compile_to_hardware(seqs,
                                   'MajorityVote/MajorityVote',
                                   tdm_seq=True)
    return metafile
