from ..PulsePrimitives import *
from ..Cliffords import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from ..tools.clifford_tools import clifford_mat, inverse_clifford
from .helpers import create_cal_seqs, cal_descriptor
from ..config import logger

import os
from csv import reader
import numpy as np
from functools import reduce


def create_RB_seqs(numQubits,
                   lengths,
                   repeats=32,
                   interleaveGate=None,
                   recovery=True):
    """
    Create a list of lists of Clifford gates to implement RB.

    Parameters
    ----------
    numQubits : int
        Number of qubits to create sequences for
    lengths : int iterable
        Array-like list of integers that denote the specific length for the
        various RB experiments.  A common examples would be powers of two
        spacing = [2, 4, 8, 16, 32, 64, ...]
    repeats : int, optional
        Number of individual randomizations for each number of lengths.
        Default = 32.
    interleaveGate : int, optional
        This is the index of a Clifford operation that can be optionally
        interleaved in the sequence if you would like to do interleaved RB.
        The index corresponds to the mapping in QGL.Cliffords.
    recovery : boolean, optional
        Optional parameter which, if false, leaves off the recovery operation

    Returns
    -------
    seq : int list of lists
        A list of lists containing integer pulses indicies based on those in
        QGL.Cliffords.

    Examples
    --------
    >>> create_RB_seqs(1, [2,4,8], repeats=2, interleaveGate=1)
    [[19, 1, 6],
     [3, 1, 0],
     [1, 1, 18, 1, 9, 1, 15],
     [11, 1, 8, 1, 19, 1, 20],
     [6, 1, 21, 1, 23, 1, 8, 1, 13, 1, 2, 1, 3, 1, 0],
     [2, 1, 8, 1, 4, 1, 10, 1, 18, 1, 20, 1, 10, 1, 19]]
    """
    if numQubits == 1:
        cliffGroupSize = 24
    elif numQubits == 2:
        cliffGroupSize = 11520
    else:
        raise Exception("Can only handle one or two qubits.")

    #Create lists of of random integers
    #Subtract one from length for recovery gate
    seqs = []
    for length in lengths:
        seqs += np.random.randint(0, cliffGroupSize,
                                  size=(repeats, length - 1)).tolist()

    #Possibly inject the interleaved gate
    if interleaveGate:
        newSeqs = []
        for seq in seqs:
            newSeqs.append(np.vstack((np.array(
                seq, dtype=np.int), interleaveGate * np.ones(
                    len(seq), dtype=np.int))).flatten(order='F').tolist())
        seqs = newSeqs

    if recovery:
        #Calculate the recovery gate
        for seq in seqs:
            if len(seq) == 1:
                mat = clifford_mat(seq[0], numQubits)
            else:
                mat = reduce(lambda x, y: np.dot(y, x),
                             [clifford_mat(c, numQubits) for c in seq])
            seq.append(inverse_clifford(mat))

    return seqs

def SingleQubitRB(qubit, seqs, cliff_type='std', purity=False, showPlot=False, add_cals=True):
    """
    Single qubit randomized benchmarking using 90 and 180 generators.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    seqs : int iterable
        list of lists of Clifford group integers produced by create_RB_seqs
    purity : boolean, optional
        If True, this create sequences for purity RB
    showPlot : boolean, optional
        Whether to plot
    add_cals : boolean, optional
        Whether to append calibration pulses to the end of the sequence

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> seqs = create_RB_seqs(1, [2,4,8], repeats=2, interleaveGate=1);
    >>> mf = SingleQubitRB(q1, seqs);
    Compiled 10 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    if cliff_type.upper() not in clifford_map.keys():
        raise ValueError(f"Unknown clifford type: must be one of {clifford.map.keys()}.")

    clifford = clifford_map[cliff_type.upper()]

    seqsBis = []
    op = [Id(qubit, length=0), Y90m(qubit), X90(qubit)]
    for ct in range(3 if purity else 1):
        for seq in seqs:
            seqsBis.append([clifford(qubit,c) for c in seq])
            #append tomography pulse to measure purity
            seqsBis[-1].append(op[ct])
            #append measurement
            seqsBis[-1].append(MEAS(qubit))

    axis_descriptor = [{
        'name': 'length',
        'unit': None,
        'points': list(map(len, seqs)),
        'partition': 1
    }]

    #Tack on the calibration sequences
    if add_cals:
        seqsBis += create_cal_seqs((qubit, ), 2)
        axis_descriptor.append(cal_descriptor((qubit,), 2))

    metafile = compile_to_hardware(seqsBis, 'RB/RB', axis_descriptor = axis_descriptor, extra_meta = {'sequences':seqs})

    if showPlot:
        plot_pulse_files(metafile)
    return metafile

def SingleQubitLeakageRB(qubit, seqs, pi2args, cliff_type='std', showPlot=False):
    """
    Single qubit randomized benchmarking using 90 and 180 generators to
    measure leakage outside the qubit subspace.
    See https://journals.aps.org/prl/supplemental/10.1103/
    PhysRevLett.123.120502/Rol_SOM.pdf for description of algorithm.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    seqs : int iterable
        list of lists of Clifford group integers produced by create_RB_seqs
    pi2args: dictionary mapping
        Arguments passed to the X90 gate for the 1 <-> 2 transition during
        calibration
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> seqs = create_RB_seqs(1, [2,4,8]);
    >>> mf = SingleQubitLeakageRB(q1, seqs, {'one': 1, 'two': 2});
    Compiled 10 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    if cliff_type.upper() not in clifford_map.keys():
        raise ValueError(f"Unknown clifford type: must be one of {clifford.map.keys()}.")

    clifford = clifford_map[cliff_type.upper()]

    seqsBis = []
    for seq in seqs:
        combined_seq = [clifford(qubit, c) for c in seq]

        # Append sequence with tomography ids and measurement
        seqsBis.append(combined_seq + [Id(qubit), Id(qubit), MEAS(qubit)])

        # Append sequence with tomography pulses and measurement
        seqsBis.append(combined_seq + [X90(qubit), X90(qubit), MEAS(qubit)])

    # Add the calibration sequences
    seqsBis.append([Id(qubit), Id(qubit), Id(qubit), Id(qubit), MEAS(qubit)])
    seqsBis.append([X90(qubit), X90(qubit), Id(qubit), Id(qubit), MEAS(qubit)])
    seqsBis.append([X90(qubit), X90(qubit), X90(qubit, **pi2args), X90(qubit, **pi2args), MEAS(qubit)])

    axis_descriptor = [
    {
        'name': 'length',
        'unit': None,
        'points': [len(s) for s in seqs for i in range(2)],
        'partition': 1
    },
    {
        'name': 'calibration',
        'unit': 'state',
        'partition': 2,
        'points': ['0', '1', '2']
    }]

    metafile = compile_to_hardware(seqsBis, 'RB/LRB', axis_descriptor = axis_descriptor, extra_meta = {'sequences':seqs})

    if showPlot:
        plot_pulse_files(metafile)
    return metafile



def TwoQubitRB(q1, q2, seqs, cliff_type='std', showPlot=False, suffix="", add_cals=True):
    """
    Two qubit randomized benchmarking using 90 and 180 single qubit generators
    and ZX90.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement RB
    qubit : Channels.LogicalChannel
        Logical channel to implement RB
    seqs : int iterable
        list of lists of Clifford group integers produced by create_RB_seqs
    showPlot : boolean, optional
        Whether to plot
    suffix : string, optional
        Suffix to add to generated files
    add_cals : boolean, optional
        Whether to append calibration pulses to the end of the sequence

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> seqs = create_RB_seqs(2, [2,4,8], repeats=2);
    >>> mf = TwoQubitRB(q1, q2, seqs);
    Compiled 14 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    seqsBis = []
    for seq in seqs:
        seqsBis.append(reduce(operator.add, [TwoQubitClifford(q2, q1, kind=cliff_type)
                                             for c in seq]))

    #Add the measurement to all sequences
    for seq in seqsBis:
        seq.append(MEAS(q1) * MEAS(q2))

    axis_descriptor = [{
        'name': 'length',
        'unit': None,
        'points': list(map(len, seqs)),
        'partition': 1
    }]

    #Tack on the calibration sequences
    if add_cals:
        seqsBis += create_cal_seqs((q1, q2), 2)
        axis_descriptor.append(cal_descriptor((q1, q2), 2))

    metafile = compile_to_hardware(seqsBis, 'RB/RB', axis_descriptor = axis_descriptor, suffix = suffix, extra_meta = {'sequences':seqs})

    if showPlot:
        plot_pulse_files(metafile)
    return metafile

def TwoQubitLeakageRB(q1, q2, meas_qubit, seqs, pi2args, cliff_type='std', showPlot=False):
    """
    Two qubit randomized benchmarking using 90 and 180 single qubit generators
    and ZX90 to measure leakage outside the qubit subspace.  See https://
    journals.aps.org/prl/supplemental/10.1103/PhysRevLett.123.120502/Rol_SOM.pdf
    for description of algorithm.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement RB
    qubit : Channels.LogicalChannel
        Logical channel to implement RB
    meas_qubit : Channels.LogicalChannel
        Qubit to measure
    seqs : int iterable
        list of lists of Clifford group integers produced by create_RB_seqs
    pi2args: dictionary mapping
        Arguments passed to the X90 gate for the 1 <-> 2 transition during
        calibration
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> seqs = create_RB_seqs(2, [2,4,8], repeats=2);
    >>> mf = TwoQubitLeakageRB(q1, q2, q1, seqs, {'one': 1, 'two': 2});
    Compiled 14 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """
    seqsBis = []
    for seq in seqs:
        combined_seq = reduce(operator.add, [TwoQubitClifford(q2, q1, c, kind=cliff_type) for c in seq])

        # Append sequence with tomography ids and measurement
        seqsBis.append(combined_seq + [Id(meas_qubit), Id(meas_qubit), MEAS(meas_qubit)])

        # Append sequence with tomography pulses and measurement
        seqsBis.append(combined_seq + [X90(meas_qubit), X90(meas_qubit), MEAS(meas_qubit)])

    # Add the calibration sequences
    seqsBis.append([Id(meas_qubit), Id(meas_qubit), Id(meas_qubit), Id(meas_qubit), MEAS(meas_qubit)])
    seqsBis.append([X90(meas_qubit), X90(meas_qubit), Id(meas_qubit), Id(meas_qubit), MEAS(meas_qubit)])
    seqsBis.append([X90(meas_qubit), X90(meas_qubit), X90(meas_qubit, **pi2args), X90(meas_qubit, **pi2args), MEAS(meas_qubit)])

    axis_descriptor = [
    {
        'name': 'length',
        'unit': None,
        'points': [len(s) for s in seqs for i in range(2)],
        'partition': 1
    },
    {
        'name': 'calibration',
        'unit': 'state',
        'partition': 2,
        'points': ['0', '1', '2']
    }]

    metafile = compile_to_hardware(seqsBis, 'RB/LRB', axis_descriptor = axis_descriptor, extra_meta = {'sequences':seqs})

    if showPlot:
        plot_pulse_files(metafile)
    return metafile

def SingleQubitRB_AC(qubit, seqs, purity=False, showPlot=False, add_cals=True):
    """
    Single qubit randomized benchmarking using atomic Clifford pulses.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    seqs : int iterable
        list of lists of Clifford group integers produced by create_RB_seqs
    purity : boolean, optional
        If True, this create sequences for purity RB: measure <Z>,<X>,<Y> of
        final state, to measure purity. See J.J. Wallman et al.,
        New J. Phys. 17, 113020 (2015)
    showPlot : boolean, optional
        Whether to plot
    add_cals : boolean, optional
        Whether to append calibration pulses to the end of the sequence

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> seqs = create_RB_seqs(1, [2,4,8], repeats=2, interleaveGate=1);
    >>> mf = SingleQubitRB_AC(q1, seqs);
    Compiled 10 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    logger.warning("This function is deprecated and may be removed in a future release of QGL! " +
        "Use `SingleQubitRB` with the `cliff_type` keyword argument instead.")

    seqsBis = []
    op = [Id(qubit, length=0), Y90m(qubit), X90(qubit)]
    for ct in range(3 if purity else 1):
        for seq in seqs:
            seqsBis.append([AC(qubit, c) for c in seq])
            #append tomography pulse to measure purity
            seqsBis[-1].append(op[ct])
            #append measurement
            seqsBis[-1].append(MEAS(qubit))

    axis_descriptor = [{
        'name': 'length',
        'unit': None,
        'points': list(map(len, seqs)),
        'partition': 1
    }]

    #Tack on the calibration sequences
    if add_cals:
        seqsBis += create_cal_seqs((qubit, ), 2)
        axis_descriptor.append(cal_descriptor((qubit,), 2))

    metafile = compile_to_hardware(seqsBis, 'RB/RB', axis_descriptor = axis_descriptor, extra_meta = {'sequences':seqs})

    if showPlot:
        plot_pulse_files(metafile)
    return metafile

def SingleQubitRB_DiAC(qubit,
                       seqs,
                       compiled=True,
                       purity=False,
                       showPlot=False,
                       add_cals=True):
    """
    Single qubit randomized benchmarking using diatomic Clifford pulses.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    seqs : int iterable
        list of lists of Clifford group integers produced by create_RB_seqs
    compiled : boolean, optional
        If True, compile Z90(m)-X90-Z90(m) to Y90(m) pulses
    purity : boolean, optional
        If True, this create sequences for purity RB: measure <Z>,<X>,<Y> of
        final state, to measure purity. See J.J. Wallman et al.,
        New J. Phys. 17, 113020 (2015)
    showPlot : boolean, optional
        Whether to plot
    add_cals : boolean, optional
        Whether to append calibration pulses to the end of the sequence

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> seqs = create_RB_seqs(1, [2,4,8], repeats=2, interleaveGate=1);
    >>> mf = SingleQubitRB_DiAC(q1, seqs);
    Compiled 10 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    logger.warning("This function is deprecated and may be removed in a future release of QGL! " +
        "Use `SingleQubitRB` with the `cliff_type` keyword argument instead.")

    seqsBis = []
    op = [Id(qubit, length=0), Y90m(qubit), X90(qubit)]
    for ct in range(3 if purity else 1):
        for seq in seqs:
            seqsBis.append([DiAC(qubit, c, compiled) for c in seq])
            #append tomography pulse to measure purity
            seqsBis[-1].append(op[ct])
            #append measurement
            seqsBis[-1].append(MEAS(qubit))

    axis_descriptor = [{
        'name': 'length',
        'unit': None,
        'points': list(map(len, seqs)),
        'partition': 1
    }]

    #Tack on the calibration sequences
    if add_cals:
        seqsBis += [[Id(qubit), MEAS(qubit)], [Id(qubit), MEAS(qubit)], [X90(qubit), X90(qubit), MEAS(qubit)], [X90(qubit), X90(qubit), MEAS(qubit)]]
        axis_descriptor.append(cal_descriptor((qubit,), 2))

    metafile = compile_to_hardware(seqsBis, 'RB_DiAC/RB_DiAC', axis_descriptor = axis_descriptor, extra_meta = {'sequences':seqs})

    if showPlot:
        plot_pulse_files(metafile)
    return metafile

def SingleQubitIRB_AC(qubit, seqFile, showPlot=False):
    """
    Single qubit interleaved randomized benchmarking using atomic Clifford
    pulses.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel to implement sequence
    seqsFiles : string
        String defining the path to the file with sequence strings
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Examples
    --------
    >>> seqs = create_RB_seqs(1, [2,4,8], repeats=2, interleaveGate=1);
    >>> mf = SingleQubitIRB_AC(q1, '/path/to/seq/strings/file');
    Compiled 10 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """
    #Setup a pulse library
    pulseLib = [AC(qubit, cliffNum) for cliffNum in range(24)]
    pulseLib.append(pulseLib[0])
    measBlock = MEAS(qubit)

    with open(seqFile, 'r') as FID:
        fileReader = reader(FID)
        seqs = []
        for pulseSeqStr in fileReader:
            seq = []
            for pulseStr in pulseSeqStr:
                seq.append(pulseLib[int(pulseStr)])
            seq.append(measBlock)
            seqs.append(seq)

    #Hack for limited APS waveform memory and break it up into multiple files
    #We've shuffled the sequences so that we loop through each gate length
    #on the inner loop
    numRandomizations = 36
    for ct in range(numRandomizations):
        chunk = seqs[ct::numRandomizations]
        chunk1 = chunk[::2]
        chunk2 = chunk[1::2]
        #Tack on the calibration scalings
        chunk1 += [[Id(qubit), measBlock], [X(qubit), measBlock]]
        metafile = compile_to_hardware(chunk1,
                                        'RB/RB',
                                        suffix='_{0}'.format(2 * ct + 1))
        chunk2 += [[Id(qubit), measBlock], [X(qubit), measBlock]]
        metafile = compile_to_hardware(chunk2,
                                        'RB/RB',
                                        suffix='_{0}'.format(2 * ct + 2))

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def SingleQubitRBT(qubit,
                   seqFileDir,
                   analyzedPulse,
                   showPlot=False,
                   add_cals=True):
    """
    Single qubit randomized benchmarking tomography using atomic Clifford
    pulses.

    This relies on specific sequence files and is here for historical
    purposes only!

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    seqFile : file containing sequence strings
    analyzedPulse : specific pulse to analyze
    showPlot : whether to plot (boolean)

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files
    """
    #Setup a pulse library
    pulseLib = [AC(qubit, cliffNum) for cliffNum in range(24)]
    pulseLib.append(analyzedPulse)
    measBlock = MEAS(qubit)

    seqs = []
    for ct in range(10):
        fileName = 'RBT_Seqs_fast_{0}_F1.txt'.format(ct + 1)
        tmpSeqs = []
        with open(os.path.join(seqFileDir, fileName), 'r') as FID:
            fileReader = reader(FID)
            for pulseSeqStr in fileReader:
                seq = []
                for pulseStr in pulseSeqStr:
                    seq.append(pulseLib[int(pulseStr) - 1])
                seq.append(measBlock)
                tmpSeqs.append(seq)
            seqs += tmpSeqs[:12] * 12 + tmpSeqs[12:-12] + tmpSeqs[-12:] * 12

    seqsPerFile = 100
    numFiles = len(seqs) // seqsPerFile

    for ct in range(numFiles):
        chunk = seqs[ct * seqsPerFile:(ct + 1) * seqsPerFile]
        #Tack on the calibration scalings
        if add_cals:
            numCals = 4
            chunk += [[Id(qubit), measBlock]] * numCals + [[X(qubit), measBlock]
                                                        ] * numCals
        metafile = compile_to_hardware(chunk,
                                        'RBT/RBT',
                                        suffix='_{0}'.format(ct + 1))

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def SimultaneousRB_AC(qubits, seqs, showPlot=False, add_cals=True):
    """
    Simultaneous randomized benchmarking on multiple qubits using atomic
    Clifford pulses.

    Parameters
    ----------
    qubits : Channels.LogicalChannel tuple
        A tuple of two logical channels to implement RB
    seqs : int iterable tuple
        A length two tuple containing list of lists of Clifford group
        integers produced by create_RB_seqs
    showPlot : boolean, optional
        Whether to plot
    add_cals : boolean, optional
        Whether to append calibration pulses to the end of the sequence

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files

    Example
    -------
    >>> seqs1 = create_RB_seqs(1, [2, 4, 8, 16])
    >>> seqs2 = create_RB_seqs(1, [2, 4, 8, 16])
    >>> SimultaneousRB_AC((q1, q2), (seqs1, seqs2), showPlot=False)
    """
    seqsBis = []
    for seq in zip(*seqs):
        seqsBis.append([reduce(operator.__mul__,
                               [AC(q, c) for q, c in zip(qubits, pulseNums)])
                        for pulseNums in zip(*seq)])

    #Add the measurement to all sequences
    for seq in seqsBis:
        seq.append(reduce(operator.mul, [MEAS(q) for q in qubits]))

    axis_descriptor = [{
        'name': 'length',
        'unit': None,
        'points': list(map(len, seqs)),
        'partition': 1
    }]

    #Tack on the calibration sequences
    if add_cals:
        seqsBis += create_cal_seqs((qubits), 2)
        axis_descriptor.append(cal_descriptor((qubits), 2))

    metafile = compile_to_hardware(seqsBis, 'RB/RB', axis_descriptor = axis_descriptor, extra_meta = {'sequences':seqs})

    if showPlot:
        plot_pulse_files(metafile)
    return metafile
