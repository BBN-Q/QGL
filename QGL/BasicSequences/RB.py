from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from ..Cliffords import clifford_seq, clifford_mat, inverse_clifford
from .helpers import create_cal_seqs

import os
from csv import reader
import numpy as np
from functools import reduce


def create_RB_seqs(numQubits, lengths, repeats=32, interleaveGate=None, recovery=True):
    """Create a list of lists of Clifford gates to implement RB. """
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

def SingleQubitRB(qubit, seqs, purity=False, showPlot=False):
    """Single qubit randomized benchmarking using 90 and 180 generators.

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    seqs : list of lists of Clifford group integers
    showPlot : whether to plot (boolean)
    """

    seqsBis = []
    op = [Id(qubit, length=0), Y90m(qubit), X90(qubit)]
    for ct in range(3 if purity else 1):
        for seq in seqs:
            seqsBis.append(reduce(operator.add, [clifford_seq(c, qubit)
                                                for c in seq]))
            #append tomography pulse to measure purity
            seqsBis[-1].append(op[ct])
            #append measurement
            seqsBis[-1].append(MEAS(qubit))

    #Tack on the calibration sequences
    seqsBis += create_cal_seqs((qubit, ), 2)

    metafile = compile_to_hardware(seqsBis, 'RB/RB')

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def TwoQubitRB(q1, q2, seqs, showPlot=False, suffix=""):
    """Two qubit randomized benchmarking using 90 and 180 single qubit generators and ZX90

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    seqs : list of lists of Clifford group integers
    showPlot : whether to plot (boolean)
    suffix : suffix to apply to sequence file names
    """
    seqsBis = []
    for seq in seqs:
        seqsBis.append(reduce(operator.add, [clifford_seq(c, q1, q2)
                                             for c in seq]))

    #Add the measurement to all sequences
    for seq in seqsBis:
        seq.append(MEAS(q1) * MEAS(q2))

    #Tack on the calibration sequences
    seqsBis += create_cal_seqs((q1, q2), 2)

    metafile = compile_to_hardware(seqsBis, 'RB/RB', suffix=suffix)

    if showPlot:
        plot_pulse_files(metafile)
    return metafile

def SingleQubitRB_AC(qubit, seqs, purity=False, showPlot=False):
    """Single qubit randomized benchmarking using atomic Clifford pulses.

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    seqFile : file containing sequence strings
    showPlot : whether to plot (boolean)
    """
    seqsBis = []
    op = [Id(qubit, length=0), Y90m(qubit), X90(qubit)]
    for ct in range(3 if purity else 1):
        for seq in seqs:
            seqsBis.append([AC(qubit, c) for c in seq])
            #append tomography pulse to measure purity
            seqsBis[-1].append(op[ct])
            #append measurement
            seqsBis[-1].append(MEAS(qubit))

    #Tack on the calibration sequences
    seqsBis += create_cal_seqs((qubit, ), 2)

    metafile = compile_to_hardware(seqsBis, 'RB/RB')

    if showPlot:
        plot_pulse_files(metafile)
    return metafile

def SingleQubitRB_DiAC(qubit, seqs, compiled=True, purity=False, showPlot=False):
    """Single qubit randomized benchmarking using diatomic Clifford pulses.

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    seqFile : file containing sequence strings
    compiled : if True, compile Z90(m)-X90-Z90(m) to Y90(m) pulses
    purity : measure <Z>,<X>,<Y> of final state, to measure purity. See J.J.
        Wallman et al., New J. Phys. 17, 113020 (2015)
    showPlot : whether to plot (boolean)
    """
    seqsBis = []
    op = [Id(qubit, length=0), Y90m(qubit), X90(qubit)]
    for ct in range(3 if purity else 1):
        for seq in seqs:
            seqsBis.append([DiAC(qubit, c, compiled) for c in seq])
            #append tomography pulse to measure purity
            seqsBis[-1].append(op[ct])
            #append measurement
            seqsBis[-1].append(MEAS(qubit))

    #Tack on the calibration sequences (using pi/2 pulses for consistency)
    seqsBis += [[Id(qubit), MEAS(qubit)], [Id(qubit), MEAS(qubit)], [X90(qubit), X90(qubit), MEAS(qubit)], [X90(qubit), X90(qubit), MEAS(qubit)]]

    metafile = compile_to_hardware(seqsBis, 'RB_DiAC/RB_DiAC')

    if showPlot:
        plot_pulse_files(metafile)
    return metafile

def SingleQubitIRB_AC(qubit, seqFile, showPlot=False):
    """Single qubit interleaved randomized benchmarking using atomic Clifford pulses.

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    seqFile : file containing sequence strings
    showPlot : whether to plot (boolean)
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
    #We've shuffled the sequences so that we loop through each gate length on the inner loop
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


def SingleQubitRBT(qubit, seqFileDir, analyzedPulse, showPlot=False):
    """	Single qubit randomized benchmarking tomography using atomic Clifford pulses.

    This relies on specific sequence files and is here for historical purposes only.

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    seqFile : file containing sequence strings
    analyzedPulse : specific pulse to analyze
    showPlot : whether to plot (boolean)
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
        numCals = 4
        chunk += [[Id(qubit), measBlock]] * numCals + [[X(qubit), measBlock]
                                                       ] * numCals
        metafile = compile_to_hardware(chunk,
                                        'RBT/RBT',
                                        suffix='_{0}'.format(ct + 1))

    if showPlot:
        plot_pulse_files(metafile)
    return metafile


def SimultaneousRB_AC(qubits, seqs, showPlot=False):
    """
    Simultaneous randomized benchmarking on multiple qubits using atomic Clifford pulses.

    Parameters
    ----------
    qubits : iterable of logical channels to implement seqs on (list or tuple)
    seqs : a tuple of sequences created for each qubit in qubits
    showPlot : whether to plot (boolean)

    Example
    -------
    >>> q1 = QubitFactory('q1')
    >>> q2 = QubitFactory('q2')
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

    #Tack on the calibration sequences
    seqsBis += create_cal_seqs((qubits), 2)

    metafile = compile_to_hardware(seqsBis, 'RB/RB')

    if showPlot:
        plot_pulse_files(metafile)
    return metafile
