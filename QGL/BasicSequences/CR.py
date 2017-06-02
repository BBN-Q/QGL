from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..ChannelLibrary import EdgeFactory
from ..PulseSequencePlotter import plot_pulse_files
from .helpers import create_cal_seqs, delay_descriptor, cal_descriptor
import numpy as np

def PiRabi(controlQ,
           targetQ,
           lengths,
           riseFall=40e-9,
           amp=1,
           phase=0,
           calRepeats=2,
           showPlot=False):
    """
	Variable length CX experiment.

	Parameters
	----------
	controlQ : logical channel for the control qubit (LogicalChannel)
	targetQ: logical channel for the target qubit (LogicalChannel)
	lengths : pulse lengths of the CR pulse to sweep over (iterable)
	riseFall : rise/fall time of the CR pulse (s)
	amp : amplitude of the CR pulse
	phase : phase of the CR pulse (rad)
	showPlot : whether to plot (boolean)
	"""

    CRchan = EdgeFactory(controlQ, targetQ)
    seqs = [[Id(controlQ),
             flat_top_gaussian(CRchan, riseFall, amp=amp, phase=phase, length=l),
             MEAS(targetQ)*MEAS(controlQ)] for l in lengths] + \
           [[X(controlQ),
             flat_top_gaussian(CRchan, riseFall, amp=amp, phase=phase, length=l),
             X(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for l in lengths] + \
              create_cal_seqs([targetQ,controlQ], calRepeats, measChans=(targetQ,controlQ))

    metafile = compile_to_hardware(seqs, 'PiRabi/PiRabi',
        axis_descriptor=[
            delay_descriptor(np.concatenate((lengths, lengths))),
            cal_descriptor((controlQ, targetQ), calRepeats)
        ])

    if showPlot:
        plot_pulse_files(metafile)

    return metafile


def EchoCRLen(controlQ,
              targetQ,
              lengths,
              riseFall=40e-9,
              amp=1,
              phase=0,
              calRepeats=2,
              showPlot=False):
    """
	Variable length CX experiment, with echo pulse sandwiched between two CR opposite-phase pulses.

	Parameters
	----------
	controlQ : logical channel for the control qubit (LogicalChannel)
	targetQ: logical channel for the target qubit (LogicalChannel)
	lengths : pulse lengths of the CR pulse to sweep over (iterable)
	riseFall : rise/fall time of the CR pulse (s)
	amp : amplitude of the CR pulse
	phase : phase of the CR pulse (rad)
	calRepeats : number of repetitions of readout calibrations for each 2-qubit state
	showPlot : whether to plot (boolean)
	"""
    seqs = [[Id(controlQ),
             echoCR(controlQ, targetQ, length=l, phase=phase, amp=amp, riseFall=riseFall),
             Id(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for l in lengths] + \
           [[X(controlQ),
             echoCR(controlQ, targetQ, length=l, phase= phase, amp=amp, riseFall=riseFall),
             X(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for l in lengths] + \
           create_cal_seqs((targetQ,controlQ), calRepeats, measChans=(targetQ,controlQ))

    metafile = compile_to_hardware(seqs, 'EchoCR/EchoCR',
        axis_descriptor=[
            delay_descriptor(np.concatenate((lengths, lengths))),
            cal_descriptor((controlQ, targetQ), calRepeats)
        ])

    if showPlot:
        plot_pulse_files(metafile)

    return metafile


def EchoCRPhase(controlQ,
                targetQ,
                phases,
                riseFall=40e-9,
                amp=1,
                length=100e-9,
                calRepeats=2,
                showPlot=False):
    """
	Variable phase CX experiment, with echo pulse sandwiched between two CR opposite-phase pulses.

	Parameters
	----------
	controlQ : logical channel for the control qubit (LogicalChannel)
	CRchan: logical channel for the cross-resonance pulse (LogicalChannel)
	phases : pulse phases of the CR pulse to sweep over (iterable)
	riseFall : rise/fall time of the CR pulse (s)
	amp : amplitude of the CR pulse
	length : duration of each of the two flat parts of the CR pulse (s)
	calRepeats : number of repetitions of readout calibrations for each 2-qubit state
	showPlot : whether to plot (boolean)
	"""
    seqs = [[Id(controlQ),
             echoCR(controlQ, targetQ, length=length, phase=ph, amp=amp, riseFall=riseFall),
             X90(targetQ)*Id(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for ph in phases] + \
           [[X(controlQ),
             echoCR(controlQ, targetQ, length=length, phase= ph, amp=amp, riseFall = riseFall),
             X90(targetQ)*X(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for ph in phases] + \
             create_cal_seqs((targetQ,controlQ), calRepeats, measChans=(targetQ,controlQ))

    axis_descriptor = [
        {
            'name': 'phase',
            'unit': 'radians',
            'points': list(phases)+list(phases),
            'partition': 1
        },
        cal_descriptor((controlQ, targetQ), calRepeats)
    ]

    metafile = compile_to_hardware(seqs, 'EchoCR/EchoCR',
        axis_descriptor=axis_descriptor)

    if showPlot:
        plot_pulse_files(metafile)

    return metafile


def EchoCRAmp(controlQ,
              targetQ,
              amps,
              riseFall=40e-9,
              length=50e-9,
              phase=0,
              calRepeats=2,
              showPlot=False):
    """
	Variable amplitude CX experiment, with echo pulse sandwiched between two CR opposite-phase pulses.

	Parameters
	----------
	controlQ : logical channel for the control qubit (LogicalChannel)
	targetQ: logical channel for the target qubit (LogicalChannel)
	amps : pulse amplitudes of the CR pulse to sweep over (iterable)
	riseFall : rise/fall time of the CR pulse (s)
	length : duration of each of the two flat parts of the CR pulse (s)
	phase : phase of the CR pulse (rad)
	calRepeats : number of repetitions of readout calibrations for each 2-qubit state
	showPlot : whether to plot (boolean)
	"""
    seqs = [[Id(controlQ),
             echoCR(controlQ, targetQ, length=length, phase=phase, riseFall=riseFall,amp=a),
             Id(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for a in amps] + \
           [[X(controlQ),
             echoCR(controlQ, targetQ, length=length, phase= phase, riseFall=riseFall,amp=a),
             X(controlQ),
             MEAS(targetQ)*MEAS(controlQ)] for a in amps] + \
           create_cal_seqs((targetQ,controlQ), calRepeats, measChans=(targetQ,controlQ))

    axis_descriptor = [
        {
            'name': 'amplitude',
            'unit': None,
            'points': list(amps)+list(amps),
            'partition': 1
        },
        cal_descriptor((controlQ, targetQ), calRepeats)
    ]

    metafile = compile_to_hardware(seqs, 'EchoCR/EchoCR',
        axis_descriptor=axis_descriptor)
    if showPlot:
        plot_pulse_files(metafile)

    return metafile
