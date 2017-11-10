from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
from .helpers import create_cal_seqs, delay_descriptor, cal_descriptor


def HahnEcho(qubit, pulseSpacings, periods=0, calRepeats=2, showPlot=False):
    """
	A single pulse Hahn echo with variable phase of second pi/2 pulse.

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel)
	pulseSpacings : pulse spacings to sweep over; the t in 90-t-180-t-180 (iterable)
	periods: number of artificial oscillations
	calRepeats : how many times to repeat calibration scalings (default 2)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""
    seqs = []
    for k in range(len(pulseSpacings)):
        seqs.append([X90(qubit), Id(qubit, pulseSpacings[k]), Y(qubit), Id(qubit,pulseSpacings[k]), \
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


def CPMG(qubit, numPulses, pulseSpacing, calRepeats=2, showPlot=False):
    """
	CPMG pulse train with fixed pulse spacing. Note this pulse spacing is centre to centre,
	i.e. it accounts for the pulse width

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel)
	numPulses : number of 180 pulses; should be even (iterable)
	pulseSpacing : spacing between the 180's (seconds)
	calRepeats : how many times to repeat calibration scalings (default 2)
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
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
            delay_descriptor(pulseSpacing * numPulses),
            cal_descriptor((qubit,), calRepeats)
        ])

    if showPlot:
        plot_pulse_files(metafile)
    return metafile
