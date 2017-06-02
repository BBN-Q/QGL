from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files

from itertools import chain
from numpy import pi


def SPAM(qubit, angleSweep, maxSpamBlocks=10, showPlot=False):
    """

	X-Y sequence (X-Y-X-Y)**n to determine quadrature angles or mixer correction.

	Parameters
	----------
	qubit : logical channel to implement sequence (LogicalChannel) 
	angleSweep : angle shift to sweep over
	maxSpamBlocks : maximum number of XYXY block to do
	showPlot : whether to plot (boolean)

	Returns
	-------
	plotHandle : handle to plot window to prevent destruction
	"""

    def spam_seqs(angle):
        """ Helper function to create a list of sequences increasing SPAM blocks with a given angle. """
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