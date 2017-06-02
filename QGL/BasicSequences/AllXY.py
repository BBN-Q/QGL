from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..PulseSequencePlotter import plot_pulse_files
import QGL.PulseShapes
from .helpers import create_cal_seqs


def AllXY(q, showPlot=False):
    firstPulses = [Id(q)] + 2 * [X(q), Y(q)] + 3 * [X90(q), Y90(q)] + [
        X(q), Y(q), X90(q), X(q), Y90(q), Y(q), X(q), Y(q), X90(q), Y90(q)
    ]
    secondPulses = [Id(q), X(q), Y(q), Y(q), X(q), Id(q), Id(q), Y90(q),
                    X90(q), Y(q), X(q), Y90(q), X90(q), X(q), X90(q), Y(q),
                    Y90(q), Id(q), Id(q), X90(q), Y90(q)]
    seqs = []
    seqs += [[firstPulses[int(np.floor(ind / 2))],
              secondPulses[int(np.floor(ind / 2))], MEAS(q)]
             for ind in range(42)]

    metafile = compile_to_hardware(seqs, 'AllXY/AllXY')

    if showPlot:
        plot_pulse_files(metafile)

    return metafile