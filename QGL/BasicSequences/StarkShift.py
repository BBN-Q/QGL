from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..ChannelLibraries import EdgeFactory
from ..PulseSequencePlotter import plot_pulse_files
from .helpers import create_cal_seqs, delay_descriptor, cal_descriptor
import sqlalchemy
import numpy as np
from collections.abc import Iterable
from itertools import product

from typing import Iterable as IterableType
from typing import Union

def StarkSpectroscopy(qubit: Channels.LogicalChannel,
                      measurement: sqlalchemy.ext.declarative.api.DeclarativeMeta, 
                      amplitude: float, 
                      delay: Union[int, float] = 200e-9, 
                      length: Union[int, float] = 1e-6, 
                      showPlot: bool = False) -> str:
    """
    Stark shift spectroscopy experiment. Applies a coherent displacement
    to the qubit readout cavity while doing pulsed spectroscopy.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel for the control qubit
    measurement : bbndb.qgl.Measurement
        Measurement channel to apply displacement pulse to
    amplitude : float
        Amplitude of the measurement pulse. Valid range: [0.0, 1.0].
    delay : float, optional
        Delay between end of spectroscopy pulse and start of MEAS (seconds)
    lengths : int/float, optional
        Total length of cavity displacement pulse (seconds).  4 ns minimum.
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files.

    Examples
    --------
    >>> mf = StarkSpectroscopy(q1, q1.measure_chan, np.linspace(0.6, 0.8, 51));
    Compiled 51 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    if not isinstance(amplitude, Iterable):
        amplitude = [amplitude]

    def stark_shift_pulse(amp):
        pump_pulse = Utheta(measurement, amp=amp, length=length)
        l1 = length - delay - qubit.pulse_params["length"] - delay
        spec_pulse = Id(qubit, length=l1)+X(qubit)+Id(qubit,length=delay)
        return spec_pulse*pump_pulse

    seqs = [[stark_shift_pulse(a), MEAS(qubit)] for a in amplitude]
    axis_descriptor = [{
        'name': 'Stark Shift Amplitude',
        'unit': None,
        'points': list(amplitude),
        'partition': 1
    }]
    metafile = compile_to_hardware(seqs, 'StarkSpec/StarkSpec', axis_descriptor=axis_descriptor)

    if showPlot:
        plot_pulse_files(metafile)

    return metafile

def StarkEcho(qubit: Channels.LogicalChannel,
              measurement: sqlalchemy.ext.declarative.api.DeclarativeMeta, 
              amplitudes: IterableType[Union[int,float]], 
              delays: IterableType[Union[int,float]],
              wait: Union[int,float] = 200e-9, 
              periods: int = 4, 
              showPlot: bool = False) -> str:
    """
    Hahn echo sequence with a coherent displacement of the qubit measurement
    cavity. Used to measure photon-induced dephasing. This sequence can cause
    a lot of cache pressure so number of points may be limited.

    TODO: Use QGL intrinsics to reduce sequence and memory cache utilization.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel for the Hahn echo
    measurement : bbndb.qgl.Measurement
        Measurement channel of the qubit
    amplitude : int/float iterable
        Amplitude(s) of cavity displacement pulse. Valid range: [0.0, 1.0].
    delays : int/float iterable
        Delay between end of spectroscopy pulse and start of MEAS (seconds)
    wait : int/float, optional
        Hahn echo delays - the t in 90-t-180-t-180 (seconds)
        (seconds).  4 ns minimum.
    periods : int, optional
        Number of artificial oscillations
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files.

    Examples
    --------
    >>> mf = StarkEcho(q1, q1.measure_chan,
                       np.linspace(0.6, 0.8, 10),
                       np.linspace(20.0e-9, 200.02e-6, 10));
    Compiled 210 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    if not isinstance(amplitudes, Iterable):
        amplitudes = [amplitudes]

    if not isinstance(delays, Iterable):
        delays = [delays]

    def echo_phase(n):
        return 2*np.pi*periods/len(delays)*n

    def echo_stark(n, amp, max_delay, meas_delay=200e-9):
        x_len = qubit.pulse_params["length"]
        max_len = 3*x_len + 2*max_delay + meas_delay
        echo_wait = max_len - (3*x_len + 2*delays[n])

        echo_seq = Id(qubit, echo_wait) + X90(qubit) + Id(qubit, delays[n]) + \
                        Y(qubit) + Id(qubit, delays[n]) + U90(qubit, echo_phase(n))

        meas_seq = Utheta(measurement, amp=amp, length=max_len)

        return echo_seq*meas_seq


    seqs = [[echo_stark(n, amp, np.max(delays)), Id(measurement, length=wait), MEAS(qubit)]
                for n, amp in product(range(len(delays)), amplitudes)]

    axis_descriptor = [delay_descriptor(delays)] * len(amplitudes)

    metafile = compile_to_hardware(seqs, 'StarkEcho/StarkEcho', axis_descriptor=axis_descriptor)

    if showPlot:
        plot_pulse_files(metafile)

    return metafile


def CavityPumpProbe(qubit: Channels.LogicalChannel, 
                    measurement: sqlalchemy.ext.declarative.api.DeclarativeMeta, 
                    offsets: IterableType[Union[int,float]], 
                    amplitude: IterableType[Union[int,float]],
                    length: Union[int,float] = 1e-6, 
                    wait: Union[int,float] = 2e-6, 
                    showPlot: bool = False) -> str:
    """
    Time resolved cavity spectroscopy. Applies a coherent displacement to qubit
    readout cavity while sweeping qubit spectroscopy pulse delay. Useful to
    measure cavity kappa and cavity population.

    Parameters
    ----------
    qubit : Channels.LogicalChannel
        Logical channel for the Hahn echo
    measurement : bbndb.qgl.Measurement
        Measurement channel of the qubit
    offsets : int/float iterable
        Spectroscopy pulse offset relative to start of cavity displacement
        pulse (seconds)
    amplitude : int/float iterable
        Amplitude(s) of cavity displacement pulse. Valid range: [0.0, 1.0].
    length : int/float, optional
        Total length of cavity displacement pulse (seconds)
    wait : int/float, optional
        Delay between end of cavity displacement pulse and start of measurement
        (seconds).  4 ns minimum.
    showPlot : boolean, optional
        Whether to plot

    Returns
    -------
    metafile : string
        Path to a json metafile with details about the sequences and paths
        to compiled machine files.

    Examples
    --------
    >>> mf = CavityPumpProbe(q1, q1.measure_chan,
                       np.linspace(20.0e-9, 200.02e-6, 10),
                       0.6);
    Compiled 210 sequences.
    >>> mf
    '/path/to/exp/exp-meta.json'
    """

    if not isinstance(offsets, Iterable):
        offsets = [offsets]

    def cavity_probe(offset):
        pump_pulse = Utheta(measurement, amp=amplitude, length=length)
        x_len = qubit.pulse_params["length"]
        if offset < -1*x_len:
            return [X(qubit), Id(qubit, length=(-x_len-offset)), pump_pulse, Id(qubit, length=wait)]
        elif offset < 0:
            total_len = length-offset
            pm = Id(measurement, length=offset)+pump_pulse
            pq = X(qubit)+Id(qubit, length=(total_len-x_len))
            return [pm*pq, Id(qubit, length=wait)]
        elif offset < length:
            pq = Id(qubit, length=offset)+X(qubit)+Id(qubit, length=(length-offset-x_len))
            return [pump_pulse*pq, Id(qubit, length=wait)]
        elif offset >= length:
            assert offset < (length+wait), f"Wait time {wait} is too short!"
            wait_len = wait - (offset-length+x_len)
            return [pump_pulse, Id(qubit, length=(offset-length)), X(qubit), Id(qubit, length=wait_len)]

    seqs = [[cavity_probe(off), MEAS(qubit)] for off in offsets]
    axis_descriptor = [delay_descriptor(offsets)]
    metafile = compile_to_hardware(seqs, 'CavityPumpProbe/CavityPumpProbe', axis_descriptor=axis_descriptor)

    if showPlot:
        plot_pulse_files(metafile)

    return metafile
