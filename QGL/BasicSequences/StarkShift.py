from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..ChannelLibraries import EdgeFactory
from ..PulseSequencePlotter import plot_pulse_files
from .helpers import create_cal_seqs, delay_descriptor, cal_descriptor
import numpy as np
from collections.abc import Iterable
from itertools import product

def StarkSpectroscopy(qubit, amplitude, delay=200e-9, length=1e-6, showPlot=False, meas_chan=None):
    """Stark shift spectroscopy experiment. Applies a coherent displacement
        to the qubit readout cavity while doing pulsed spectroscopy.

        Args:
            qubit: Qubit channel to apply spectroscopy pulse to.

            amplitude: Measurement pulse amplitude(s)

            delay: Delay between end of spectroscopy pulse and start of MEAS(qubit). Default 200 ns.

            length: Total length of cavity displacement pulse. Default 1 us.

            meas_chan: Measurement channel. If None, assumes qubit's measurement channel.

        Returns:
            metafile: Path to compiled sequence metafile.
    """

    if not isinstance(amplitude, Iterable):
        amplitude = [amplitude]

    if meas_chan is None:
        meas_chan = qubit.measure_chan

    def stark_shift_pulse(amp):
        pump_pulse = Utheta(meas_chan, amp=amp, length=length)
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

def StarkEcho(qubit, amplitudes, delays, wait=200e-9, periods=4, showPlot=False, meas_chan=None):
    """Hahn echo sequence with a coherent displacement of the qubit measurement cavity.
        Used to measure photon-induced dephasing. This sequence can cause a lot of cache pressure
        so number of points may be limited.

        TODO: Use QGL intrinsics to reduce sequence and memory cache utilization.

        Args:
            qubit: Qubit channel for Hahn echo.

            amplitudes: Amplitude(s) of cavity displacement pulse.

            delays: Hahn echo delays - the t in 90-t-180-t-180.

            wait: Delay between end of cavity displacement pulse and start of MEAS(qubit).

            periods: Number of artificial oscillations.

            meas_chan: Measurement channel. If None, assumes qubit's measurement channel.


        Returns:
            metafile: Path to compiled sequence metafile.
    """

    if not isinstance(amplitudes, Iterable):
        amplitudes = [amplitudes]

    if not isinstance(delays, Iterable):
        delays = [delays]

    if meas_chan is None:
        meas_chan = qubit.measure_chan

    def echo_phase(n):
        return 2*np.pi*periods/len(delays)*n

    def echo_stark(n, amp, max_delay, meas_delay=200e-9):
        x_len = qubit.pulse_params["length"]
        max_len = 3*x_len + 2*max_delay + meas_delay
        echo_wait = max_len - (3*x_len + 2*delays[n])

        echo_seq = Id(qubit, echo_wait) + X90(qubit) + Id(qubit, delays[n]) + \
                        Y(qubit) + Id(qubit, delays[n]) + U90(qubit, echo_phase(n))

        meas_seq = Utheta(meas_chan, amp=amp, length=max_len)

        return echo_seq*meas_seq


    seqs = [[echo_stark(n, amp, np.max(delays)), Id(measurement, length=wait), MEAS(qubit)]
                for n, amp in product(range(len(delays)), amplitudes)]

    axis_descriptor = [delay_descriptor(delays)] * len(amplitudes)

    metafile = compile_to_hardware(seqs, 'StarkEcho/StarkEcho', axis_descriptor=axis_descriptor)

    if showPlot:
        plot_pulse_files(metafile)

    return metafile


def CavityPumpProbe(qubit, offsets, amplitude, length=1e-6, wait=1e-6, showPlot=False, meas_chan=None):
    """Time resolved cavity spectroscopy. Applies a coherent displacement to qubit
        readout cavity while sweeping qubit spectroscopy pulse delay. Useful to measure
        cavity kappa and cavity population.

        Args:
            qubit: Qubit channel for spectroscopy.

            offsets: Spectroscopy pulse offset relative to start of cavity displacement pulse.

            amplitude: Measurement pulse amplitude.

            length: Total length of cavity displacement pulse.

            wait: Delay between end of cavity displacement pulse and start of MEAS(qubit).

            meas_chan: Measurement channel. If None, assumes qubit's measurement channel.

        Returns:
            metafile: Path to compiled sequence metafile.
    """

    if not isinstance(offsets, Iterable):
        offsets = [offsets]

    if meas_chan is None:
        meas_chan = qubit.measure_chan

    def cavity_probe(offset):
        pump_pulse = Utheta(meas_chan, amp=amplitude, length=length)
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

def StarkRamsey(qubit, relax_delays, ramsey_delays, tppi_freq=4e6, meas_amp=0.05, wait_time=100e-9, meas_chan=None,
                prep_state=0, showPlot=False):
    """Measure photon decay in qubit cavity using Ramsey sequence. The sequence is:
        PREP - MEAS --- relax_delay ---- X90 --- ramsey_delay --- X90 --- wait_time --- MEAS
        You can use auspex.analysis.qubit_fits.PhotonNumberFit to fit this data.
        See: D.T. McClure, "Rapid Driven Reset of a Qubit Readout Resonator". PRA 5, 011001 (2016)

        Args:
            qubit: Qubit channel for spectroscopy.

            relax_delays: Delay to wait before start of ramsey sequence.

            ramsey_delays: Delays for Ramsey sequence.

            tppi_freq: Artificial detuning of Ramsey sequence.

            wait_time: Wait time between end of Ramsey and start of final measurement.

            meas_chan: Measurement channel. If None, assumes qubit's measurement channel.

            prep_state: Prepare qubit in ground or excited state.

        Returns:
            metafile: Path to compiled sequence metafile.
    """
    if meas_chan is None:
        meas_chan = qubit.measure_chan

    if prep_state == 1:
        prep = X(qubit)
    else:
        prep = Id(qubit)

    seqs = [[prep, MEAS(qubit, amp=meas_amp, dig_trig=False), Id(qubit, tw), X90(qubit), Id(qubit,tr),
                U90(qubit,phase = 2*np.pi*tppi_freq*tr), Id(qubit, wait_time), MEAS(qubit)]
                    for tw, tr in product(relax_delays, ramsey_delays)]

    seqs += create_cal_seqs((qubit,), 2)

    ## TODO: Fix this.
    #axis_descriptor = [delay_descriptor(ramsey_delays)]*len(relax_delays) + [cal_descriptor((qubit,), 2)]

    metafile = compile_to_hardware(seqs, 'StarkRamsey/StarkRamsey')

    if showPlot:
        plot_pulse_files(metafile)

    return metafile
