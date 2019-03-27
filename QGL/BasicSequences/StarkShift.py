from ..PulsePrimitives import *
from ..Compiler import compile_to_hardware
from ..ChannelLibraries import EdgeFactory
from ..PulseSequencePlotter import plot_pulse_files
from .helpers import create_cal_seqs, delay_descriptor, cal_descriptor
import numpy as np
from collections.abc import Iterable

def StarkShiftSpectroscopy(qubit, measurement, amplitude, 
							delay=200e-9, length=1e-6, showPlot=False):
	"""Stark shift spectroscopy experiment. Applies a coherent displacement
		to the qubit readout cavity while doing pulsed spectroscopy.

		Args:
			qubit: Qubit channel to apply spectroscopy pulse to.

			measurement: Measurement channel to apply displacement pulse to.

			amplitude: Measurement pulse amplitude(s)

			delay: Delay between end of spectroscopy pulse and start of MEAS(qubit).

			length: Total length of cavity displacement pulse.

		Returns:
			metafile: Path to compiled sequence metafile.
	"""

	if not isinstance(amplitude, Iterable):
		amplitude = [amplitude]

	def stark_shift_pulse(amp):
		pump_pulse = Utheta(measurement, amp=amp, length=length)
	    l1 = length - delay - qubit.pulse_params["length"] - delay
	    spec_pulse = Id(qubit, length=l1)+X(qubit)+Id(qubit,length=delay)
	    return spec_pulse*pump_pulse

    seqs = [[stark_shift_pulse(a), MEAS(qubit)] for a in amplitude]
    axis_descriptor = {
        'name': 'Stark Shift Amplitude',
        'unit': None,
        'points': amplitude,
        'partition': 1
    }
    metafile = compile_to_hardware(seqs, 'StarkSpec/StarkSpec', axis_descriptor=axis_descriptor)

    if showPlot:
    	plot_pulse_files(metafile)

    return metafile

def StarkEcho(qubit, measurement, amplitudes, delays, 
							wait=200e-9, periods=4, showPlot=False):
	"""Hahn echo sequence with a coherent displacement of the qubit measurement cavity.
		Used to measure photon-induced dephasing. This sequence can cause a lot of cache pressure
		so number of points may be limited.

		Args:
			qubit: Qubit channel for Hahn echo.

			measurement: Measurement channel of qubit.
			
			amplitudes: Amplitude(s) of cavity displacement pulse. 

			delays: Hahn echo delays - the t in 90-t-180-t-180. 

			wait: Delay between end of cavity displacement pulse and start of MEAS(qubit).

			periods: Number of artificial oscillations.

		Returns:
			metafile: Path to compiled sequence metafile.
	"""

	if not isinstance(amplitudes, Iterable):
		amplitudes = [amplitudes]

	if not isinstance(delays, Iterable):
		delays = [delays]

	for k in range(len(pulseSpacings)):
        seqs.append([X90(qubit), Id(qubit, pulseSpacings[k]), Y(qubit), Id(qubit,pulseSpacings[k]), \
         U90(qubit,phase=2*pi*periods/len(pulseSpacings)*k), MEAS(qubit)])

    def echo_phase(n):
    	return 2*np.pi*periods/len(delays)*n

	def echo_stark(n, amp, max_delay, meas_delay=200e-9):
		x_len = qubit.pulse_params["length"]
		max_len = 3*x_len + 2*max_delay + meas_delay
		echo_wait = total_len - (3*x_len + 2*delays[n]) 

		echo_seq = Id(qubit, echo_wait) + X90(qubit) + Id(qubit, delays[n]) + \
						Y(qubit) + Id(qubit, delays[n]) + U90(qubit, echo_phase(n))

		meas_seq = Utheta(measurement, amp=amp, length=max_len)

		return echo_seq*meas_seq


    seqs = [[echo_stark(n, amp, np.max(delays)), Id(measurement, length=wait), MEAS(qubit)] 
    			for n, amp in product(range(len(delays)), amplitudes)]

    axis_descriptor = delay_descriptor(delays) * len(amplitudes)

    metafile = compile_to_hardware(seqs, 'StarkEcho/StarkEcho', axis_descriptor=axis_descriptor)

    if showPlot:
    	plot_pulse_files(metafile)

    return metafile


def CavityPumpProbe(qubit, measurement, offsets, amplitude, 
							length=1e-6, wait=1e-6, showPlot=False):
	"""Time resolved cavity spectroscopy. Applies a coherent displacement to qubit
		readout cavity while sweeping qubit spectroscopy pulse delay. Useful to measure
		cavity kappa and cavity population.

		Args:
			qubit: Qubit channel for spectroscopy.

			measurement: Measurement channel of qubit.
			
			offsets: Spectroscopy pulse offset relative to start of cavity displacement pulse.

			amplitude: Measurement pulse amplitude.

			length: Total length of cavity displacement pulse.

			wait: Delay between end of cavity displacement pulse and start of MEAS(qubit).

		Returns:
			metafile: Path to compiled sequence metafile.
	"""

	if not isinstance(offsets, Iterable):
		offsets = [offsets]

	def cavity_probe(offset):
		pump_pulse = Utheta(measurement, amp=amp, length=length)
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
    axis_descriptor = delay_descriptor(offsets)
    metafile = compile_to_hardware(seqs, 'CavityPumpProbe/CavityPumpProbe', axis_descriptor=axis_descriptor)

    if showPlot:
    	plot_pulse_files(metafile)

    return metafile



