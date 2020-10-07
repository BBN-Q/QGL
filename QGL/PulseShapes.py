'''
All generic pulse shapes are defined here.
'''

import numpy as np
from math import pi, sin, cos, acos, sqrt


def gaussian(amp=1, length=0, cutoff=2, sampling_rate=1e9, **params):
    '''
    A simple gaussian shaped pulse.
    cutoff is how many sigma the pulse goes out
    '''
    if length == 0:
        raise ValueError("gaussian() got 0 length")
    if sampling_rate == 0:
        raise ValueError("gaussian() got 0 sampling_rate")
    #Round to how many points we need
    numPts = int(np.round(length * sampling_rate))
    xPts = np.linspace(-cutoff, cutoff, numPts)
    xStep = xPts[1] - xPts[0]
    nextPoint = np.exp(-0.5 * ((xPts[0] - xStep)**2))
    #Rescale so that the maximum equals amp
    amp = (amp / (1 - nextPoint))
    return (amp * (np.exp(-0.5 * (xPts**2)) - np.exp(-0.5 * (
        (xPts[-1] + xStep)**2)))).astype(np.complex)


def delay(length=0, sampling_rate=1e9, **params):
    '''
    A delay between pulses.
    '''
    return constant(0, length, sampling_rate)


def constant(amp=1, length=0, sampling_rate=1e9, **params):
    '''
    A constant section.
    '''
    numPts = int(np.round(length * sampling_rate))
    return amp * np.ones(numPts, dtype=np.complex)

# square is deprecated but alias square to constant
square = constant

def drag(amp=1,
         length=0,
         cutoff=2,
         drag_scaling=0.5,
         sampling_rate=1e9,
         **params):
    '''
    A gaussian pulse with a drag correction on the quadrature channel.
    '''
    #Create the gaussian along x and the derivative along y
    numPts = int(np.round(length * sampling_rate))
    xPts = np.linspace(-cutoff, cutoff, numPts)
    xStep = xPts[1] - xPts[0]
    nextPoint = np.exp(-0.5 * ((xPts[0] - xStep)**2))
    #Rescale so that the maximum in IQuad equals amp
    amp = (amp / (1 - nextPoint))
    IQuad = np.exp(-0.5 * (xPts**2)) - np.exp(-0.5 * ((xPts[0] - xStep)**2))
    #The derivative needs to be scaled in terms of AWG points from the normalized xPts units.
    #The pulse length is 2*cutoff xPts
    derivScale = 1 / (length / 2 / cutoff * sampling_rate)
    QQuad = drag_scaling * derivScale * xPts * np.exp(-0.5 * (xPts**2))
    return amp * (IQuad + 1j * QQuad)

def alt_gaussian(pulse_length, sampling_rate, cutoff=2):
    num_points = int(round(pulse_length*sampling_rate))

    # Cut off the gaussian at some number of std devs away
    x_points = np.linspace(-cutoff, cutoff, num_points)

    # Compute the gaussian itself
    gaussian_points = np.exp(-0.5 * (x_points**2))

    # Shift it down so that the cutoff point becomes zero
    gaussian_points -= np.min(gaussian_points)

    # Scale to an amp of 1
    gaussian_points /= np.max(gaussian_points)

    # Pad with zeroes to match the vector length and return
    return gaussian_points

def gng(amp=1,
        length=0,
        sampling_rate=1e9,
        secondary_gaussian_width=0,
        cutoff=2,
        notch_freq=0,
        **params):
    # Create the original gaussian pulse and the second term which will be upconverted
    print('Making gng with amp={:.2f} length={:f} sampling_rate={:f} sgw={:f} cutoff={:d} notch_freq={:f}'.format(amp,length,sampling_rate,secondary_gaussian_width, cutoff, notch_freq))
    g = alt_gaussian(length, sampling_rate, cutoff)

    if notch_freq == 0 or secondary_gaussian_width == 0:
        return g

    secondary_g = alt_gaussian(secondary_gaussian_width, sampling_rate, cutoff)

    # Create the upconverion
    secondary_time_points = np.linspace(0, secondary_gaussian_width, int(round(secondary_gaussian_width*sampling_rate)))
    shift = np.exp(2*np.pi*notch_freq*secondary_time_points)

    # Upconvert the second gaussian and calculate the spectrum of both
    shifted_fft = np.fft.fft(list(secondary_g*shift) + [0]*(len(g) - len(secondary_g)))
    fft = np.fft.fft(list(g) + [0]*(len(secondary_g) - len(g)))

    # Scale and rotate the shifted gaussian so that at the notch, the frequency component disappears
    idx_notch = np.argmin(np.abs(np.fft.fftfreq(len(g), d=1/sampling_rate) - notch_freq))
    scale_factor = (np.abs(fft[idx_notch])/np.abs(shifted_fft[idx_notch]))
    phase_factor = np.exp(1j*(np.angle(fft[idx_notch]) - np.angle(shifted_fft[idx_notch]) + np.pi))

    # Calculate the new spectrum with the desired tone cancelled
    fft += shifted_fft*scale_factor*phase_factor

    # Return the time samples scaled to an absolute magnitude of 1
    ifft = np.fft.ifft(fft)
    return ifft * amp / np.max(np.abs(ifft))

def notched_gaussian(amp=1,
                     length=0,
                     sampling_rate=1e9,
                     cutoff=2,
                     notch_freq=0,
                     notch_width=0,
                     **params):
    notch_freqs = [notch_freq]
    notch_widths = [notch_width]
    g = alt_gaussian(length, sampling_rate, cutoff)
    fft = np.fft.fft(g)
    fft_freqs = np.fft.fftfreq(fft.size, d=1/sampling_rate)

    # Find the index of the frequency closest to the edges of the desired notches
    idx_notch_starts = [np.argmin(np.abs(fft_freqs - c - (notch_widths[idx]/2))) for idx,c in enumerate(notch_freqs)]
    idx_notch_stops = [np.argmin(np.abs(fft_freqs - c + (notch_widths[idx]/2))) for idx,c in enumerate(notch_freqs)]

#     print(str(pulse_length) + ' ' + str(idx_notch_starts[0]) + str(idx_notch_stops[0]))

    for i in range(len(notch_freqs)):
        zeros = np.zeros(abs(idx_notch_stops[i] - idx_notch_starts[i]), dtype=fft.dtype)
        if(idx_notch_stops[i] > idx_notch_starts[i]):
            fft[idx_notch_starts[i]:idx_notch_stops[i]] = zeros
        else:
            fft[idx_notch_stops[i]:idx_notch_starts[i]] = zeros

    ifft = np.fft.ifft(fft)
    return ifft * amp / np.max(np.abs(ifft))

def gaussOn(amp=1, length=0, cutoff=2, sampling_rate=1e9, **params):
    '''
    A half-gaussian pulse going from zero to full
    '''
    #Round to how many points we need
    numPts = int(np.round(length * sampling_rate))
    xPts = np.linspace(-cutoff, 0, numPts)
    #Pull the edge down to zero so there is no big step
    #i.e. find the shift such that the next point in the pulse would be zero
    xStep = xPts[1] - xPts[0]
    nextPoint = np.exp(-0.5 * ((xPts[0] - xStep)**2))
    #Rescale so that it still goes to amp
    amp = (amp / (1 - nextPoint))
    return (amp * (np.exp(-0.5 * (xPts**2)) - nextPoint)).astype(np.complex)

def cosOn(amp=1, length=0, sampling_rate=1e9, **params):
    '''
    A half-cosine pulse going from zero to full
    '''
    numPts = int(np.round(length * sampling_rate))
    times = np.linspace(0, length, numPts)
    vals = amp*(1 - np.cos((times * np.pi / length)))/2
    return vals.astype(np.complex)

def cosOff(amp=1, length=0, sampling_rate=1e9, **params):
    '''
    A half-cosine pulse going from full to zero
    '''
    numPts = int(np.round(length * sampling_rate))
    times = np.linspace(0, length, numPts)
    vals = amp*(1 - np.cos((times * np.pi / length) + np.pi))/2
    return vals.astype(np.complex)


def gaussOff(amp=1, length=0, cutoff=2, sampling_rate=1e9, **params):
    '''
    A half-gaussian pulse going from full to zero
    '''
    #Round to how many points we need
    numPts = int(np.round(length * sampling_rate))
    xPts = np.linspace(0, cutoff, numPts)
    #Pull the edge down to zero so there is no big step
    #i.e. find the shift such that the next point in the pulse would be zero
    xStep = xPts[1] - xPts[0]
    nextPoint = np.exp(-0.5 * ((xPts[-1] + xStep)**2))
    #Rescale so that it still goes to amp
    amp = (amp / (1 - nextPoint))
    return (amp * (np.exp(-0.5 * (xPts**2)) - nextPoint)).astype(np.complex)


def dragGaussOn(amp=1,
                length=0,
                cutoff=2,
                drag_scaling=0.5,
                sampling_rate=1e9,
                **params):
    '''
    A half-gaussian pulse with drag correction going from zero to full
    '''
    numPts = int(np.round(length * sampling_rate))
    xPts = np.linspace(-cutoff, 0, numPts)
    xStep = xPts[1] - xPts[0]
    IQuad = np.exp(-0.5 * (xPts**2)) - np.exp(-0.5 * ((xPts[0] - xStep)**2))
    derivScale = 1 / (length / 2 / cutoff * sampling_rate)
    QQuad = drag_scaling * derivScale * xPts * IQuad
    return amp * (IQuad + 1j * QQuad)


def dragGaussOff(amp=1,
                 length=0,
                 cutoff=2,
                 drag_scaling=0.5,
                 sampling_rate=1e9,
                 **params):
    '''
    A half-gaussian pulse with drag correction going from full to zero
    '''
    numPts = int(np.round(length * sampling_rate))
    xPts = np.linspace(0, cutoff, numPts)
    xStep = xPts[1] - xPts[0]
    IQuad = np.exp(-0.5 * (xPts**2)) - np.exp(-0.5 * ((xPts[-1] + xStep)**2))
    derivScale = 1 / (length / 2 / cutoff * sampling_rate)
    QQuad = drag_scaling * derivScale * xPts * IQuad
    return amp * (IQuad + 1j * QQuad)


def tanh(amp=1, length=0, sigma=0, cutoff=2, sampling_rate=1e9, **params):
    '''
    A rounded constant shape from the sum of two tanh shapes.
    '''
    if length == 0.0:
        return np.empty(shape=(0,)).astype(np.complex)
    else:
        numPts = int(np.round(length * sampling_rate))
        xPts = np.linspace(-length / 2, length / 2, numPts)
        x1 = -length / 2 + cutoff * sigma
        x2 = +length / 2 - cutoff * sigma
        assert x1 < 0 and x2 > 0, (f'Pulse length must be greater than'
                                   f'2 * cutoff (={cutoff}) * sigma '
                                   f'(={sigma}s).  Consider '
                                   f'using a Gaussian pulse instead.')
        return amp * 0.5 * (np.tanh((xPts - x1) / sigma) + np.tanh(
            (x2 - xPts) / sigma)).astype(np.complex)


def exp_decay(amp=1, length=0, sigma=0, sampling_rate=1e9, steady_state=0.4, **params):
    """
    An exponentially decaying pulse to try and populate the cavity as quickly as possible.
    But then don't overdrive it.
    """
    numPts = int(np.round(length * sampling_rate))
    timePts = (1.0 / sampling_rate) * np.arange(numPts)
    return amp * ((1-steady_state) * np.exp(-timePts / sigma) + steady_state).astype(np.complex)

# def composite(amp=1, resolution=1e6, band_start=0, band_stop=0.5e9, sampling_rate=1e9, **params):
#     """

#     """
#     numPts = int((1/resolution)*sampling_rate)
#     freqs = np.zeros(numPts, dtype=np.complex)
#     idx_start = int(round(numPts * band_start / sampling_rate))
#     idx_stop = int(round(numPts * band_stop / sampling_rate))
#     freqs[idx_start : idx_stop] = 1
#     samples = np.ifft(



def CLEAR(amp=1, length=0, sigma=0, sampling_rate=1e9, **params):
    """
    Pulse shape to quickly deplete the cavity at the end of a measurement.
    measPulse followed by 2 steps of length step_length and amplitudes amp1, amp2.
    """
    if 'amp1' not in params:
        params['amp1'] = 0
    if 'amp2' not in params:
        params['amp2'] = 0
    if 'step_length' not in params:
        params['step_length'] = 100e-9
    timePts = (1.0 / sampling_rate) * np.arange(np.round((length-2*params['step_length']) * sampling_rate))
    flat_step = amp * (0.6 * np.exp(-timePts / sigma) + 0.4).astype(np.complex)
    numPts_clear_step = int(np.round(params['step_length'] * sampling_rate))
    clear_step_one = amp * params['amp1'] * np.ones(numPts_clear_step, dtype=np.complex)
    clear_step_two = amp * params['amp2'] * np.ones(numPts_clear_step, dtype=np.complex)
    return np.append(flat_step, [clear_step_one, clear_step_two])

def autodyne(frequency=10e6, baseShape=constant, **params):
    '''
    A pulse with modulation at a particular frequency baked in.
    '''
    if isinstance(baseShape,str):
        shape = globals()[baseShape](**params)
    else:
        shape = baseShape(**params)
    # Apply the autodyne frequency
    timePts = np.linspace(0, params['length'], len(shape))
    shape *= np.exp(-1j * 2 * np.pi * frequency * timePts)
    return shape


def arb_axis_drag(nutFreq=10e6,
                  rotAngle=0,
                  polarAngle=0,
                  aziAngle=0,
                  length=0,
                  drag_scaling=0.5,
                  sampling_rate=1e9,
                  **params):
    """
    Single-qubit arbitrary axis pulse implemented with phase ramping and frame change.
    For now we assume gaussian shape.

    Parameters
    ----------
    nutFreq: effective nutation frequency per unit of drive amplitude (Hz)
    rotAngle : effective rotation rotAngle (radians)
    polarAngle : polar angle of rotation axis (radians)
    aziAngle : azimuthal (radians)
    """

    if length > 0:
        #Start from a gaussian shaped pulse
        gaussPulse = gaussian(amp=1,
                              length=length,
                              sampling_rate=sampling_rate,
                              **params)

        #Scale to achieve to the desired rotation
        calScale = (rotAngle / 2 / pi) * sampling_rate / sum(gaussPulse)

        #Calculate the phase ramp steps to achieve the desired Z component to the rotation axis
        phaseSteps = -2 * pi * cos(
            polarAngle) * calScale * gaussPulse / sampling_rate

        #Calculate Z DRAG correction to phase steps
        #beta is a conversion between XY drag scaling and Z drag scaling
        beta = drag_scaling / sampling_rate
        instantaneousDetuning = beta * (2 * pi * calScale * sin(polarAngle) *
                                        gaussPulse)**2
        phaseSteps = phaseSteps + instantaneousDetuning * (1.0 / sampling_rate)
        #center phase ramp around the middle of the pulse time steps
        phaseRamp = np.cumsum(phaseSteps) - phaseSteps / 2

        frameChange = sum(phaseSteps)

        shape = (1.0 / nutFreq) * sin(polarAngle) * calScale * gaussPulse * \
                np.exp(1j * phaseRamp)

    elif abs(polarAngle) < 1e-10:
        #Otherwise assume we have a zero-length Z rotation
        frameChange = -rotAngle
        shape = np.array([], dtype=np.complex128)
    else:
        raise ValueError(
            'Non-zero transverse rotation with zero-length pulse.')

    return shape
