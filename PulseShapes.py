'''
All generic pulse shapes are defined here.
'''

import numpy as np
from math import pi, sin, cos, acos, sqrt

def gaussian(amp=1, length=0, cutoff=2, samplingRate=1e9, **params):
    '''
    A simple gaussian shaped pulse. 
    cutoff is how many sigma the pulse goes out
    '''
    #Round to how many points we need
    numPts = np.round(length*samplingRate)
    xPts = np.linspace(-cutoff, cutoff, numPts)
    xStep = xPts[1] - xPts[0]
    return (amp * (np.exp(-0.5*(xPts**2)) - np.exp(-0.5*((xPts[-1]+xStep)**2)))).astype(np.complex)
        
def square(amp=1, length=0, samplingRate=1e9, **params):
    '''
    A simple rectangular shaped pulse. 
    '''
    #Round to how many points we need
    numPts = np.round(length*samplingRate)
    return (amp * np.ones(numPts)).astype(np.complex)
        
def delay(length=0, samplingRate=1e9, **params):
    '''
    A delay between pulses.
    '''
    return constant(0, length, samplingRate)

def constant(amp=1, length=0, samplingRate=1e9, **params):
    '''
    A constant section.
    '''
    numPts = np.round(length*samplingRate)
    return amp*np.ones(numPts, dtype=np.complex)

def drag(amp=1, length=0, cutoff=2, dragScaling=0.5, samplingRate=1e9, **params):
    '''
    A gaussian pulse with a drag correction on the quadrature channel.
    '''
    #Create the gaussian along x and the derivative along y
    numPts = np.round(length*samplingRate)
    xPts = np.linspace(-cutoff, cutoff, numPts)
    xStep = xPts[1] - xPts[0]
    IQuad = np.exp(-0.5*(xPts**2)) - np.exp(-0.5*((xPts[0]-xStep)**2))
    #The derivative needs to be scaled in terms of AWG points from the normalized xPts units.
    #The pulse length is 2*cutoff xPts
    derivScale = 1/(length/2/cutoff*samplingRate)
    QQuad = dragScaling*derivScale*xPts * np.exp(-0.5*(xPts**2))
    return amp * (IQuad+1j*QQuad)
        
def gaussOn(amp=1, length=0, cutoff=2, samplingRate=1e9, **params):
    '''
    A half-gaussian pulse going from zero to full
    '''
    #Round to how many points we need
    numPts = np.round(length*samplingRate)
    xPts = np.linspace(-cutoff, 0, numPts)
    #Pull the edge down to zero so there is no big step
    #i.e. find the shift such that the next point in the pulse would be zero
    xStep = xPts[1] - xPts[0]
    nextPoint = np.exp(-0.5*((xPts[0]-xStep)**2))
    #Rescale so that it still goes to amp
    amp = (amp/(1-nextPoint))
    return (amp * (np.exp(-0.5*(xPts**2)) - nextPoint)).astype(np.complex)

def gaussOff(amp=1, length=0, cutoff=2, samplingRate=1e9, **params):
    '''
    A half-gaussian pulse going from full to zero
    '''
    #Round to how many points we need
    numPts = np.round(length*samplingRate)
    xPts = np.linspace(0, cutoff, numPts)
    #Pull the edge down to zero so there is no big step
    #i.e. find the shift such that the next point in the pulse would be zero
    xStep = xPts[1] - xPts[0]
    nextPoint = np.exp(-0.5*((xPts[-1]+xStep)**2))
    #Rescale so that it still goes to amp
    amp = (amp/(1-nextPoint))
    return (amp * (np.exp(-0.5*(xPts**2)) - nextPoint)).astype(np.complex)

def dragGaussOn(amp=1, length=0, cutoff=2, dragScaling=0.5, samplingRate=1e9, **params):
    '''
    A half-gaussian pulse with drag correction going from zero to full
    '''
    numPts = np.round(length*samplingRate)
    xPts = np.linspace(-cutoff, 0, numPts)
    xStep = xPts[1] - xPts[0]
    IQuad = np.exp(-0.5*(xPts**2)) - np.exp(-0.5*((xPts[0]-xStep)**2))
    derivScale = 1/(length/2/cutoff*samplingRate)
    QQuad = dragScaling*derivScale*xPts*IQuad
    return amp * (IQuad+1j*QQuad)

def dragGaussOff(amp=1, length=0, cutoff=2, dragScaling=0.5, samplingRate=1e9, **params):
    '''
    A half-gaussian pulse with drag correction going from full to zero
    '''
    numPts = np.round(length*samplingRate)
    xPts = np.linspace(0, cutoff, numPts)
    xStep = xPts[1] - xPts[0]
    IQuad = np.exp(-0.5*(xPts**2)) - np.exp(-0.5*((xPts[-1]+xStep)**2))
    derivScale = 1/(length/2/cutoff*samplingRate)
    QQuad = dragScaling*derivScale*xPts*IQuad
    return amp * (IQuad+1j*QQuad)

def tanh(amp=1, length=0, sigma=0, cutoff=2, samplingRate=1e9, **params):
    '''
    A rounded square shape from the sum of two tanh shapes. 
    '''
    numPts = np.round(length*samplingRate)
    xPts = np.linspace(-length/2, length/2, numPts)
    x1 = -length/2 + cutoff*sigma
    x2 = +length/2 - cutoff*sigma
    return amp*0.5*(np.tanh((xPts-x1)/sigma) + np.tanh((x2-xPts)/sigma)).astype(np.complex)

def measPulse(amp=1, length=0, sigma=0, samplingRate=1e9, **params):
    """
    An exponentially decaying pulse to try and populate the cavity as quickly as possible.
    But then don't overdrive it.
    """
    numPts = np.round(length*samplingRate)
    timePts = (1.0/samplingRate)*np.arange(numPts)
    return amp*(0.8*np.exp(-timePts/sigma) + 0.2).astype(np.complex)

def autodyne(frequency=10e6, baseShape=constant, **params):
    '''
    A pulse with modulation at a particular frequency baked in.
    '''
    shape = baseShape(**params)
    # Apply the autodyne frequency
    timePts = np.linspace(0, params['length'], len(shape))
    shape *= np.exp(-1j*2*np.pi*frequency*timePts)
    return shape

def arb_axis_drag(nutFreq=10e6, rotAngle=0, polarAngle=0, aziAngle=0, length=0, dragScaling=0.5, samplingRate=1e9, **params):
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
        gaussPulse = gaussian(amp=1, length=length, samplingRate=samplingRate, **params)

        #Scale to achieve to the desired rotation
        calScale = (rotAngle/2/pi)*samplingRate/sum(gaussPulse)

        #Calculate the phase ramp steps to achieve the desired Z component to the rotation axis
        phaseSteps = -2*pi*cos(polarAngle)*calScale*gaussPulse/samplingRate

        #Calculate Z DRAG correction to phase steps
        #beta is a conversion between XY drag scaling and Z drag scaling
        beta = dragScaling/samplingRate
        instantaneousDetuning = beta*(2*pi*calScale*sin(polarAngle)*gaussPulse)**2
        phaseSteps = phaseSteps + instantaneousDetuning*(1.0/samplingRate)
        #center phase ramp around the middle of the pulse time steps
        phaseRamp = np.cumsum(phaseSteps) - phaseSteps/2

        frameChange = sum(phaseSteps);

        shape = (1.0/nutFreq)*sin(polarAngle)*calScale*gaussPulse*np.exp(1j*phaseRamp)

    elif abs(polarAngle) < 1e-10:
        #Otherwise assume we have a zero-length Z rotation
        frameChange = -rotAngle;
        shape = np.array([], dtype=np.complex128)
    else:
        raise ValueError('Non-zero transverse rotation with zero-length pulse.')

    return shape
