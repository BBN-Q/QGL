'''
All generic pulse shapes are defined here.
'''

import numpy as np

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
    #Return a single point at 0
    # return np.zeros(1, dtype=np.complex)
    numPts = np.round(length*samplingRate)
    return np.zeros(numPts, dtype=np.complex)


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
    QQuad = dragScaling*derivScale*xPts*IQuad
    return amp * (IQuad+1j*QQuad)
        
def gaussOn(amp=1, length=0, cutoff=2, samplingRate=1e9, **params):
    '''
    A half-gaussian pulse going from zero to full
    '''
    #Round to how many points we need
    numPts = np.round(length*samplingRate)
    xPts = np.linspace(-cutoff, 0, numPts)
    xStep = xPts[1] - xPts[0]
    return (amp * (np.exp(-0.5*(xPts**2)) - np.exp(-0.5*((xPts[0]-xStep)**2)))).astype(np.complex)

def gaussOff(amp=1, length=0, cutoff=2, samplingRate=1e9, **params):
    '''
    A half-gaussian pulse going from zero to full
    '''
    #Round to how many points we need
    numPts = np.round(length*samplingRate)
    xPts = np.linspace(0, cutoff, numPts)
    xStep = xPts[1] - xPts[0]
    return (amp * (np.exp(-0.5*(xPts**2)) - np.exp(-0.5*((xPts[-1]+xStep)**2)))).astype(np.complex)

def tanh(amp=1, length=0, cutoff=2, samplingRate=1e9, **params):
    '''
    A rounded square shape from the sum of two tanh shapes. 
    '''
    numPts = np.round(length*samplingRate)
    xPts = np.arange(numPts)-cutoff
    tanhPts = amp*np.tanh(xPts)
    return (0.5*(tanhPts + tanhPts[::-1])).astype(np.complex)
