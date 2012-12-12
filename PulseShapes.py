'''
All generic pulse shapes are defined here.
'''

import numpy as np

def gaussian(params, AWGFreq):
    '''
    A simple gaussian shaped pulse. 
    cutoff is how many sigma the pulse goes out
    '''
    #Round to how many points we need
    numPts = np.round(params.pulseLength*AWGFreq)
    xPts = np.linspace(-params.cutoff, params.cutoff, numPts)
    xStep = xPts[1] - xPts[0]
    return np.exp(-0.5*(xPts**2)) - np.exp(-0.5*((xPts[-1]+xStep)**2))
        
def square(params, AWGFreq):
    '''
    A simple rectangular shaped pulse. 
    '''
    #Round to how many points we need
    numPts = round(params.pulseLength*AWGFreq)
    return np.ones(numPts)
        
def delay(params, AWGFreq):
    '''
    A delay between pulses.
    '''
    #Return a single point at 0
    # return np.zeros(1, dtype=np.complex)
    numPts = np.round(params.pulseLength*AWGFreq)
    return np.zeros(numPts, dtype=np.complex)


def drag(params, AWGFreq):
    '''
    A gaussian pulse with a drag correction on the quadrature channel.
    '''
    #Create the gaussian along x and the derivative along y
    numPts = np.round(params.pulseLength*AWGFreq)
    xPts = np.linspace(-params.cutoff, params.cutoff, numPts)
    xStep = xPts[1] - xPts[0]
    IQuad = np.exp(-0.5*(xPts**2)) - np.exp(-0.5*((xPts[0]-xStep)**2))
    #The derivative needs to be scaled in terms of AWG points from the normalized xPts units.
    #The pulse length is 2*cutoff xPts
    derivScale = 1/(params.pulseLength/2/params.cutoff*AWGFreq)
    QQuad = params.dragScaling*derivScale*xPts*IQuad
    return IQuad+1j*QQuad
        
def gaussOn(params, AWGFreq):
    '''
    A half-gaussian pulse going from zero to full
    '''
    #Round to how many points we need
    numPts = np.round(params.pulseLength*AWGFreq)
    xPts = np.linspace(-params.cutoff, 0, numPts)
    xStep = xPts[1] - xPts[0]
    return np.exp(-0.5*(xPts**2)) - np.exp(-0.5*((xPts[0]-xStep)**2))

def gaussOff(params, AWGFreq):
    '''
    A half-gaussian pulse going from zero to full
    '''
    #Round to how many points we need
    numPts = np.round(params.pulseLength*AWGFreq)
    xPts = np.linspace(0, params.cutoff, numPts)
    xStep = xPts[1] - xPts[0]
    return np.exp(-0.5*(xPts**2)) - np.exp(-0.5*((xPts[-1]+xStep)**2))