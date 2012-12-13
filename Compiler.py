'''
functions for compiling lists of pulses/pulseBlocks down to the hardware level.
'''

import numpy as np

def TekChannels():
    '''
    The set of empty channels for a Tektronix AWG
    '''
    return {'ch1':[], 'ch2':[], 'ch3':[], 'ch4':[], 'ch1m1':[], 'ch1m2':[], 'ch2m1':[], 'ch2m2':[], 'ch3m1':[], 'ch3m2':[] , 'ch4m1':[], 'ch4m2':[]}

def APSChannels():
    '''
    The set of empty channels for a BBN APS.
    '''
    return {chanStr:{'LLs':[], 'WFLibrary':{0:np.zeros(1)}} for chanStr in  ['ch1','ch2','ch3','ch4']}

class LLElement(object):
    def __init__(self):
        self.key = None
        self.length = 0
        self.repeat = 1
        self.isTimeAmp = False
        self.hasTrigger = False
        self.triggerDelay = 0
        self.linkListRepeat = 0

def create_padding_LL():
    tmpLL = LLElement()
    tmpLL.isTimeAmp = True
    tmpLL.key = 'TAZKey'
    return tmpLL



def compile_sequence(seq, wfLib = {} ):
    '''
    Main function to convert sequences to miniLL's and waveform libraries.
    '''
    # normalize sequence to PulseBlocks
    seq = [p.promote() for p in seq]

    #Find the set of logical channels used here and initialize them
    channels = set([])
    for step in seq:
        channels |= set(step.pulses.keys())

    logicalLLs = {}        
    for chan in channels:
        logicalLLs[chan] = []
        if chan not in wfLib:
            wfLib[chan] = {"TAZKey":  np.zeros(1, dtype=np.complex)}


    for block in seq:
        #Align the block 
        blockLength = block.maxPts
        for chan in channels:
            if chan in block.pulses.keys():
                # add aligned LL entry
            else:
                # add identity

        