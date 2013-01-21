import numpy as np
import Compiler
from warnings import warn
from APSPattern import MIN_ENTRY_LENGTH
import pdb

def delay(linkList, delay, samplingRate):
    '''
    Delays a mini link list by the given delay. Postives delays
    shift right, negative delays shift left.
    '''
    # we assume that link lists are constructed with a padding element at the beginning that can be adjusted
    assert linkList[0][0].isZero, 'Delay error: link list does not start with a padding element'
    sampShift = int(round(delay * samplingRate))
    for miniLL in linkList:
        miniLL[0].length += sampShift
        assert miniLL[0].length > 0, 'Delay error: Negative length padding element after delay.'

def modulate(linkList, SSBFreq):
    return linkList

def align(linkList, mode, length):
    for miniLL in linkList:
        miniLL_length = sum([entry.length*entry.repeat for entry in miniLL])
        paddingEntry = Compiler.create_padding_LL(length - miniLL_length)
        if mode == 'left':
            miniLL.append(paddingEntry)
        elif mode == 'right':
            miniLL.insert(0, paddingEntry)
        else:
            raise NameError("Unknown aligment mode")

def correctMixer(wfLib, T):
    for k, v in wfLib.items():
        # To get the broadcast to work in numpy, need to do the multiplication one row at a time
        iqWF = np.vstack((np.real(v), np.imag(v)))
        wfLib[k] = T[0,:].dot(iqWF) + 1j*T[1,:].dot(iqWF)

def split_multiple_triggers():
	'''
	Split entries with multiple triggers into two entries.
	'''
	pass

def create_gate_seqs(linkList, gateBuffer=0, gateMinWidth=0, samplingRate=1.2e9):
    '''
    Helper function that takes a set of analog channel LL and gates each miniLL
    '''

    # convert times into samples
    gateBuffer = int(round(gateBuffer * samplingRate))
    gateMinWidth = int(round(gateMinWidth * samplingRate))
    
    #Initialize list of sequences to return
    gateSeqs = []

    # Time from end of previous LL entry that trigger needs to go
    # high to gate pulse
    startDelay = gateBuffer
    for miniLL in linkList:
        #Initialize a zero length padding sequence
        gateSeqs.append([Compiler.create_padding_LL(0)])
        # we need to pad the miniLL with an extra entry if the last entry is not a zero
        if not miniLL[-1].isZero:
            miniLL.append(Compiler.create_padding_LL(max(MIN_ENTRY_LENGTH, gateBuffer)))
        if miniLL[-1].length < gateBuffer:
            miniLL[-1].length = gateBuffer


        #Initialze the state low and the sample count and a carry for delays longer than the entry
        state = 0 # 0 = low, 1 = high
        sampct = 0
        carry = 0
        #Loop over miniLL entries
        for curEntry, nextEntry in zip(miniLL[:-1], miniLL[1:]):
            entryWidth = curEntry.length * curEntry.repeat - carry
            carry = 0
            assert entryWidth>0
            # If current state is low and next linkList is pulse, then
            # we go high in this entry.
            # If current state is high and next entry is TAZ then go low
            # in this one (but check the gateMinWidth)
            # Otherwise maintain current state by adding to last entry's length
            if state == 0 and not nextEntry.isZero:
                #Finish up the current entry
                gateSeqs[-1][-1].length += entryWidth - startDelay
                #Create a new entry with the high state
                gateSeqs[-1].append(Compiler.create_padding_LL(startDelay, high=True))
                state = 1
            elif state == 1 and nextEntry.isZero and (((nextEntry.length*nextEntry.repeat) > gateMinWidth) or nextEntry==miniLL[-1]):
                # Time from beginning of pulse LL entry that trigger needs to go
                # low to end gate pulse
                endDelay = np.fix(entryWidth + gateBuffer);
                if endDelay < 0:
                    endDelay = 0
                    warn("gatePulses warning: fixed buffer pulse to start of pulse")
                gateSeqs[-1][-1].length += endDelay
                if endDelay < entryWidth:
                    gateSeqs[-1].append(Compiler.create_padding_LL(endDelay-entryWidth))
                else:
                    gateSeqs[-1].append(Compiler.create_padding_LL(0))
                    carry = endDelay-entryWidth
                state = 0
            else:
                gateSeqs[-1][-1].length += entryWidth

        # end loop through miniLL
    # end loop through link lists
    return gateSeqs

def add_marker_pulse(LL, startPt, length):
    '''
    Helper function to add a marker pulse to a LL from a given startPt and with a given length
    '''
    pass


def add_digitizer_trigger(awgData, trigChan, delays):
    '''
    Add the digitizer trigger.  For now hardcoded but should be loaded from config file.
    '''
    awgName, awgChan = trigChan.physicalChannel.name.split('-')

    #Assume awgChan comes in form xmx e.g. 2m1
    analogChan = int(awgChan[0])
    IQKey = 'ch12' if analogChan < 3 else 'ch34'
    awgModel = trigChan.physicalChannel.AWG.model
    if awgModel == 'Tek5000':
        mrksPerChan = 2
    elif awgModel == 'BBNAPS':
        mrksPerChan = 1
    else:
        raise NameError('Unknown AWG model: {0}'.format(awgModel))

    mrkChan = np.mod(analogChan-1,2)*mrksPerChan + int(awgChan[-1])

    for LL, delay in zip(awgData[IQKey][LL], delays):
        delayPts = round(delay*trigChan.physicalChannel.AWG.samplingRate/APSPattern.ADDRESS_UNIT)
        add_marker_pulse(LL, delayPts)


def add_slave_trigger(LL):
    pass


    
        



