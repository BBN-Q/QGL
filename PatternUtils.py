import numpy as np
import Compiler
from warnings import warn
from APSPattern import MIN_ENTRY_LENGTH

def delay(linkList, delay):
    return linkList

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
    return wfLib

def split_multiple_triggers():
	'''
	Split entries with multiple triggers into two entries.
	'''
	pass

def gatePulses(linkList, delay=0, gateBuffer=0, bufferReset=0, samplingRate=1.2e9, ch=1):
    # convert times into samples
    delay = round(delay * samplingRate)
    gateBuffer = round(gateBuffer * samplingRate)
    bufferReset = round(bufferReset * samplingRate)
    
    if ch == 1:
        markerStr = 'markerDelay1'
    elif ch == 2:
        markerStr = 'markerDelay2'
    else:
        raise NameError("Unknown marker channel {0}".format(ch))

    # Time from end of previous LL entry that trigger needs to go
    # high to gate pulse
    startDelay = round(gateBuffer - delay)
    for miniLL in linkList:
        # we need to pad the miniLL with an extra entry if the last entry is not a zero
        if not miniLL[-1].isZero:
            miniLL.append(Compiler.create_padding_LL(max(MIN_ENTRY_LENGTH, gateBuffer + delay)))

        state = 0 # 0 = low, 1 = high
        for ct in range(len(miniLL)-1):
            entryWidth = miniLL[ct].length * miniLL[ct].repeat
            # If current state is low and next linkList is pulse, then
            # we go high in this entry.
            # If current state is high and next entry is TAZ then go low
            # in this one (but check bufferReset)
            if state == 0 and not miniLL[ct+1].isZero:
                setattr(miniLL[ct], markerStr, entryWidth - startDelay)
                state = 1
            elif state == 1 and miniLL[ct+1].isZero and (miniLL[ct+1].length*miniLL[ct+1].repeat) > bufferReset:
                # Time from beginning of pulse LL entry that trigger needs to go
                # low to end gate pulse
                endDelay = np.fix(entryWidth + gateBuffer + delay);
                if endDelay < 0:
                    endDelay = 0
                    warn("gatePulses warning: fixed buffer pulse to start of pulse")
                setattr(miniLL[ct], markerStr, endDelay)
                state = 0
        # end loop through miniLL
    # end loop through link lists
    return linkList

        