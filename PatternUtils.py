import numpy as np
import create_padding_LL from Compiler

def delay(linkList, delay):
    pass

def modulate(linkList, SSBFreq):
    pass

def align(linkList, mode, length):
    for miniLL in linkList:
        miniLL_length = sum([entry.lengthU*entry.repeat for entry in miniLL])
        paddingEntry = create_padding_LL(length - miniLL_length)
        if mode == 'left':
            miniLL.append(paddingEntry)
        elif mode == 'right':
            miniLL.insert(0, paddingEntry)
        else:
            raise NameError("Unknown aligment mode")

def correctMixer(wfLib, T):
    pass

def split_multiple_triggers():
	'''
	Split entries with multiple triggers into two entries.
	'''
	pass