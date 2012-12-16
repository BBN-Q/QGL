import numpy as np
import Compiler

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