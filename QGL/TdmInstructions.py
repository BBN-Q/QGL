'''
Copyright 2018 Raytheon BBN Technologies

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

# Instructions for the the prototype APS2-TDM hardware.

class CustomInstruction(object):

    def __init__(self, name, in_addr, out_addr, **kwargs):
        self.instruction = name
        self.in_addr = in_addr
        self.out_addr = out_addr
        self.kwargs = kwargs
        self.length = 0

    def promote(self, ptype):
        return self

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self == other


def MajorityVote(in_addr, out_addr, nmeas): # alternatively, append the loadcmpvram instruction when compiling (see STOREMEAS)
    return [LoadCmpVramInstruction('LOADCMPVRAM', 1, in_addr, 2**nmeas-1, True), CustomInstruction('MAJORITY', in_addr, out_addr)]

def MajorityMask(value):
    return [WriteAddrInstruction('INVALIDATE', None, 1, 0, 0x0, True), WriteAddrInstruction('WRITEADDR', None, 0, 0, value, True), LoadCmpVramInstruction('LOADCMPVRAM', 1, 0, 0xffff, True), CustomInstruction('MAJORITYMASK', 0, 1)]

def Decode(in_addr, out_addr, nmeas):
    return [LoadCmpVramInstruction('LOADCMPVRAM', 1, in_addr, 2**nmeas-1, True), CustomInstruction('TSM', in_addr, out_addr)]

def DecodeSetRounds(in_addr, out_addr, value):
    return [WriteAddrInstruction('INVALIDATE', None, 1, in_addr, 0x0, True), WriteAddrInstruction('WRITEADDR', None, 0, in_addr, value, True), LoadCmpVramInstruction('LOADCMPVRAM', 1, in_addr, 0xffff, True), CustomInstruction('TSM_SET_ROUNDS', in_addr, out_addr)]

# TODO: the rest of the CUSTOM instructions

class WriteAddrInstruction(object):

    def __init__(self, name, channel, modifier, addr, value, tdm, **kwargs):
        self.instruction = name
        self.channel = channel
        self.invalid = modifier
        self.addr = addr
        self.value = value
        self.kwargs = kwargs
        self.tdm = tdm
        self.length = 0

    def promote(self, ptype):
        return self

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self == other

def WriteAddr(addr, value, channel=None, tdm=True):
    return WriteAddrInstruction('WRITEADDR', channel, 0, addr, value, tdm)

def Invalidate(addr, nmeas, channel=None, tdm=True):
    return WriteAddrInstruction('INVALIDATE', channel, 1, addr, 2**nmeas-1, tdm)

def CrossBar(addr, value, channel=None, tdm=True): # should this be a high-level instruction though?
    return WriteAddrInstruction('CROSSBAR', channel, 3, addr, value, tdm)

def StoreMeas(addr, value, channel=None, tdm=True):
    return WriteAddrInstruction('STOREMEAS', channel, 5, addr, value, tdm)

class LoadCmpVramInstruction(object):

    def __init__(self, name, use_vram, addr, mask, tdm):
        # TODO: sanity checks on input values
        self.instruction = name
        self.use_vram = use_vram
        self.mask = mask
        self.addr = addr
        self.tdm = tdm
        self.length = 0

    def promote(self, ptype):
        return self

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self == other


# def LoadCmpVram(addr, mask):
#     return LoadCmpVramInstruction('LOADCMPVRAM', 1, addr, mask)
