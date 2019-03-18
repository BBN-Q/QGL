'''
Copyright 2019 Raytheon BBN Technologies

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

from .TdmInstructions import CustomInstruction, WriteAddr, Invalidate, LoadCmpVram
from .ControlFlow import LoadCmp, Call

##custom APS OP_CODES
#APS_CLIFFORD_SET_SEED:         Set 32 bit seed for random number generator
#APS_CLIFFORD_SET_OFFSET:       Stores jump table starting offset
#APS_CLIFFORD_SET_SPACING:      Stores jump table spacing (integer to shift left),
#                               including interal address representation
#APS_CLIFFORD_INVERSE_RESET:    Reset inverse tracker
#APS_CLIFFORD_RAND:             Choose randomly from set of 24 waveforms,
#                               returns starting address of jump table entry from RAM
#APS_CLIFFORD_INVERSE:          Store waveform address of inverse waveform to be jumped to

APS2_CUSTOM = {"APS_RAND": 0,
                     "APS_CLIFFORD_RAND": 2,
                     "APS_CLIFFORD_INVERSE": 6,
                     "APS_CLIFFORD_INVERSE_RESET": 7,
                     "APS_CLIFFORD_SET_SEED": 3,
                     "APS_CLIFFORD_SET_OFFSET": 4,
                     "APS_CLIFFORD_SET_SPACING": 5}

APS2_CUSTOM_DECODE = {v: k for k, v in APS2_CUSTOM.items()}

##Note that the expected call pattern for the randomizer is:

#   APS_CLIFFORD_RAND / APS_CLIFFORD_INVERSE
#   LOADCMP
#   CALL

def RandomCliffordSetOffset(addr, offset):
    return [Invalidate(addr, 0, tdm=False),
            WriteAddr(addr, offset, tdm=False),
            LoadCmpVram(addr, 0xFFFFFFFF, tdm=False), #
            CustomInstruction("APS_CLIFFORD_SET_OFFSET", addr, 0xA)]

def RandomCliffordSetSpacing(addr, spacing):
    return [Invalidate(addr, 0, tdm=False),
            WriteAddr(addr, spacing, tdm=False),
            LoadCmpVram(addr, 0xFFFFFFFF, tdm=False),
            CustomInstruction("APS_CLIFFORD_SET_SPACING", addr, 0xA)]

def RandomClifford(target, addr):
    return [Invalidate(addr, 0, tdm=False),
            CustomInstruction("APS_CLIFFORD_RAND", 0x0, addr),
            LoadCmpVram(addr, 0xFFFFFFFF, tdm=False),
            Call(target, indirect=True)]

def RandomCliffordInverse(target, addr):
    return [Invalidate(addr, 0, tdm=False),
            CustomInstruction("APS_CLIFFORD_INVERSE", 0x0, addr),
            LoadCmpVram(addr, 0xFFFFFFFF, tdm=False),
            Call(target, indirect=True)]

def RandomCliffordInverseReset(addr): #for now don't use addr
    return CustomInstruction("APS_CLIFFORD_INVERSE_RESET", 0x0, 0xF)

def RandomCliffordSeed(seed):
    return [Invalidate(0xB, 0, tdm=False),
        WriteAddr(0xB, seed, tdm=False),
        LoadCmpVram(0xB, 0xFFFFFFFF, tdm=False),
        CustomInstruction("APS_CLIFFORD_SET_SEED", 0xB, 0xA)]
