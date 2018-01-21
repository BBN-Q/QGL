
import QGL.drivers
from QGL import *
import json
import numpy as np
import pprint
import time
import aps2_reader

from QGL.drivers.APS2TDMPattern import Instruction
import QGL.drivers.APS2TDMPattern

# import aps2

ChannelLibrary(blank=True)

_q1 = Qubit(label='q1')
_q2 = Qubit(label='q2')

def setUp():
    # Copied from the unittests (CompileUtils).
    # Probably not valid, but OK for placeholders.

    """
    q1gate = Channels.LogicalMarkerChannel(label='q1-gate')
    q1 = Qubit(label='q1', gate_chan=q1gate)
    q1.pulse_params['length'] = 30e-9

    q2gate = Channels.LogicalMarkerChannel(label='q2-gate')
    q2 = Qubit(label='q2', gate_chan=q2gate)
    q2.pulse_params['length'] = 30e-9

    trigger = Channels.LogicalMarkerChannel(label='trigger')
    measq1 = Channels.Measurement(label='M-q1', meas_type='autodyne')
    measq1.trig_chan = trigger

    measq2 = Channels.Measurement(label='M-q2', meas_type='autodyne')
    measq2.trig_chan = trigger
    """

    cl = ChannelLibrary(library_file="./meas.yml")

    # ChannelLibrary(blank=True) # Create a blank ChannelLibrary
    """
    ChannelLibraries.channelLib.channelDict = {
            'q1': q1,
            'q2': q2,
            'M-q1': measq1,
            'M-q2': measq2
    }
    """
    ChannelLibraries.channelLib.build_connectivity_graph()

setUp()


q1 = QubitFactory('q1')
q2 = QubitFactory('q2')

seq = [
        X90(q1),
        # X90(q2),
        MEASA(q1, maddr=(10, 3)),
        WriteAddr(1, 7, channel=q1),
        Invalidate(addr=4, mask=0xfff),
        MajorityMask(1, 0),
        MajorityVote(10, 9),
        WriteAddr(12, 13, channel=q1)
        ]

aps_metafile = compile_to_hardware([seq], '/tmp/f')
tdm_instr = QGL.drivers.APS2TDMPattern.get_tdm_instructions()

aps_metadata = json.loads(open(aps_metafile).read())
print(aps_metadata)

for key in aps_metadata['instruments']:
    print('INSTRUMENT %s' % str(key))
    instructions = aps2_reader.raw_instructions(aps_metadata['instruments'][key])
    # print('TYPE %s' % str(type(instructions)))
    # aps2_reader.display_decompiled_instructions(instructions)

    print('')
    for i in range(len(instructions)):
        instr_bits = instructions[i]
        instr_txt = str(Instruction.unflatten(instr_bits))
        print('%5d: 0x%.16x - %s' % (i, instr_bits, instr_txt))

print('INSTRUMENT tdm')
for i in range(len(tdm_instr)):
    instr_bits = tdm_instr[i]
    instr_txt = str(Instruction.unflatten(instr_bits))
    print('%5d: 0x%.16x - %s' % (i, instr_bits, instr_txt))



