
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
        Id(q1),

        WriteAddr(1, 7),
        MajorityMask(1, 0),

        Invalidate(addr=10, mask=0x7),

        MEASA(q1, maddr=(10, 0)),
        # MEASA(q1, maddr=(10, 1)),
        # MEASA(q1, maddr=(10, 2)),
        MEASA(q1, maddr=(20, 0)),

        LoadCmpTdm(0xfedc, 0x1234678),

        MajorityVote(10, 11),
        ]

aps_metafile = compile_to_hardware([seq], '/tmp/f')
tdm_instr = QGL.drivers.APS2TDMPattern.get_tdm_instructions()

aps_metadata = json.loads(open(aps_metafile).read())
print(aps_metadata)

for key in aps_metadata['instruments']:
    print('')
    print('INSTRUMENT %s' % str(key))
    instructions = aps2_reader.raw_instructions(aps_metadata['instruments'][key])
    # print('TYPE %s' % str(type(instructions)))
    # aps2_reader.display_decompiled_instructions(instructions)

    for i in range(len(instructions)):
        instr_bits = instructions[i]
        instr_txt = str(Instruction.unflatten(instr_bits))
        print('%5d: 0x%.16x - %s' % (i, instr_bits, instr_txt))

print('')
print('INSTRUMENT tdm')
for i in range(len(tdm_instr)):
    instr_bits = tdm_instr[i]
    instr_txt = str(Instruction.unflatten(instr_bits))
    print('%5d: 0x%.16x - %s' % (i, instr_bits, instr_txt))

# TODO: insert TDM instructions into the output file



