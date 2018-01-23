# Example program for creating TDM instruction files

import copy
import json
import numpy
import os
import shutil

import QGL.drivers
from QGL import *
import aps2_reader

from QGL.drivers.APS2TDMPattern import Instruction
import QGL.drivers.APS2TDMPattern

def pp_instructions(name, instructions):
    print('')
    print('INSTRUMENT: ' + name)
    for i in range(len(instructions)):
        instr_bits = instructions[i]
        instr_txt = str(Instruction.unflatten(instr_bits))
        print('%5d: 0x%.16x - %s' % (i, instr_bits, instr_txt))

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

ChannelLibrary(blank=True)
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
        # MEASA(q2, maddr=(20, 0)),

        # LoadCmpTdm(0xfedc, 0x1234678),

        MajorityVote(10, 11),
        ]

# First, compile for the APS units.  As a side effect,
# this creates the TDM instructions, but does NOT
# put them into the APS output file.  We retrieve
# the TDM instructions, and then recompile the instructions
# to create a template TDM file, and then insert the
# TDM instructions into the template.
#
# Note that the steps of compiling for the TDM and
# inserting the instructions into the file are FRAGILE
# because they require that there be an instrument named
# "BBNAPS1" in the machine config.  This name has special
# meaning.

# IMPORTANT: compilation is destructive: it modifies
# the input sequences (and possibly the instances in those
# sequences.  So, running compiler_to_hardware on the
# same sequence twice can FAIL (or give weird results).
# So copy the input sequence each time we use it...

aps_metafile = compile_to_hardware([copy.copy(seq)], '/tmp/aps')
tdm_instr = QGL.drivers.APS2TDMPattern.get_tdm_instructions()
aps_metadata = json.loads(open(aps_metafile).read())

tdm_metafile = compile_to_hardware([copy.copy(seq)], '/tmp/tdm')
tdm_metadata = json.loads(open(tdm_metafile).read())
aps2_reader.replace_instructions(
        tdm_metadata['instruments']['BBNAPS1'],
        numpy.array(tdm_instr, dtype=np.uint64))

for key in aps_metadata['instruments']:
    instructions = aps2_reader.raw_instructions(
            aps_metadata['instruments'][key])
    pp_instructions(str(key), instructions)

# Read the TDM instructions from file, just to make sure
#
tdm_instr_from_file = aps2_reader.raw_instructions(
        tdm_metadata['instruments']['BBNAPS1'])
pp_instructions('tdm', tdm_instr_from_file)

