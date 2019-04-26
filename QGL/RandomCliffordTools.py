'''
Functions for dealing with random clifford pulse sequences

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

from . import config
from functools import reduce
from . import PatternUtils
from .PatternUtils import flatten, has_gate
from . import Channels
from . import ChannelLibraries
from . import PulseShapes
from . import PulsePrimitives
from . import Compiler
from .PulsePrimitives import Id, clear_pulse_cache
from .PulseSequencer import Pulse, PulseBlock, CompositePulse
from . import ControlFlow
from . import BlockLabel
from . import TdmInstructions # only for APS2-TDM
from .APS2CustomInstructions import *
from .Cliffords import C1, inverse_clifford, clifford_multiply

default_clifford_options = {"offset": 0x0, "spacing": 0x1, "seed": 0x31}

VALID_CLIFFORD_TYPES = ('RandomAC',)

def generate_clifford_jump_table(cliff_wires, jt_label = None):
    """Generate the jump table that will be used to call into the clifford set"""

    if jt_label is None:
        jt_label = BlockLabel.BlockLabel("JT", jump_table=True, table_size=48)
    if not isinstance(jt_label, BlockLabel.BlockLabel):
        raise ValueError("Jump table label must be a BlockLabel.")

    cliff_labels = [BlockLabel.BlockLabel(f"C{n}") for n in range(24)]
    jump_table = [jt_label]
    for cl in cliff_labels:
        jump_table.append(ControlFlow.Call(cl))
        jump_table.append(ControlFlow.Return())

    for cliff, cl in zip(cliff_wires, cliff_labels):
        cliff.insert(0, cl)
        cliff.append(ControlFlow.Return())

    cliff_wires.insert(0, jump_table)

    return jt_label

def insert_clifford_calls(seqs, jt_label=None, cliff_addr=0x3, add_inv = True,
                            inv_addr=0x4, inv_values=[],
                            clifford_options=default_clifford_options):

    if jt_label is None:
        jt_label = BlockLabel.BlockLabel("JT")
    if not isinstance(jt_label, BlockLabel.BlockLabel):
        raise ValueError("Jump table label must be a BlockLabel.")

    #Insert defaults
    for k,v in default_clifford_options.items():
        if k not in clifford_options.keys():
            cliford_options[k] = v

    for idx, seq in enumerate(seqs[:]):
        if not PatternUtils.contains_runtime_pulses(seq):
            continue

        new_seq = []

        for pulse in seq:
            if isinstance(pulse, Pulse) and pulse.isRunTime \
                and pulse.label in VALID_CLIFFORD_TYPES:
                has_random_cliff = True
                new_seq.extend(RandomClifford(jt_label, cliff_addr))
                #print("Inserting clifford pulse!")
            else:
                new_seq.append(pulse)

        # if isinstance(seq[0], BlockLabel.BlockLabel):
        #     new_seq[1:1] = info_seqs
        # else:
        #     new_seq[0:0] = info_seqs

        if add_inv:
            #insert reset after first wait
            w_idx = next(i for i, v in enumerate(new_seq) if isinstance(v, ControlFlow.Wait))
            new_seq.insert(w_idx+1, RandomCliffordInverseReset(0x0))
            #insert at end of sequence or before last GOTO
            if isinstance(new_seq[-1], ControlFlow.Goto):
                  new_seq[-1:-1] = RandomCliffordInverse(jt_label, inv_addr)
            else:
                  new_seq.extend(RandomCliffordInverse(jt_label, inv_addr))

        seqs[idx] = new_seq

    info_seqs = []
    info_seqs.extend(RandomCliffordSetOffset(0x1, clifford_options["offset"]))
    info_seqs.extend(RandomCliffordSetSpacing(0x2, clifford_options["spacing"]))
    info_seqs.extend(RandomCliffordSeed(clifford_options["seed"]))
    seqs.insert(0, info_seqs)


def randomize_clifford_sequences(qubit, seqs, clifford_type=PulsePrimitives.RandomAC):

    c_seqs = [[clifford_type(qubit) for _ in seq] for seq in seqs]
    inverses = [inverse_clifford(C1[reduce(clifford_multiply, seq)]) for seq in seqs]
    return c_seqs, inverses
