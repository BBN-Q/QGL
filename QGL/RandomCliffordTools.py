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
from . import PatternUtils
from .PatternUtils import flatten, has_gate
from . import Channels
from . import ChannelLibraries
from . import PulseShapes
from . import PulsePrimitives
from .PulsePrimitives import Id, clear_pulse_cache
from .PulseSequencer import Pulse, PulseBlock, CompositePulse
from . import ControlFlow
from . import BlockLabel
from . import TdmInstructions # only for APS2-TDM
from . import APS2CustomInstructions

def generate_clifford_jump_table(cliff_wires, jt_label = None):
    """Generate the jump table that will be used to call into the clifford set"""

    if jt_label is None:
        jt_label = BlockLabel.BlockLabel("JT")
    if not isinstance(jt_label, BlockLabel.BlockLabel):
        raise ValueError("Jump table label must be a BlockLabel.")

    cliff_labels = [f"C{n}" for n in range(24)]
    jump_table = [jt_label]
    for cl in cliff_labels:
        jump_table.append(ControlFlow.Call(cl))
        jump_table.append(ControlFlow.Return())

    for cliff, cl in zip(cliff_wires, cliff_labels):
        cliff.insert(0, BlockLabel.BlockLabel(cl))
        cliff.append(ControlFlow.Return())

    cliff_wires.insert(0, jump_table)
