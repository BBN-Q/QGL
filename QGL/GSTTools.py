'''
Various tools to interface with pyGSTi for running GST experiments.

Created on May 16, 2018

Original Author: Guilhem Ribeill

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
from .PulsePrimitives import *
from .Cliffords import *
from .BasicSequences.helpers import create_cal_seqs
from .Compiler import compile_to_hardware
from itertools import chain
from random import choices

PYGSTI_PRESENT = False
try:
   from pygsti.objects.circuit import Circuit
   PYGSTI_PRESENT = True
except:
   pass

#Default mapping from pyGSTi naming convention to QGL gates.
gst_gate_map = {"Gx": X90,
                "Gy": Y90,
                "Gi": Id}

def gst_map_1Q(gst_list, qubit, qgl_map=gst_gate_map, append_meas=True):
    """
    Helper function that takes an arbitrarily nested list of pygsti gatestrings
    and converts them into QGL sequences, keeping the same nesting of lists.

    Inputs:
        gst_list: GateString to convert, or possibly nested list of pyGSTi GateStrings.
        qubit: QGL qubit to apply the sequence to
        qgl_map: Dictionary that maps between pyGSTi "Gx" string to QGL pulse
        append_meas: Append a measurement to each sequence.
    Returns:
        QGL sequences, preserving the input list nesting (as a generator)
    """
    if isinstance(gst_list, Circuit):
        gst_list = [gst_list]
    for item in gst_list:
        if isinstance(item, Circuit):
            mapped = map(lambda x: qgl_map[str(x)](qubit), item.tup)
            if append_meas:
                yield list(chain(mapped, [MEAS(qubit)]))
            else:
                yield list(mapped)
        elif isinstance(item, list):
            yield list(gst_map_1Q(item, qubit, qgl_map=qgl_map, append_meas=append_meas))

def gst_map_2Q(gst_list, qubits, qgl_map=None, append_meas=False):
    """
    Helper function that takes an arbitrarily nested list of pygsti gatestrings
    and converts them into QGL sequences, keeping the same nesting of lists.

    Inputs:
        gst_list: GateString to convert, or possibly nested list of pyGSTi GateStrings.
        qubit: QGL qubit to apply the sequence to
        qgl_map: Dictionary that maps between pyGSTi "Gx" string to QGL pulse
        append_meas: Append a measurement to each sequence.
    Returns:
        QGL sequences, preserving the input list nesting (as a generator)
    """
    if isinstance(gst_list, GateString):
        gst_list = [gst_list]
    for item in gst_list:
        if isinstance(item, GateString):
            mapped = map(lambda x: qgl_map[x], item.tup)
            if append_meas:
                yield list(chain(mapped, [reduce(lambda x,y: x*y, map(MEAS, qubits))]))
            else:
                yield list(mapped)
        elif isinstance(item, list):
            yield list(gst_map_2Q(item, qubit, qgl_map=qgl_map, append_meas=append_meas))

def create_gst_sequence_from_pygsti(gst_list, qubit, gate_map=gst_gate_map):
    """ Returns list of QGL sequences from a pyGSTi GateString list. See gst_map_1Q.
        The return value is a list of sequences that can be complied by QGL.
    """
    return list(gst_map_1Q(gst_list, qubit, qgl_map=gate_map, append_meas=True))

def pygsti_to_cliffords(gst_seq):

    #Map from GST convention to cliffords
    cliff_map = {"{}": 0,
                 "Gi": 1,
                 "Gx": 2,
                 "Gy": 5}
    #convert to dictionary of lambdas for compatibility with gst_map_1Q
    lambda_map = {k: lambda x, v=v: v for k, v in cliff_map.items()}

    return list(gst_map_1Q(gst_seq, None, qgl_map=lambda_map,
                                    append_meas=False))

def pauli_rand_clifford_circuit(gst_seq):

    def seqreduce(s):
        if not s:
            return 0
        else:
            return reduce(lambda x,y: clifford_multiply(x,y), s)

    def inv_cliff(c):
        return inverse_clifford(clifford_mat(c, 1))

    c_ps = [0, 2, 5, 8]

    c_seqs = pygsti_to_cliffords(gst_seq)

    r_seqs = []
    for seq in c_seqs:
        if not seq:
            r_seqs.append([])
        else:
            rand_pauli = choices(c_ps, k=len(seq))
            inter = 0
            bare = 0
            rseq = []
            for j in range(len(seq)):
                inter = clifford_multiply(clifford_multiply(inter, rand_pauli[j]), seq[j])
                bare = clifford_multiply(bare, seq[j])
                rseq.append(clifford_multiply(rand_pauli[j], seq[j]))
            recovery = clifford_multiply(inv_cliff(inter), bare)
            rseq[-1] = clifford_multiply(rseq[-1], recovery)
            r_seqs.append(rseq)

    all_ok = all((r == i for r, i in zip(map(seqreduce, r_seqs), map(seqreduce, c_seqs))))
    assert all_ok, "Something went wrong when Pauli-frame randomizing!"

    return r_seqs

def SingleQubitCliffordGST(qubit, pygsti_seq, pulse_library="Standard", randomized=False, num_cals=100, diac_compiled=True):

    pulse_library = pulse_library.upper()

    # QGL pulse libraries handle the Id pulse differently.  In the standard
    # case, the Id is of finite length equal to all the other one-pulse
    # elements of the library.  In the Atomic and DiAtomic cases, the ID is
    # of length 0 by default.  In GST, we need access to both types of the ID
    # gate with the first experiment in any GST experiment  equal to {} =
    # Id(length = 0).  All other Id gates in the sequence should be of finite
    # length.  So we'll modify the Clifford indexing here to make Id(length=0)
    # the first element in the library and Id(length=length) the second.
    if pulse_library == "STANDARD":
        #clifford_pulse = lambda x: clifford_seq(x, qubit)
        clifford_pulse = [clifford_seq(i, qubit) for i in range(24)]
        clifford_pulse.insert(0, Id(qubit, length=0.0))
    elif pulse_library == "DIAC":
        #clifford_pulse = lambda x: DiAC(qubit, x, diac_compiled)
        clifford_pulse = [AC(qubit, i, diac_compiled) for i in range(24)]
        clifford_pulse.insert(1, Id(qubit))
    elif pulse_library == "AC":
        #clifford_pulse = lambda x: AC(qubit, x)
        clifford_pulse = [AC(qubit, i) for i in range(24)]
        clifford_pulse.insert(1, Id(qubit))
        raise ValueError("Pulse library must be one of 'standard', 'diac', or 'ac'. Got {} instead".format(pulse_library))

    if randomized:
        seqs = pauli_rand_clifford_circuit(pygsti_seq)
    else:
        seqs = pygsti_to_cliffords(pygsti_seq)

    qgl_seqs = []

    for seq in seqs:
        qgl_seqs.append([clifford_pulse[c] for c in seq])
        qgl_seqs[-1].append(MEAS(qubit))

    if num_cals != 0:
        qgl_seqs += create_cal_seqs((qubit, ), abs(num_cals))

    metafile = compile_to_hardware(qgl_seqs, 'GST/GST')
    return metafile
