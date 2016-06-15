# Getting Started with QGL

To produce a sequence file usable with an AWG (e.g. BBN APS1, APS2, or Tek
5014), QGL relies upon a mapping from *logical* resources (qubits, measurements,
markers) to *physical* resources (specific AWG channels). So, your first step is
to create a basic channel library with the appropriate logical and physical
channels. An example minimal library might include these *logical* channels:

* q1 (Qubit)
* M-q1 (Measurement)
* slaveTrig (LogicalMarker)
* digitizerTrig (LogicalMarker)

The naming convention for measurement channels of M-**qubitname** *must* be
followed. In addition, the `slaveTrig` and `digitizerTrig` channels are required
to compile a sequence.

To create and manage this library, you can use the `ExpSettingsGUI` from
[PyQLab](https://github.com/BBN-Q/PyQLab). The first tab, labeled "Channels" has
two sub-panels for logical and physical channels, respectively. The fastest way
to get started is to go to the "Instruments" tab, select the "AWG's" sub-tab,
and create one or more AWGs of the appropriate type. Then, from the File menu
choose "auto populate physical channels".  Finally, create a set of logical
channels (e.g. `q1` and `M-q1`) and select the appropriate physical channel for
each.

If you choose to create physical channels manually, note that these channel
names must follow the convention *AWGName*-**channel**, where *AWGName*
corresponds to the name of the corresponding AWG in the instrument library. The
available **channel** names depend on the type of AWG.

APS: `12, 34, 1m1, 2m1, 3m1, 4m1`  
APS2: `12, 12m1, 12m2, 12m3, 12m4`  
Tek5014: `12, 34, 1m1, 1m2, 2m1, 2m2, 3m1, 3m2, 4m1, 4m2`

You'll notice that QGL explicitly groups pairs of analog output channels into quadrature pairs for I/Q modulation of a carrier waveform. Support for single-channel real output is still a TODO.

## Example QGL usage

### Ramsey

A simple Ramsey experiment with 51 steps between 0 and 100us:
```python
from QGL import *
output_file() # or output_notebook() from a jupyter notebook

q1 = QubitFactory('q1')
seqs = [[X90(q1), Id(q1, length=d), X90(q1), MEAS(q1)] for d in np.linspace(0,100e-6, 51)]
# get a view of the 10th sequence grouped by logical channels
show(seqs[9])
# compile
compile_to_hardware(seqs)
```

The `QubitFactory` method is a convenience method to create a `Qubit` object
with parameters read from the pulse parameters file. The string `'q1'`
identifies the qubit label to look up in the parameter file. You can also create
`Qubit`s on the fly by passing in the parameters.

To define the QGL program, the example uses a python list comprehension to
express a list of sequences in a single line.

### Multi-qubit

Multi-qubit sequences are represented with notation similar to a tensor product.
To get this to work, you must add to your channel library a logical channel
representing the coupling between qubits. In QGL terminology, this is known as
an `Edge`, and is a *directed* edge in the connectivity graph of your device.
QGL uses directed edges because certain two-qubit interactions have a preferred
ordering of the interaction. For instance, a cross resonance gate has a
preferred sign of the qubit-qubit detuning. By storing directed edges, we can
write two-qubit primitives that emit different pulses depending on whether the
(control, target) pair is aligned or anti-aligned with the underlying
interaction Hamiltonian.

To see this in action, add two qubits to your channel library and one edge
**q1-q2** connecting them. Then try the following program:

```python
from QGL import *
output_file() # or output_notebook() from a jupyter notebook

q1 = QubitFactory('q1')
q2 = QubitFactory('q2')
seq = [X90(q1)*Y(q2), CNOT_CR(q1,q2), Y(q1), CNOT(q1,q2), X90(q1)*Y(q2), MEAS(q1)*MEAS(q2)]
show(seq)
compile_to_hardware(seq)
```

In this example you can see the use of the two-qubit primitive `CNOT_CR`, which
emits a sequence of 1- and 2-qubit gates necessary to compose a CNOT operation
from single-qubit gates and a ZX90. The exact sequence will depend on the
(source, target) order you selected in creating the **q1-q2** `Edge`.
