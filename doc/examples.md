## Ramsey Sequence

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
`Qubit`s on the fly by passing in parameters to the `Qubit` constructor.

To define the QGL program, the Ramsey example uses a python list comprehension
to express a list of sequences in a single line.

## Multi-qubit Sequence

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
seq = [X90(q1)*Y(q2), CNOT(q1,q2), Y(q1), CNOT(q1,q2), X90(q1)*Y(q2), MEAS(q1)*MEAS(q2)]
show(seq)
compile_to_hardware(seq)
```

In this example you can see the use of the two-qubit primitive `CNOT`, which by
default emits a sequence of 1- and 2-qubit gates necessary to compose a CNOT
operation from single-qubit gates and a ZX90, as is appropriate for a
cross-resonance interaction. The exact sequence will depend on the (source,
target) order you selected in creating the **q1-q2** `Edge`. You can select a
different default CNOT implementation by modifying the 'cnot_implementation' key
in QGL's config.json configuration file.

## Built-in Basic Sequences

QGL provides many pre-defined methods for sequences commonly used to
characterize a quantum device. These methods are defined in QGL's
`BasicSequences` package and include:

* `RabiAmp`
* `RabiWidth`
* `PulsedSpec`
* `InversionRecovery`
* `Ramsey`
* `HahnEcho`
* `CPMG`
* `SingleQubitRB`
* `TwoQubitRB`

Usage of each is defined in its respective doc string. For instance, at an
ipython prompt, you may type
```shell
In [1]: ?RabiAmp
```
to learn about the `RabiAmp` function. We encourage users to peruse the methods
defined in `BasicSequences` for templates that may be useful in writing their
own QGL programs.
