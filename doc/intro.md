# Getting Started with QGL

Quantum Gate Language (QGL) is a domain specific language embedded in python for
specifying gate sequences on quantum processors. It is a low-level language in
the sense that users write programs at the level of gates on physical qubits.
While the QGL developers pay particular attention to the requirements of
superconducting qubit systems, its structure is generic to any qubit system with
dipole-coupled qubits in a rotating frame. In such systems, the rotation axis in
the X-Y plane is determined by the pulse phase, and Z-axis rotations may be
achieved through *frame updates*.

## Channels

Many early quantum processors require non-uniform control parameters to achieve
high-fidelity gates across all qubits in the device. To deal with this, QGL
provides a number of *channel* objects to store the *individual* parameters
needed to control or measure the qubits.

Gates in QGL operate on resources known as `LogicalChannels`. For control, these
channels are either `Qubits` or `Edges`. The `Qubit` channels encode properties
specific to manipulating individual qubits in quantum processor, while `Edges`
encode the connectivity of the device. Since some 2-qubit gates have a preferred
directionality due to the physical parameters of the device, `Edges` correspond
to a *directed* edge in the qubit connectivity graph. Qubit measurements instead
act upon a third `LogicalChannel` type which is the `Measurement` channel. A
final logical resource type is the `LogicalMarkerChannel` which is used to carry
ancillary pulse information for things such as event triggers.

All `LogicalChannels` must be mapped to `PhysicalChannels` in order for the QGL
compiler to produce sequence files for the target hardware. The setup of this
mapping is described later in the section on **Channel Library Setup**.

While setup of these channels is important for final sequence compilation, QGL
programs typically refer only to `Qubit` channels. Actions on other channel
types may be implied by the operation. For example, to create a `Qubit` object
in QGL, one can write:
```python
q1 = QubitFactory("q1")
```

The `QubitFactory` method returns a `Qubit` object with the properties defined
by the name `q1` if found in the channel library. If the name is not found, then
the users gets a `Qubit` object with the default properties.

## Gate Primitives

The underlying representation of all QGL operations is a `Pulse` object.
However, users are not expected to create `Pulses` directly, but instead
interact with various pre-defined one- and two-qubit primitives.

### Single-Qubit Operations

QGL provides the following single-qubit gates:
```python
# generic rotation angle and phase
Utheta(q, angle, phase)

# generic rotations about a specific axis (phase)
Xtheta(q, angle)
Ytheta(q, angle)
Ztheta(q, angle)

# generic rotations of a specific angle
U(q, phase)
U90(q, phase)

# rotations of particular angle and phase
# X (phase = 0)
X(q)    # rotates by pi (180 degrees)
X90(q)  # rotates by +pi/2
X90m(q) # rotates by -pi/2

# Y (phase = pi/2)
Y(q)
Y90(q)
Y90m(q)

# just frame-updates
Z(q)
Z90(q)
Z90m(q)

# identity (delay or no-op)
Id(q, length) # length parameter is optional

# measurement
MEAS(q)
```

Due to the utility of Clifford-group operations in characterizing gate
performance, QGL also directly provides a primitive to implement the 24-element
single-qubit Clifford group:
```python
# atomic Clifford operation on 1-qubit
AC(q, n)
```

This method is "atomic" because it implements the full 1-qubit Clifford group
with one pulse per element, as opposed to requiring a sequence of the primitives
already given above. We known of no canonical way to specify the elements of the
Clifford group; consequently, `AC` identifies which Clifford by a numerical
index (0-23). See the definition of `AC` in `PulsePrimitives.py` or the
definition of `C1` in `Cliffords.py` to find our enumeration of the group.

### Two-qubit Operations

QGL provides only one high-level two-qubit primitives, `CNOT`. The implementation
of CNOT may be chosen by specifying the `cnot_implementation` key in QGL's
config file.

```python
# high-level primitives
CNOT(q1, q2)

# mid-level primitives
CNOT_simple(q1, q2) # a single-pulse gate on Edge(q1, q2)
CNOT_CR(q1, q2)     # an "echoed" cross-resonance CNOT gate on Edge(q1, q2)
ZX90_CR(q1, q2)     # a ZX90 on Edge(q1, q2) implemented with "echoed"
                    # cross-resonance pulses

# lowest-level primitives
echoCR(q1, q2)  # A "echoed" cross-resonance pulse
```

## Sequences and Simultaneous Operations

Programs in QGL are specified using python lists. For example,
```python
q1 = QubitFactory("q1")
seq = [X90(q1), X(q1), Y(q1), X90(q1), MEAS(q1)]
```

The elements of the list provide a time-ordered sequence of pulses to execute.
Users express simultaneity in QGL using the `*` operator. For instance,
```python
q1 = QubitFactory("q1")
q2 = QubitFactory("q2")
seq = [X90(q1)*X90(q2), MEAS(q1)*MEAS(q2)]
```

would execute the same sequence on `Qubit`s `q1` and `q2`. If the gate durations
differ between `q1` and `q2`, the QGL compiler injects delays to create aligned
`PulseBlock`s. By default, simultaneous pulses are "left-aligned", meaning that
the leading pulse edges are aligned and padding delays are injected on the
trailing edge. However, the user may change this behavior with the `align`
method:
```python
seq = [align(X90(q1)*X90(q2)), align(MEAS(q1)*MEAS(q2), mode="right")]
```

`align` takes a `mode` argument ("left", "right", or default "center") to
specify a particular pulse alignment within a `PulseBlock`.


## Channel Library Setup

To produce a sequence file usable with an AWG (e.g. BBN APS1, APS2, or Tek
5014), QGL relies upon a mapping from *logical* resources (qubits, measurements,
markers) to *physical* resources (specific instrument channels). So, if you
intend to produce program outputs usable by real hardware, you must first create
a basic channel library with the appropriate logical and physical channels. An
example minimal library might include these *logical* channels:

* q1 (Qubit)
* M-q1 (Measurement)
* slaveTrig (LogicalMarker)

The naming convention for measurement channels of M-**qubitname** *must* be
followed. The `slaveTrig` channel is used in master-slave AWG configurations
where one master AWG is used to trigger the start of the slave AWGs.

To create and manage this library, one can use the `ExpSettingsGUI` from
[PyQLab](https://github.com/BBN-Q/PyQLab). The first tab, labeled "Channels" has
two sub-panels for logical and physical channels, respectively. The fastest way
to get started is to go to the "Instruments" tab, select the "AWG's" sub-tab,
and create one or more AWGs of the appropriate type. Then, on the "Physical"
sub-tab of "Channels", click the "Auto" button to populate a list of physical
channels.  Finally, create a set of logical channels (e.g. `q1` of type `Qubit`,
`M-q1` of type `Measurement`, and `slaveTrig` of type `LogicalMarker`) and
select the appropriate physical channel for each.

If you choose to create physical channels manually, note that these channel
names must follow the convention *AWGName*-**channel**, where *AWGName*
corresponds to the name of the corresponding AWG in the instrument library. The
available **channel** names depend on the type of AWG.

APS: `12, 34, 1m1, 2m1, 3m1, 4m1`  
APS2: `12, 12m1, 12m2, 12m3, 12m4`  
Tek5014: `12, 34, 1m1, 1m2, 2m1, 2m2, 3m1, 3m2, 4m1, 4m2`

You'll notice that QGL explicitly groups pairs of analog output channels into
quadrature pairs for I/Q modulation of a carrier waveform. Support for
single-channel real output is still a TODO.

## Configuration Options

Lorem ipsum...

## Pulse Shapes

Lorem ipsum...

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
`Qubit`s on the fly by passing in parameters to the `Qubit` constructor.

To define the QGL program, the Ramsey example uses a python list comprehension
to express a list of sequences in a single line.

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
