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

## Pulse Shapes

Lorem ipsum...
