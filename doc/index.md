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
mapping is [described later](config.md#channel-library-setup).

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
[config file](config.md#configuration-options).

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

### Additional Parameters

All QGL pulse primitives accept an arbitrary number of additional keyword
arguments. In particular, any QGL primitive accepts a `length` keyword to modify
the length of the resulting operation. These additional parameters are passed to
the [shape function](#pulse-shapes-and-waveforms) when the QGL compiler
constructs waveforms from `Pulse` objects.

## Sequences and Simultaneous Operations

Programs in QGL are specified using python lists. For example,
```python
q1 = QubitFactory("q1")
seq = [X90(q1), X(q1), Y(q1), X90(q1), MEAS(q1)]
```

The elements of the list provide a time-ordered sequence of pulses to execute.
Using the python list to describe sequences allows for the use of python's
powerful list comprehension syntax to describe sequence variations. For
instance, you can concisely write a scan over a rotation angle or delay in a
list comprehension such as:
```python
seq = [[X90(q1), Id(q1, length=d), X90(q1), MEAS(q1)] for d in np.linspace(0, 10e-6, 11)]
```
QGL's compiler assumes that such lists of lists represent a series of related
experiments and schedules them to occur sequentially in the AWG output.

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

## Composite Pulses

Occasionally one wants to construct a sequence of pulses and treat them as if
the entire sequence were a single pulse. For this, QGL allows pulses to be
joined with the `+` operator. This allows, for example, us to define
```python
def hadamard(q):
    return Y90(q) + X(q)
```
and then use `hadmard(q)` just like any other pulse primitive, even though it is
composed of a sequence of two pulses.

## Pulse Shapes and Waveforms

The QGL compiler constructs waveforms to implement the desired quantum
operations. To do this, each pulse has a `shapeFun` (shape function) that is
called with its `shapeParams`. A number of commonly used shapes are defined in
the `PulseShapes` module including:

* `constant` - i.e. a "square" pulse with constant amplitude
* `tanh` - essentially a square pulse with rounded edges
* `gaussian` - a truncated Gaussian shape
* `drag` - the DRAG pulse gives a Gaussian shape with its derivative on the opposite quadrature.
* `gaussOn` - the first half of a truncated Gaussian shape
* `gaussOff` - the second half of a truncated Gaussian shape

The default pulse shape is determined by properties in the [channel
library](config.md#channel-library-setup). However, the QGL programmer may
override the default shape with a keyword argument. For example, to force the
use of square pulse shape we may write:
```python
seq = [X(q1, shapeFun=PulseShapes.constant), MEAS(q1)]
```

One common use case for specifying a shape function is in the construction of
composite pulses. For instance, you may want a square pulse shape with Gaussian
edges rather than those given by the `tanh` function. To do this you might write:
```python
seq = [X(q1, shapeFun=PulseShapes.gaussOn) +\
       X(q1, shapeFun=PulseShapes.constant) +\
       X(q1, shapeFun=PulseShapes.gaussOff),
       MEAS(q1)]
```

Shape functions can be an arbitrary piece of python code that returns a NumPy
array of complex values. Shape functions must accept **all** of their arguments
as keyword arguments. The only arguments that are guaranteed to exist are
`samplingRate` and `length`. The pulse length is always assumed to be units of
seconds; it is up to the shape function to use the passed sampling rate to
convert from time into number of points/samples. As an example, we could define
a ramp shape with
```python
def ramp(length=0, samplingRate=1e9, **kwargs):
    numPts = int(np.round(length * samplingRate))
    return np.linspace(0, 1, numPts)
```

Then use it with any pulse primitive, e.g.:
```python
seq = [X(q1, shapeFun=ramp)]
```

If your custom shape function requires additional arguments, you must either
arrange for these parameters to exist in the `LogicalChannel`'s `shapeParams`
dictionary, or pass them at the call site. For instance,
```python
def foo(length=0, samplingRate=1e9, bar=1, **kwargs):
    numPts = int(np.round(length * samplingRate))
    # something involving bar...

seq = [X(q1, bar=0.5, shapeFun=foo)] # bar is passed as a keyword arg
```

See the `PulseShapes` module for further examples.

## Compiling and Plotting

To reduce a pulse sequence to AWG vendor-specific hardware instructions, use the
`compile_to_hardware()` method, e.g.:
```python
seq = [[X90(q1), Id(q1, length=d), X90(q1), MEAS(q1)] for d in np.linspace(0, 10e-6, 11)]
meta_info = compile_to_hardware(seq, 'test/ramsey')
```

This code snippet will create a folder called `test` inside
[`AWGDir`](config.md#configuration-options) and produce sequence files for each
AWG targeted by the `PhysicalChannels` associated with the QGL program. For
instance, if the `q1` channel targeted an AWG named `APS1` and the `M-q1`
channel targeted `APS2`, then the above call to `compile_to_hardware` would
produce two files: `ramsey-APS1.h5` and `ramsey-APS2.h5` in the `test` folder.
It would also produce a *meta information* file `ramsey-meta.json` which
contains data about the QGL program that may be useful for executing the
program in an instrument control platform such as
[Auspex](https://github.com/BBN-Q/Auspex). `compile_to_hardware` returns the
path to this meta info file.

The `plot_pulse_files()` creates a visual representation of the pulse sequence
created by a QGL program. For example,
```python
plot_pulse_files(meta_info)
```
will create an interactive plot where each line represents a physical output
channel of an AWG referenced by the QGL program.

You may also view a QGL program prior to the logical -> physical mapping with
`show()`. For example,
```python
seq = [X90(q1), Id(q1, length=100e-9), X90(q1), MEAS(q1)]
show(seq)
```
will create a plot grid where each subplot shows the operations on individual
`LogicalChannels`.

## Axis Descriptors

TODO
