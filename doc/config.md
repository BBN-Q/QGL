## Channel Library Setup

To produce a sequence file usable with an AWG (e.g. BBN APS1, APS2, or Tek
5014), QGL relies upon a mapping from *logical* resources (qubits, measurements,
markers) to *physical* resources (specific instrument channels). So, if you
intend to produce program outputs usable by real hardware, you must first create
a basic channel library with the appropriate logical and physical channels. An
example minimal library might include these *logical* channels:

* q1 (Qubit)
* M-q1 (Measurement)
* slave_trig (LogicalMarker)

The naming convention for measurement channels of M-**qubitname** *must* be
followed. The `slave_trig` channel is used in master-slave AWG configurations
where one master AWG is used to trigger the start of the slave AWGs.

To create and manage this library, one can use the `ExpSettingsGUI` from
[PyQLab](https://github.com/BBN-Q/PyQLab). The first tab, labeled "Channels" has
two sub-panels for logical and physical channels, respectively. The fastest way
to get started is to go to the "Instruments" tab, select the "AWG's" sub-tab,
and create one or more AWGs of the appropriate type. Then, on the "Physical"
sub-tab of "Channels", click the "Auto" button to populate a list of physical
channels.  Finally, create a set of logical channels (e.g. `q1` of type `Qubit`,
`M-q1` of type `Measurement`, and `slave_trig` of type `LogicalMarker`) and
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

QGL depends on a number of global configuration options defined in `config.json`
in the QGL package directory. If this file cannot be found when the package is
loaded, it is created and populated with these default parameters
```json
{
	"AWGDir": "/my/path/to/awg/",
	"ChannelLibraryFile": "/my/path/to/ChannelParams.json",
	"PulsePrimitivesLibrary": "standard",
	"cnot_implementation": "CNOT_CR",
	"PlotBackground": "#EAEAF2",
	"GridColor": "white"
}
```

These parameters control:

* `AWGDir` - the directory where sequence file outputs of the compiler should be
    saved.
* `ChannelLibraryFile` - path to the channel parameters JSON file
* `PulsePrimitivesLibrary` - `standard` or `all90`; the latter uses only 90 rotations
    to create all operations (e.g. `X = X90 + X90`)
* `cnot_implementation` - names another method to delegate implementation of
    `CNOT`; usually `CNOT_simple` or `CNOT_CR`
* `PlotBackground` - hex string controlling background color of plots
* `GridColor` - hex string controlling color of grid lines in plots
