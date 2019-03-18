import unittest
import numpy as np
import time
import os.path
import struct
from bbndb.qgl import PhysicalChannel, LogicalChannel, Measurement
from QGL import *

# Waveform numpy assert_allclose Test for QGL
#
# Usage add sequences to the *TestCases classes
# Run *TestCases.show() to view waveforms to make sure they are correct
# Run *TestCases.write() to write out the valid waveform files
# Add the new files with git if desired
# Run *TestCases.validate() to validate outside of unit test framework
#
# There are currently two classes which subclass unittest.TestCase
# Which will run the validation as part of a unit test


class SequenceTestCases(object):
    # base class for sequence test cases
    testFileDirectory = './tests/test_data/'
    fileHeader = ''

    def __init__(self):
        super(SequenceTestCases, self).__init__()
        self.sequences = {}
        self.waveforms = {}
        self.validWaveforms = {}

        # default state generates, compiles, and loads
        # all available sequences (which are defined in subclasses)
        self.generate()
        self.compile()
        self.load()  

    def generate(self):
        # this function should be overridden in a subclass to define
        # the sequences
        pass

    def compile(self):
        # compiles all available sequences
        for name, seq in self.sequences.items():
            self.waveforms[name] = build_waveforms(seq)

    def show(self):
        # display all sequences for human verification
        # of correctness
        for name, waveform in self.waveforms.items():
            plot_waveforms(waveform, name + ': ')
            # give bokeh time to generate the plot
            time.sleep(1)

    def build_filename(self, name):
        # utility for reading and writing
        return self.testFileDirectory + self.fileHeader + '-' + name + '.aps2'

    def write(self):
        # writes each sequence to a file in the test data directory a file header
        # which may be overridden in the subclasses
        for name, waveform in self.waveforms.items():

            fileName = self.build_filename(name)
            if os.path.isfile(fileName):
                os.remove(fileName)
            with open(fileName, 'wb') as FID:
                FID.write(b'TEST')                     # target 
                FID.write(np.float32(1.0).tobytes())   # Version
                FID.write(np.float32(1.0).tobytes())   # minimum firmware version
                FID.write(np.uint16(len(waveform)).tobytes())      # number of channels

                channels = waveform.keys()
                for channel in channels:
                    channel.label.ljust(32,"#").encode("utf-8")
                    FID.write(channel.label.ljust(32,"#").encode("utf-8"))
                    FID.write(np.uint64(waveform[channel].size).tobytes())
                    FID.write(np.complex128(waveform[channel]).tobytes()) # waveform data length for channel


    def load(self):
        # attempts to load valid waveforms for each of the sequences test cases
        for caseName in self.sequences:
            validWaveform = {}
            fileName = self.build_filename(caseName)
            if not os.path.isfile(fileName):
                print(
                    'Warning: valid waveform file for {0} not found at: {1}'.format(
                        caseName, fileName))
                continue
            # -----
            # print( "\n\rDBG::Calling \"with open({0}, 'wb')\"...".format( fileName) )
            # -----
            # The following pulls in Git Large File Storage (LFS) data files
            # from cached signature references;  if the open call fails
            # with an OSError, double-check git-lfs library installation (in
            # addition to git) -- <https://git-lfs.github.com/>
            #
            # Where git-lfs was NOT installed the open call was observed
            # returning:
            #
            # OSError: Unable to open file (file signature not found))
            #
            with open(fileName, 'rb') as FID:
                target_hw    = FID.read(4).decode('utf-8')
                file_version = struct.unpack('<f', FID.read(4))[0]
                min_fw       = struct.unpack('<f', FID.read(4))[0]
                num_chans    = struct.unpack('<H', FID.read(2))[0]
                assert target_hw=="TEST", f"Cannot compare to file type {target_hw}. Expecting TEST"

                for i in range(num_chans):
                    wf_name = FID.read(32).decode('utf-8').replace('#','')
                    wf_len  = struct.unpack('<Q', FID.read(8))[0]
                    wf_dat  = np.frombuffer(FID.read(16*wf_len), dtype=np.complex128)
                    validWaveform[wf_name] = wf_dat
            self.validWaveforms[caseName] = validWaveform

    def validate(self):
        for caseName in self.sequences:
            self.validate_case(caseName)

    def validate_case(self, caseName):
        # validates each sequences by using numpy assert_allclose for each channel
        assert (caseName in self.validWaveforms)
        validWaveform = self.validWaveforms[caseName]
        for channel, waveform in self.waveforms[caseName].items():
            print('Validating {0} Case {1} Channel {2}'.format(
                self.__class__.__name__, caseName, str(channel)))
            
            assert (channel.label in validWaveform)
            np.testing.assert_allclose(waveform,
                                       validWaveform[channel.label],
                                       rtol=1e-5,
                                       atol=0)


class SingleQubitTestCases(SequenceTestCases):
    # Single Qubit Sequence Test Cases

    fileHeader = 'single'

    def newQ1(self):
        cl = ChannelLibrary(db_resource_name=":memory:")
        cl.clear()
        q1 = cl.new_qubit("q1")
        q1.pulse_params['length'] = 30e-9
        q1.measure_chan = Measurement(label="test_meas")
        q1.phys_chan = PhysicalChannel(label="test_phys_c", sampling_rate=1.2e9) 
        q1.measure_chan.phys_chan = PhysicalChannel(label="test_phys", sampling_rate=1.2e9)
        return q1

    def generate(self):
        q1 = self.newQ1()
        self.sequences['ramsey'] = [[X90(q1), Id(q1, delay), X90(q1)]
                                    for delay in np.linspace(0.0, 1e-6, 11)]

        self.sequences['repeat'] = [X90(q1)] + repeat(5, Y(q1)) + [X90(q1)]

class MultiQubitTestCases(SequenceTestCases):
    # Multi Qubit Sequence Test Cases

    fileHeader = 'multi'

    def newQubits(self):
        # Create an in-memory blank channel library
        cl = ChannelLibrary(db_resource_name=":memory:")
        cl.clear()
        q1 = cl.new_qubit("q1")
        q2 = cl.new_qubit("q2")
        q1.pulse_params['length'] = 30e-9
        q1.measure_chan = Measurement(label="test_meas1")
        q1.phys_chan = PhysicalChannel(label="test_phys_c1", sampling_rate=1.2e9) 
        q1.measure_chan.phys_chan = PhysicalChannel(label="test_phys1", sampling_rate=1.2e9) 
        q2.pulse_params['length'] = 30e-9
        q2.measure_chan = Measurement(label="test_meas2")
        q2.phys_chan = PhysicalChannel(label="test_phys_c2", sampling_rate=1.2e9) 
        q2.measure_chan.phys_chan = PhysicalChannel(label="test_phys2", sampling_rate=1.2e9) 
        q1q2 =cl.new_edge(q1, q2)
        q1q2.phys_chan = PhysicalChannel(label="test_phys_12", sampling_rate=1.2e9) 
        return (q1, q2)

    def generate(self):
        q1, q2 = self.newQubits()
        self.sequences['operators'] = [X90(q1), X(q1) * Y(q2), CNOT_simple(q1, q2),
                                       Xm(q2), Y(q1) * X(q2)]

        self.sequences['align'] = [align('right',
            X90(q1), Xtheta(q2, amp=0.5, length=100e-9)), Y90(q1) * Y90(q2)]

        flipFlop = X(q1) + X(q1)
        self.sequences['composite'] = [flipFlop*Y(q2), Y90(q1)]

# unittest classes


class SingleQubit(unittest.TestCase):
    # Single Qubit  Unit Test

    def setUp(self):
        self.sequences = SingleQubitTestCases()

    def test_sequences(self):
        self.sequences.validate()


class MultiQubit(unittest.TestCase):
    # Multi Qubit  Unit Test
    def setUp(self):
        self.sequences = MultiQubitTestCases()

    def test_sequences(self):
        self.sequences.validate()


### Utilty functions which are intended to be run from a repl
def show_test_cases():
    SingleQubitTestCases().show()
    MultiQubitTestCases().show()


def write_test_cases():
    SingleQubitTestCases().write()
    MultiQubitTestCases().write()


def validate():
    SingleQubitTestCases().validate()
    MultiQubitTestCases().validate()


if __name__ == "__main__":
    unittest.main()
