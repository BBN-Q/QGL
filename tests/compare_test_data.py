import os
import sys
import glob
import h5py
from builtins import input

from QGL import *
import QGL

BASE_AWG_DIR = QGL.config.AWGDir
BASE_TEST_DIR = './test_data/awg/'

def compare_sequences():
    test_subdirs = ['TestAPS1', 'TestAPS2']
    for subdir in test_subdirs:
        testdirs = glob.glob(os.path.join(BASE_TEST_DIR, subdir, '*'))
        for test in testdirs:
            # build up subdirectory name
            _,name = os.path.split(test)
            testfiles = glob.glob(os.path.join(test, '*'))
            # recurse into subdirectories
            while len(testfiles) == 1 and os.path.isdir(testfiles[0]):
                _,subname = os.path.split(testfiles[0])
                name = os.path.join(name, subname)
                testfiles = glob.glob(os.path.join(testfiles[0], '*'))
            newpath = os.path.join(BASE_AWG_DIR, subdir, name)
            print("{0} comparing to {1}".format(test, newpath))
            newfiles = glob.glob(os.path.join(newpath, '*'))
            #filter py27 look for py27 versions
            testfiles = filter_py27(testfiles)
            newfiles = filter_py27(newfiles)
            PulseSequencePlotter.plot_pulse_files_compare(testfiles, newfiles)
            c = input('Enter to continue (q to quit): ')
            if c == 'q':
                break

def filter_py27(filenames):
    py27_files = [f for f in filenames if f[-8:] == "_py27.h5"]
    if py27_files:
        if sys.version_info[0] > 2:
            #strip them out for python3 (and above?)
            filenames = [f for f in filenames if f[-8:] != "_py27.h5"]
        else:
            #otherwise strip the non "_py27" version
            filenames = [f for f in filenames if f.replace(".h5", "_py27.h5") not in py27_files]
    return filenames

def update_test_files():
    for device in ['APS1', 'APS2']:
        testdirs = glob.glob(os.path.join(BASE_TEST_DIR, 'Test'+device, '*'))
        for test in testdirs:
            testfiles = glob.glob(os.path.join(test, '*'))
            # recurse into subdirectories
            while len(testfiles) == 1 and os.path.isdir(testfiles[0]):
                testfiles = glob.glob(os.path.join(testfiles[0], '*'))
            for tfile in testfiles:
                FID = h5py.File(tfile)
                FID['/'].attrs['target hardware'] = device
                FID.close()

if __name__ == '__main__':
    # run the following line if you are comparing to older h5 files that don't
    # have the 'target hardware' attribute
    # update_test_files()
    output_file()
    compare_sequences()
