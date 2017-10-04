# Quantum Gate Language

[![Build Status](https://travis-ci.org/BBN-Q/QGL.svg?branch=master)](https://travis-ci.org/BBN-Q/QGL) [![Coverage Status](https://coveralls.io/repos/BBN-Q/QGL/badge.svg?branch=master)](https://coveralls.io/r/BBN-Q/QGL)

Quantum Gate Language (QGL) is a domain specific language embedded in python for
specifying pulse sequences.

Read the [online documentation](https://bbn-q.github.io/QGL/) and see example
usage in this [Jupyter
notebook](https://github.com/BBN-Q/QGL/blob/master/doc/QGL-demo.ipynb).

## Setup instructions

The most straightforward way to get up and running is to use the [Anaconda
Python distribution](http://continuum.io/downloads). This includes nearly all
the dependencies. The remaining dependencies can be installed from the terminal
or Anaconda Command Prompt on Windows.

### Python 2.7

```bash
conda install atom future
pip install watchdog
```

### Python 3.4+

```bash
conda install future
pip install watchdog
conda install -c ecpy atom
```

For Python 3.6+ you may need to install watchdog from conda forge:
```
conda install -c conda-forge watchdog
```

If the `QGL` module is not installed, the repository path needs to be in the
`PYTHONPATH`. On Windows machines, you add/modify this environment variable by
going to System -> Advanced Settings -> Environment variables. On Mac/Linux
machines add the following line to your .bashrc or .bash_profile: ``` export
PYTHONPATH=/path/to/QGL/repo:$PYTHONPATH```

The QGL config file will be created the first time you run `import QGL` or `from QGL import *`.

## Dependencies
* Python 2.7 or 3.4+
* JSONLibraryUtils (https://github.com/BBN-Q/JSONLibraryUtils, integrated as a Git submodule)
* Numpy/Scipy
* Nucleic atom (from ecpy channel for Python 3)
* h5py
* watchdog
* Bokeh 0.11
* networkx
* iPython/Jupyter 4.0 (only for Jupyter notebooks)
