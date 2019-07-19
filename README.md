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
or Anaconda Command Prompt on Windows. While QGL can be run on windows, our 
experiment control software [Auspex](https://github.com/BBN-Q/auspex) relies on linux
when running qubit experiments.

### Python 3.6+

```bash
cd QGL/
pip install .
```
Alternatively, if you plan to modify the source code it will be easier to perform a
developer install using:
```bash
pip install -e .
```
If you'd like to use some of the built-in gate-set-tomography functionality,
you can grab the PyGSTi package during the install:
```
pip install '.[gst]'
```
If the `QGL` module is not installed, the repository path needs to be in the
`PYTHONPATH`. On Windows machines, you add/modify this environment variable by
going to System -> Advanced Settings -> Environment variables. On Mac/Linux
machines add the following line to your .bashrc or .bash_profile: ``` export
PYTHONPATH=/path/to/QGL/repo:$PYTHONPATH```

## Usage
QGL is designed to be utilized alongside Auspex, and most of the usage example, 
including how to define a channel library, can be found in the [Auspex documentation](https://auspex.readthedocs.io/en/develop/qubits.html)

## Dependencies
* Python 3.6+
* Numpy/Scipy
* networkx 2.0
* iPython/Jupyter 4.0 (only for Jupyter notebooks)
* bbndb
