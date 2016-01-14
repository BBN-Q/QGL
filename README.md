# Quantum Gate Language

[![Build Status](https://travis-ci.org/BBN-Q/QGL.svg?branch=develop)](https://travis-ci.org/BBN-Q/QGL) [![Coverage Status](https://coveralls.io/repos/BBN-Q/QGL/badge.svg?branch=develop)](https://coveralls.io/r/BBN-Q/QGL)

Quantum Gate Language (QGL) is a domain specific language embedded in python for
specifying pulse sequences.

See example usage in this [Jupyter notebook](https://github.com/BBN-Q/PyQLab/blob/develop/doc/QGL-demo.ipynb).

## Setup instructions

The most straightforward way to get up and running is to use the [Anaconda
Python distribution](http://continuum.io/downloads). This includes nearly all
the dependencies. The few remaining can be installed from the terminal or
Anaconda Command Prompt on Windows

```bash
git clone https://github.com/BBN-Q/JSONLibraryMigrators
pip install watchdog
```

Use of the `QGL` module requires modification of your `PYTHONPATH`. On windows machines, you add/modify this environment variable by going to System -> Advanced Settings -> Environment variables. On Mac/Linux machines add the following line to your .bashrc or .bash_profile:
```
export PYTHONPATH=/path/to/QGL:/path/to/JSONLibraryMigrators:$PYTHONPATH
```

The QGL config file will be created the first time you run `import QGL` or `from QGL import *`.

## Dependencies
* Python 2.7
* JSONLibraryMigrators (https://github.com/BBN-Q/JSONLibraryMigrators)
* Numpy/Scipy
* Nucleic atom
* h5py
* watchdog
* Bokeh 0.7
* iPython 3.0 (only for Jupyter notebooks)
