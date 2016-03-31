# Quantum Gate Language

[![Build Status](https://travis-ci.org/BBN-Q/QGL.svg?branch=master)](https://travis-ci.org/BBN-Q/QGL) [![Coverage Status](https://coveralls.io/repos/BBN-Q/QGL/badge.svg?branch=master)](https://coveralls.io/r/BBN-Q/QGL)

Quantum Gate Language (QGL) is a domain specific language embedded in python for
specifying pulse sequences.

See example usage in this [Jupyter notebook](https://github.com/BBN-Q/PyQLab/blob/develop/doc/QGL-demo.ipynb).

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
pip install cppy
pip install git+https://github.com/nucleic/atom.git@1.0.0-dev
```

Use of the `QGL` module requires modification of your `PYTHONPATH`. On windows machines, you add/modify this environment variable by going to System -> Advanced Settings -> Environment variables. On Mac/Linux machines add the following line to your .bashrc or .bash_profile:
```
export PYTHONPATH=/path/to/QGL:$PYTHONPATH
```

The QGL config file will be created the first time you run `import QGL` or `from QGL import *`.

## Dependencies
* Python 2.7 or 3.4+
* JSONLibraryUtils (https://github.com/BBN-Q/JSONLibraryUtils, integrated as a Git submodule)
* Numpy/Scipy
* Nucleic atom (1.0 branch for Python 3)
* h5py
* watchdog
* Bokeh 0.11
* networkx
* future
* iPython/Jupyter 4.0 (only for Jupyter notebooks)
