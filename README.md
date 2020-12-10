# Quantum Gate Language (QGL)  

![Build Status](https://github.com/gribeill/QGL/workflows/Python%20Package%20using%20Conda/badge.svg?branch=develop)

# Overview

Quantum Gate Language (QGL) is a domain specific language embedded in python for
specifying pulse sequences.

Our "getting started" documentation is published [online](https://gribeill.github.io/QGL/readme/index.html) from the local 
[file](doc/readme.md). This documentation includes dependency, 
installation and basic programming information. The documentation is organized into the following
sections:   

1. What is QGL
1. Dependencies
1. Installation
1. Examples
1. Channels and Qubits
1. Gate Primitives
1. Sequences and Concurrent Operations
1. Pulse Shapes and Waveforms
1. Compiling and Plotting
1. Built-in Basic Sequences

## Usage

There are a number of QGL example Jupyer notebooks in the QGL/doc 
[folder](doc/):

1. ex1_basic_QGL.ipynb: Basic setup of 'qubit' objects, defining sequences of pulses on qubits, and visualizing these pulse sequences.
1. ex2_single_qubit_sequences.ipynb: Simple spectroscopy and coherence experiments on a single qubit.
1. ex3_two_qubit_sequences.ipynb: Examples of two-qubit sequences, including CR gates.

Obviously, we suggest that you start with ex1_basic_QGL.   

QGL requires the installation and use of [bbndb](https://github.com/BBN-Q/bbndb). bbndb is a 
shared, versioned, means of storing instrument, qubit, and other configuration information. 
It is based on the SQLAlchemy framework.

QGL is typically used with Auspex -- an experiment management framework. More sophisticated uses of bbndb, 
especially usage of a channel library, can be found in the 
[Auspex documentation](https://auspex.readthedocs.io/en/develop/qubits.html)

## Dependencies

* Python 3.6+
* [bbndb](https://github.com/BBN-Q/bbndb)

Note additional setup information in [setup.py](setup.py). This file is typically used by pip and other package managers. 
