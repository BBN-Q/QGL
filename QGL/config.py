#Package configuration information

import os.path
import sys
import re
import importlib

# Where to store AWG data
AWGDir         = None

# The measurement file
meas_file      = None

# plotting options
plotBackground = '#EAEAF2'
gridColor      = None

# select pulse library (standard or all90)
pulse_primitives_lib = "standard"

# select a CNOT implementation (a name of a Pulse function that implements
# CNOT in your gate set, e.g. CNOT_simple or CNOT_CR)
cnot_implementation  = "CNOT_simple"

def load_config():
    if os.getenv('BBN_CONFIG_FILE'):
        cfg = os.getenv("BBN_CONFIG_FILE")
        sys.path.append(os.path.dirname(cfg))
        importlib.import_module(os.path.splitext(os.path.basename(cfg))[0])
    else:
        raise Exception("Could not find the measurement file in the environment variables or the auspex globals.")
