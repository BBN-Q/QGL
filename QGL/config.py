#Package configuration information

import os.path
import sys
import re
import importlib

# Where to store AWG data
AWGDir           = None

# The db file, where the channel libraries are stored
db_resource_name = None

# The config file (executed upon channel library loading)
# config_file      = None

# plotting options
plotBackground = '#EAEAF2'
gridColor      = None

# select pulse library (standard or all90)
pulse_primitives_lib = "standard"

# select a CNOT implementation (a name of a Pulse function that implements
# CNOT in your gate set, e.g. CNOT_simple or CNOT_CR)
cnot_implementation  = "CNOT_simple"

def load_config():
    global config_file
    if os.getenv('BBN_CONFIG'):
        try:
            config_file = os.getenv("BBN_CONFIG")
            sys.path.append(os.path.dirname(config_file))
            importlib.import_module(os.path.splitext(os.path.basename(config_file))[0])
        except:
            raise Exception(f"Could not import/execute the BBN_CONFIG {os.getenv('BBN_CONFIG')}")
            
def load_db():
    global db_resource_name
    if os.getenv('BBN_DB'):
        db_resource_name = os.getenv("BBN_DB")
    return db_resource_name