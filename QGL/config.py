#Package configuration information

import os.path
import re
import yaml
from pony.orm import *

# Here is a placeholder db
db             = Database()

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

class LoaderMeta(type):
    def __new__(metacls, __name__, __bases__, __dict__):
        """Add include constructer to class."""
        # register the include constructor on the class
        cls = super().__new__(metacls, __name__, __bases__, __dict__)
        cls.add_constructor('!include', cls.construct_include)
        return cls
class Loader(yaml.Loader, metaclass=LoaderMeta):
    """YAML Loader with `!include` constructor."""
    def __init__(self, stream):
        """Initialise Loader."""
        try:
            self._root = os.path.split(stream.name)[0]
        except AttributeError:
            self._root = os.path.curdir
        super().__init__(stream)
        self.add_implicit_resolver(
            u'tag:yaml.org,2002:float',
            re.compile(u'''^(?:
             [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
            |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
            |\\.[0-9_]+(?:[eE][-+][0-9]+)?
            |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
            |[-+]?\\.(?:inf|Inf|INF)
            |\\.(?:nan|NaN|NAN))$''', re.X),
            list(u'-+0123456789.'))
        self.filenames = [os.path.abspath(stream.name)]
    def construct_include(self, node):
        """Include file referenced at node."""
        filename = os.path.abspath(os.path.join(
            self._root, self.construct_scalar(node)
        ))
        extension = os.path.splitext(filename)[1].lstrip('.')
        self.filenames.append(filename)
        with open(filename, 'r') as f:
            if extension in ('yaml', 'yml'):
                return yaml.load(f, Loader)
            else:
                return ''.join(f.readlines())

def find_meas_file():
    if os.getenv('BBN_MEAS_FILE'):
        return os.getenv('BBN_MEAS_FILE')
    raise Exception("Could not find the measurement file in the environment variables or the auspex globals.")

def load_config(filename=None):
    global meas_file, AWGDir, plotBackground, gridColor, pulse_primitives_lib, cnot_implementation

    if filename:
        meas_file = filename
    else:
        meas_file = find_meas_file()

    with open(meas_file, 'r') as FID:
        # cfg = yaml.load(f)
        loader = Loader(FID)
        try:
            cfg = loader.get_single_data()
        finally:
            loader.dispose()

    # pull out the variables
    # abspath allows the use of relative file names in the config file
    if 'AWGDir' in cfg['config'].keys():
        AWGDir = os.path.abspath(cfg['config']['AWGDir'])
    else:
        raise KeyError("Could not find AWGDir in the YAML config section")

    plotBackground       = cfg['config'].get('PlotBackground', '#EAEAF2')
    gridColor            = cfg['config'].get('GridColor', None)
    pulse_primitives_lib = cfg['config'].get('PulsePrimitivesLibrary', 'standard')
    cnot_implementation  = cfg['config'].get('cnot_implementation', 'CNOT_simple')

    return meas_file