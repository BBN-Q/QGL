#Package configuration information

import json
import os.path

import qgl_config_loc

# Load the configuration from the json file
# and populate the global configuration dictionary

QGLCfgFile = qgl_config_loc.get_config_path()
print('Note: using QGLCfgFile [%s]' % QGLCfgFile)

if not os.path.isfile(QGLCfgFile):
    rootFolder = os.path.dirname(os.path.abspath(__file__))
    rootFolder = rootFolder.replace('\\', '/')  # use unix-like convention

    # build a config file from the template
    templateFile = os.path.join(rootFolder, 'config.example.json')
    ifid = open(templateFile, 'r')
    ofid = open(QGLCfgFile, 'w')
    for line in ifid:
        ofid.write(line.replace('/my/path/to', rootFolder))
    ifid.close()
    ofid.close()

    print('Note: created QGLCfgFile using template [%s]' % templateFile)

with open(QGLCfgFile, 'r') as f:
    QGLCfg = json.load(f)

#pull out the variables
#abspath allows the use of relative file names in the config file
AWGDir = os.path.abspath(QGLCfg['AWGDir'])
configFile = os.path.abspath(QGLCfg['ConfigurationFile'])

# plotting options
plotBackground = QGLCfg.get('PlotBackground', '#EAEAF2')
gridColor = QGLCfg.get('GridColor', None)

# select pulse library (standard or all90)
pulse_primitives_lib = QGLCfg.get('PulsePrimitivesLibrary', 'standard')

# select a CNOT implementation (a name of a Pulse function that implements
# CNOT in your gate set, e.g. CNOT_simple or CNOT_CR)
cnot_implementation = QGLCfg.get('cnot_implementation', 'CNOT_simple')
