#Package configuration information

import json
import os.path
import sys

#Load the configuration from the json file and populate the global configuration dictionary
rootFolder = os.path.dirname( os.path.abspath(__file__) )
rootFolder = rootFolder.replace('\\', '/') # use unix-like convention
QGLCfgFile = os.path.join(rootFolder, 'config.json')
if not os.path.isfile(QGLCfgFile):
	# build a config file from the template
	templateFile = os.path.join(rootFolder, 'config.example.json')
	ifid = open(templateFile, 'r')
	ofid = open(QGLCfgFile, 'w')
	for line in ifid:
		ofid.write(line.replace('/my/path/to', rootFolder))
	ifid.close()
	ofid.close()


with open(QGLCfgFile, 'r') as f:
	QGLCfg = json.load(f)

#pull out the variables
#abspath allows the use of relative file names in the config file
AWGDir = os.path.abspath(QGLCfg['AWGDir'])
channelLibFile = os.path.abspath(QGLCfg['ChannelLibraryFile'])

# plotting options
plotBackground = QGLCfg['PlotBackground'] if 'PlotBackground' in QGLCfg else '#EAEAF2'
gridColor = QGLCfg['GridColor'] if 'GridColor' in QGLCfg else None
