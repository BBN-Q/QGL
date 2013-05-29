'''
Created on Jan 8, 2012

@author: cryan

A simple GUI for plotting pulse sequences for visual inspection

Copyright 2013 Raytheon BBN Technologies

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

import sys

import matplotlib
matplotlib.rcParams['backend.qt4']='PySide'

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

from PySide import QtCore, QtGui

import numpy as np

from Libraries import instrumentLib
from instruments.AWGs import APS, Tek5014, Tek7000
from APSPattern import read_APS_file
from TekPattern import read_Tek_file

import argparse
import os.path

class PulseSeqPlotWindow(QtGui.QWidget):
    
    def __init__(self, AWGWFs=None):
        super(PulseSeqPlotWindow, self).__init__()
        
        self.AWGWFs = AWGWFs
        
        numSeqs = max([len(tmpChan) for awg in self.AWGWFs.values() for tmpChan in awg.values()])
        
        #Create the GUI
        self.resize(1000,700)
        self.center()
        self.setWindowTitle('Pulse Sequence Plotter')
        
        # generate the plot
        self.fig = Figure(figsize=(12,6), dpi=72)
        self.ax = self.fig.add_subplot(111)

        # generate the canvas to display the plot
        self.canvas = FigureCanvas(self.fig)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self)
        
        #Create a slider to move through the different sequences
        slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        slider.setRange(0, numSeqs-1)
        slider.setTickInterval(1)
        slider.setSingleStep(1)
        sliderLabel = QtGui.QLabel('Sequence Num.')
        if numSeqs > 1:
            numDigits = int(np.ceil(np.log10(numSeqs-1)))
        else:
            numDigits = 1
        sliderLCD = QtGui.QLCDNumber(numDigits)
        slider.valueChanged.connect(sliderLCD.display)
        slider.valueChanged.connect(self.update_plot)
        self.slider = slider
        
        #A tree view to decide what to plot
        plotCheckBoxes = QtGui.QTreeWidget()
        plotCheckBoxes.setMinimumWidth(150)
        plotCheckBoxes.setHeaderLabel('Channel Plot')
        plotCheckBoxes.itemClicked.connect(self.new_plot)
        for awgName, awg in self.AWGWFs.items():
            item = QtGui.QTreeWidgetItem([awgName])
            for chanName in sorted(awg.keys()):
                childItem = QtGui.QTreeWidgetItem([chanName])
                #Default to not checked if there is nothing going on
                if all([np.all(seq == 0) for seq in self.AWGWFs[awgName][chanName]]):
                    childItem.setCheckState(0, QtCore.Qt.Unchecked)
                else:
                    childItem.setCheckState(0,QtCore.Qt.Checked)
                item.addChild(childItem)
            plotCheckBoxes.addTopLevelItem(item)
            item.setExpanded(True)
        self.plotCheckBoxes = plotCheckBoxes
        
        
        delayShiftCheckBox = QtGui.QCheckBox('Undo Delays')
        delayShiftCheckBox.stateChanged.connect(self.update_plot)
        delayShiftCheckBox.setEnabled(False)
        
        #Lay everything out
        hboxSlider = QtGui.QHBoxLayout()
        hboxSlider.addStretch(1)
        hboxSlider.addWidget(sliderLabel)
        hboxSlider.addWidget(sliderLCD)
        hboxSlider.addWidget(slider)

        vboxOptions = QtGui.QVBoxLayout()
        vboxOptions.addWidget(plotCheckBoxes)
        vboxOptions.addWidget(delayShiftCheckBox)
        vboxOptions.addStretch(1)

        hboxTop = QtGui.QHBoxLayout()
        hboxTop.addWidget(self.canvas)
        hboxTop.addLayout(vboxOptions)
        
        hboxBottom = QtGui.QHBoxLayout()
        hboxBottom.addWidget(self.mpl_toolbar)
        hboxBottom.addStretch(1)
        hboxBottom.addLayout(hboxSlider)
        
        vboxTot = QtGui.QVBoxLayout()
        vboxTot.addLayout(hboxTop)
        vboxTot.addLayout(hboxBottom)
        
        self.setLayout(vboxTot) 
        self.lines = {}
        self.new_plot()


    def new_plot(self):
        '''
        Starts a new plot if the channels have changed
        '''
        self.ax.clear()

        #Get the current segment number
        curSegNum = self.slider.sliderPosition()
        
        self.lines = {}
        vertShift = 0
        for itemct in range(self.plotCheckBoxes.topLevelItemCount()):
            tmpItem = self.plotCheckBoxes.topLevelItem(itemct)
            awgName = str(tmpItem.text(0))
            for childct in range(tmpItem.childCount()):
                child = tmpItem.child(childct)
                if child.checkState(0) == QtCore.Qt.Checked:
                    if curSegNum < len(self.AWGWFs[awgName][str(child.text(0))]):
                        chanStr = ''.join([awgName,str(child.text(0))]) 
                        self.lines[chanStr] = {}
                        self.lines[chanStr]['handle'], = self.ax.plot(self.AWGWFs[awgName][str(child.text(0))][curSegNum] + vertShift)
                        self.lines[chanStr]['vertShift'] = vertShift
                        self.lines[chanStr]['awgData'] = self.AWGWFs[awgName][str(child.text(0))]
                    self.ax.text(0, vertShift,awgName+'-'+str(child.text(0)), fontsize=8)
                    vertShift += 2
        self.ax.set_ylim((-1, vertShift))                    
        self.canvas.draw()
        

    def update_plot(self):
        '''
        Just updates the y-date for the sequence number changes
        '''
        #Get the current segment number
        curSegNum = self.slider.sliderPosition()

        for line in self.lines.values():
            if curSegNum < len(line['awgData']):
                line['handle'].set_data(np.arange(line['awgData'][curSegNum].size), line['awgData'][curSegNum] + line['vertShift'])
        
        self.canvas.draw()
        
    def center(self):
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())         
        

def plot_pulse_seqs(AWGWFs):
    '''
    Helper function to plot direct awg waveforms.
    Expects a dictionary keyed on AWG's.
    '''
    
    #Look to see if iPython's event loop is running
    app = QtCore.QCoreApplication.instance()
    if app is None:
        app = QtGui.QApplication(sys.argv)

    plotterWindow = PulseSeqPlotWindow(AWGWFs)
    plotterWindow.show()

    try: 
        from IPython.lib.guisupport import start_event_loop_qt4
        start_event_loop_qt4(app)
    except ImportError:
        sys.exit(app.exec_())
       
    #Need to a keep a reference to the window alive.
    return plotterWindow

def plot_pulse_files(AWGFileNames):
    '''
    Helper function to plot AWG files
    '''
    assert isinstance(AWGFileNames, list), 'Please pass in a list of filenames.'
    AWGWFs = {}
    #Load each of the files using helper functions from APS/TekPattern
    for tmpFile in AWGFileNames:
        #Assume a naming convenction path/to/file/SequenceName-AWGName.h5
        AWGName = (os.path.split(os.path.splitext(tmpFile)[0])[1]).split('-')[1]
        try:
            AWGName = AWGName[:AWGName.index('_')]
        except ValueError:
            pass
        #Look up the appropriate model in the instrumentLib
        if isinstance(instrumentLib[AWGName], Tek5014) or isinstance(instrumentLib[AWGName], Tek7000):
            AWGWFs[AWGName] = read_Tek_file(tmpFile)
        elif isinstance(instrumentLib[AWGName], APS):
            AWGWFs[AWGName] = read_APS_file(tmpFile)
        else:
            raise NameError('Unknown AWG Type for {0}: we currently only handle TekAWG and BBNAPS.'.format(tmpFile))

    #Need to a keep a reference to the window alive.
    return plot_pulse_seqs(AWGWFs)

if __name__ == '__main__':
    
    #See if we have been passed AWG files
    parser = argparse.ArgumentParser()
    parser.add_argument('--AWGFiles', action='store', dest='AWGFiles',  nargs='+', default=None)    
    options =  parser.parse_args(sys.argv[1:])
    if options.AWGFiles:
        plot_pulse_files(options.AWGFiles)
        
    else:
        AWGWFs = {}
        AWGWFs['TekAWG1'] = {}
        AWGWFs['TekAWG1']['ch1'] = []
        AWGWFs['TekAWG1']['ch2'] = []
        AWGWFs['TekAWG1']['ch1m1'] = []
        
        for ct in range(1,10):
            AWGWFs['TekAWG1']['ch1'].append(np.sin(np.linspace(0,ct*np.pi,10000)))
            AWGWFs['TekAWG1']['ch2'].append(np.cos(np.linspace(0,ct*np.pi,10000)))
            AWGWFs['TekAWG1']['ch1m1'].append(np.abs(AWGWFs['TekAWG1']['ch1'][-1]) < 0.5)
    
        plot_pulse_seqs(AWGWFs)