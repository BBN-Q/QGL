'''
Channels is where we store information for mapping virtual (qubit) channel to real channels.

Created on Jan 19, 2012

@author: cryan


'''

import sys
import json
import PulseShapes

from PySide import QtGui, QtCore
from operator import itemgetter
from math import tan,cos,pi

class ChannelTypes(object):
    '''
    Enumerate the possible types:
    direct - goes straight to something i.e. no modulated carrier
    digital - a logical digital channel usually assigned to a marker channel
    amplitudeMod - an amplitude modulated carrier
    quadratureMod - a quadrature modulated carrier
    '''
    (direct, marker, amplitudeMod, quadratureMod) = range(4)
    (logical, physical) = range(2)
    
class LogicalChannel(object):
    '''
    The main class from which we will generate sequences.  At some point it needs to be assigned to a physical channel.
    '''
    def __init__(self, name=None, channelType=None, physicalChannel=None):
        self.name = name
        self.channelType = channelType
        self.physicalChannel = physicalChannel
        
    @property
    def isLogical():
        return True
        
    @property
    def isPhysical():
        return False

    # add __eq__ method to make hashable
    def __eq__(self, other):
        return id(self) == id(other)

class PhysicalChannel(object):
    '''
    The main class for actual AWG channels.
    '''
    def __init__(self, name=None, AWGName=None, channelType=None, samplingRate=1.0e9):
        self.name = name
        self.channelType = channelType
        self.AWGName = AWGName
        self.samplingRate = samplingRate

    @property
    def isLogical():
        return False
        
    @property
    def isPhysical():
        return True

    # add __eq__ method to make hashable
    def __eq__(self, other):
        return id(self) == id(other)


class PhysicalMarkerChannel(PhysicalChannel):
    '''
    An digital output channel on an AWG.
    '''
    def __init__(self, name=None, AWGName=None, channel=None, channelShift=0.0, **kwargs):
        super(PhysicalMarkerChannel, self).__init__(name=name, AWGName=AWGName, channelType=ChannelTypes.marker)
        self.channelType = ChannelTypes.marker
        self.channelShift = channelShift
        self.channel = channel
        
class QuadratureChannel(PhysicalChannel):
    '''
    Something used to implement a standard qubit channel with two analog channels and a microwave gating channel.
    '''
    def __init__(self, name=None, AWGName=None, carrierGen=None, IChannel=None, QChannel=None, channelShift=0.0, ampFactor=1.0, phaseSkew=0.0, **kwargs):
        super(QuadratureChannel, self).__init__(name=name, AWGName=AWGName, channelType=ChannelTypes.quadratureMod)
        self.carrierGen = carrierGen
        self.IChannel = IChannel
        self.QChannel = QChannel
        self.channelShift = channelShift
        self.ampFactor = ampFactor
        self.phaseSkew = phaseSkew

    @property
    def correctionT(self):
        return [[self.ampFactor, self.ampFactor*tan(self.phaseSkew*pi/180)], [0, 1/cos(self.phaseSkew*pi/180)]]
                
class LogicalMarkerChannel(LogicalChannel):
    '''
    A class for digital channels for gating sources or triggering other things.
    '''
    def __init__(self, name=None, physicalChannel=None, **kwargs):
        super(LogicalMarkerChannel, self).__init__(name=name, channelType=ChannelTypes.marker, physicalChannel=physicalChannel)        
    
    def gatePulse(self, length, delay=0):
        tmpBlock = PulseSequencer.PulseBlock()
        if delay>0:
            tmpBlock.add_pulse(PatternGen.QId(delay), self)
        tmpBlock.add_pulse(PatternGen.Square(length, amp=1), self)
        return tmpBlock
        
class Qubit(LogicalChannel):
    '''
    The main class for generating qubit pulses.  
    '''
    def __init__(self, name=None, physicalChannel=PhysicalChannel(), freq=None, piAmp=0.0, pi2Amp=0.0, shapeFun=PulseShapes.gaussian, pulseLength=0.0, bufferTime=0.0, dragScaling=0, cutoff=2, **kwargs):
        super(Qubit, self).__init__(name=name, channelType=ChannelTypes.quadratureMod, physicalChannel=physicalChannel)
        self.shapeFun = shapeFun
        self.pulseLength = pulseLength
        self.bufferTime = bufferTime
        self.piAmp = piAmp
        self.pi2Amp = pi2Amp
        self.dragScaling = dragScaling
        self.cutoff = cutoff        
        
class Generator(object):
    '''
    Although not quite a channel, it is tightly linked to channels.
    '''
    def __init__(self, name=None, gateChannel=None, gateBuffer=0.0, gateMinWidth=0.0, gateChannelShift=0.0):
        self.name = name
        self.gateChannel = gateChannel
        self.gateBuffer = gateBuffer
        self.gateMinWidth = gateMinWidth
        self.gateChannelShift = gateChannelShift
        
        
def save_channel_info(channelDict, fileName=None):
    '''
    Helper function to save a channelInfo dictionary to a JSON file or string.
    '''
    #Convert the channel into a dictionary
    if fileName is None:
        return json.dumps(channelDict, sort_keys=True, indent=2)
    else:
        with open(fileName,'w') as FID:
            json.dump(channelDict, FID, sort_keys=True, indent=2)
    
def load_channel_dict(fileName=None):
    '''
    Helper function to convert back from a JSON file to a channel dictionary
    '''
    with open(fileName,'r') as FID:
        return json.load(FID)
        
def load_channel_info(fileName=None):
    '''
    Helper function to load a channel info file into channel and channelInfo objects
    '''
    channels = {}
    #First load the file into a dictionary
    channelDicts = load_channel_dict(fileName)
    for channelName, channelDict in channelDicts.items():
        channelType = channelDict['channelType']
        #Deal with logical channels
        if channelDict['isLogical']:
            if channelType == 'quadratureMod':
                #Create the qubit channel
                channelFunc = QubitChannel
            elif channelType == 'marker':
                #Create the marker channel
                channelFunc = LogicalMarkerChannel
        else:
            if channelType == 'quadratureMod':
                channelFunc = QuadratureChannel
            elif channelType == 'marker':
                channelFunc = PhysicalMarkerChannel

        channels[channelName] = channelFunc(**channelDict)

    return channels, channelDicts
    
'''    
*****************************************************************************
GUI Stuff.
*****************************************************************************
'''
class ChannelInfoView(QtGui.QMainWindow):
    def __init__(self, fileName):
        super(ChannelInfoView, self).__init__()
        
        #Load the channel information from the file
        self.fileName = fileName
        self.channelDict = load_channel_dict(fileName)
        
        #Create an item view for the logical channels
        self.logicalChannelListModel = QtGui.QStringListModel([tmpKey for tmpKey in self.channelDict.keys() if self.channelDict[tmpKey]['isLogical']])
        self.logicalChannelListModel.sort(0)
        self.logicalChannelListView = QtGui.QListView()
        self.logicalChannelListView.setModel(self.logicalChannelListModel)
        self.logicalChannelListView.clicked.connect(lambda(idx): self.update_channelView(idx, self.logicalChannelListModel))
        
        #Create an item view for the physical channels
        self.physicalChannelListModel = QtGui.QStringListModel([tmpKey for tmpKey in self.channelDict.keys() if self.channelDict[tmpKey]['isPhysical']])
        self.physicalChannelListModel.sort(0)
        self.physicalChannelListView = QtGui.QListView()
        self.physicalChannelListView.setModel(self.physicalChannelListModel)
        self.physicalChannelListView.clicked.connect(lambda(idx): self.update_channelView(idx, self.physicalChannelListModel))

        #Create an item view for the physical channels
        self.generatorListModel = QtGui.QStringListModel([tmpKey for tmpKey in self.channelDict.keys() if self.channelDict[tmpKey]['isGenerator']])
        self.generatorListModel.sort(0)
        self.generatorListView = QtGui.QListView()
        self.generatorListView.setModel(self.generatorListModel)
        self.generatorListView.clicked.connect(lambda(idx): self.update_channelView(idx, self.generatorListModel))

        tmpWidget = QtGui.QWidget()
        vBox = QtGui.QVBoxLayout(tmpWidget)
        
        self.tabWidget = QtGui.QTabWidget()
        self.tabWidget.addTab(self.logicalChannelListView, 'Logical')        
        self.tabWidget.addTab(self.physicalChannelListView, 'Physical')        
        self.tabWidget.addTab(self.generatorListView, 'Generators')        
        vBox.addWidget(self.tabWidget)

        #Add the buttons for adding/deleting channels
        hBox = QtGui.QHBoxLayout()
        addChanButton = QtGui.QPushButton('Add')
        addChanButton.clicked.connect(self.add_channel)
        hBox.addWidget(addChanButton)
        deleteChanButton = QtGui.QPushButton('Delete')
        deleteChanButton.clicked.connect(self.delete_channel)
        hBox.addWidget(deleteChanButton)
        hBox.addStretch(1)
        vBox.addLayout(hBox)                

        #Setup the main splitter between the channel lists and channel properties
        hSplitter = QtGui.QSplitter()
        hSplitter.addWidget(tmpWidget)
        self.channelWidgets = {}
        for tmpChanName, tmpChan in self.channelDict.items():
            self.channelWidgets[tmpChanName] = ChannelView(tmpChan)
            hSplitter.addWidget(self.channelWidgets[tmpChanName])
            self.channelWidgets[tmpChanName].hide()
        
        self.setCentralWidget(hSplitter)

        #Setup the toolbar buttons for loading and saving files
        loadAction = QtGui.QAction('Load', self)
        loadAction.setShortcut('Ctrl+L')
        loadAction.setStatusTip('Load parameter file')
        
        saveAction = QtGui.QAction('Save',self)
        saveAction.setShortcut('Ctrl+S')
        saveAction.setStatusTip('Save parameter file')
        saveAction.triggered.connect(self.save_to_file)
        
        saveAsAction = QtGui.QAction('Save As',self)
        saveAsAction.setStatusTip('Save parameter file to new file')

        exitAction = QtGui.QAction('Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        
        #Setup the menus
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(exitAction)
        
        #Setup the toolbar
        self.toolbar = self.addToolBar('Exit')
        self.toolbar.addAction(saveAction)
        self.toolbar.addAction(saveAsAction)
        self.toolbar.addAction(loadAction)
        self.toolbar.addAction(exitAction)

        self.statusBar()
        
        self.setGeometry(300,300,550,300)
        self.setWindowTitle('Channel Info.')
        

    def update_channelView(self, index, model):
        for tmpWidget in self.channelWidgets.values():
            tmpWidget.hide()
        
        tmpChan = model.data(index, QtCore.Qt.DisplayRole)
        self.channelWidgets[tmpChan].show()
        for tabct in range(3):
            if tabct != self.tabWidget.currentIndex():
                self.tabWidget.widget(tabct).clearSelection()
        
    def save_to_file(self):
        #Update the dictionary from the GUI fields
        for tmpChanName, tmpChan in self.channelDict.items():
            self.channelWidgets[tmpChanName].update_from_view()
        save_channel_info(self.channelDict, self.fileName)  
        
    def add_channel(self):
        msgBox = QtGui.QMessageBox()
        msgBox.setText("Not implemented yet!")
        msgBox.exec_()
        
    def delete_channel(self):
        msgBox = QtGui.QMessageBox()
        msgBox.setText("Not implemented yet!")
        msgBox.exec_()
        
        
class ChannelView(QtGui.QWidget):
    def __init__(self, channel):
        super(ChannelView, self).__init__()
        self.channel = channel
        
        skipFields = ['channelType', 'name', 'isLogical', 'isPhysical', 'isGenerator', 'correctionT']

        #Create the layout as a vbox of hboxes
        form = QtGui.QFormLayout()
        self.GUIhandles = {}
        #Do the channelType on top
        form.addRow('channelType', QtGui.QLabel(channel['channelType']))
        
        for key,value in sorted(channel.items(), key=itemgetter(0)):
            if key not in skipFields:
                if isinstance(value, basestring):
                    tmpWidget = QtGui.QLineEdit(value)
                else:
                    tmpWidget = QtGui.QLineEdit(str(value))
                    tmpWidget.setValidator(QtGui.QDoubleValidator())
                form.addRow(key, tmpWidget)
                self.GUIhandles[key] = tmpWidget
        self.setLayout(form)
            
    def update_from_view(self):
        '''
        Update the channel dictionary with the current values entered
        '''
        for key,tmpWidget in self.GUIhandles.items():
            if isinstance(self.channel[key], basestring):
                self.channel[key] = tmpWidget.text()
            else:
                self.channel[key] = float(tmpWidget.text())

        #Calculate the correction T if necessary
        if self.channel['channelType'] == 'quadratureMod' and self.channel['isPhysical']:
            self.channel['correctionT'] = [[self.channel['ampFactor'], self.channel['ampFactor']*tan(self.channel['phaseSkew']*pi/180)], [0, 1/cos(self.channel['phaseSkew']*pi/180)]]
    
    

if __name__ == '__main__':
    channelDict = {}
    channelDict['q1'] = {'name':'q1', 'channelType':'quadratureMod', 'isLogical':True, 'isPhysical':False, 'isGenerator':False, 'piAmp':1.0, 'pi2Amp':0.5, 'pulseType':'drag', 'pulseLength':40e-9, 'bufferTime':2e-9, 'dragScaling':1, 'physicalChannel':'TekAWG2-12', 'frequency':5}
    channelDict['q2'] = {'name':'q2', 'channelType':'quadratureMod', 'isLogical':True, 'isPhysical':False, 'isGenerator':False, 'piAmp':1.0, 'pi2Amp':0.5, 'pulseType':'drag', 'pulseLength':40e-9, 'bufferTime':2e-9, 'dragScaling':1, 'physicalChannel':'TekAWG2-34', 'frequency':5}
    channelDict['CR'] = {'name':'CR', 'channelType':'quadratureMod', 'isLogical':True, 'isPhysical':False, 'isGenerator':False, 'piAmp':1.0, 'pi2Amp':0.5, 'pulseType':'drag', 'pulseLength':40e-9, 'bufferTime':2e-9, 'dragScaling':1, 'physicalChannel':'TekAWG1-12', 'frequency':5}

    channelDict['measChannel'] = {'name':'measChannel', 'channelType':'marker', 'isLogical':True, 'isPhysical':False, 'isGenerator':False, 'physicalChannel':'TekAWG1-ch3m1' }
    channelDict['digitizerTrig'] = {'name':'digitizerTrig','channelType':'marker', 'isLogical':True, 'isPhysical':False, 'isGenerator':False, 'physicalChannel':'TekAWG1-ch3m2'}

    channelDict['TekAWG1-12'] = {'name':'TekAWG1-12', 'channelType':'quadratureMod', 'isLogical':False, 'isPhysical':True, 'isGenerator':False, 'AWGName':'TekAWG1', 'IChannel':'ch1', 'QChannel':'ch2', 'channelShift':0e-9,  'correctionT':[[1,0],[0,1]], 'ampFactor':1.0, 'phaseSkew':0.0, 'carrierGen':'QPC1-1691'}
    channelDict['TekAWG1-34'] = {'name':'TekAWG1-34', 'channelType':'quadratureMod', 'isLogical':False, 'isPhysical':True, 'isGenerator':False, 'AWGName':'TekAWG1', 'IChannel':'ch3', 'QChannel':'ch4', 'channelShift':0e-9,  'correctionT':[[1,0],[0,1]], 'ampFactor':1.0, 'phaseSkew':0.0, 'carrierGen':'Agilent1'}
    channelDict['TekAWG1-ch1m1'] = {'name':'TekAWG1-ch1m1', 'channelType':'marker', 'isLogical':False, 'isPhysical':True, 'isGenerator':False, 'AWGName':'TekAWG1', 'channel':'ch1m1', 'channelShift':0e-9 }    
    channelDict['TekAWG1-ch2m1'] = {'name':'TekAWG1-ch2m1', 'channelType':'marker', 'isLogical':False, 'isPhysical':True, 'isGenerator':False, 'AWGName':'TekAWG1', 'channel':'ch2m1', 'channelShift':0e-9 }
    channelDict['TekAWG1-ch3m1'] = {'name':'TekAWG2-ch3m1', 'channelType':'marker', 'isLogical':False, 'isPhysical':True, 'isGenerator':False, 'AWGName':'TekAWG2', 'channel':'ch3m1', 'channelShift':0e-9 }    
    channelDict['TekAWG1-ch3m2'] = {'name':'TekAWG2-ch3m2', 'channelType':'marker', 'isLogical':False, 'isPhysical':True, 'isGenerator':False, 'AWGName':'TekAWG2', 'channel':'ch3m2', 'channelShift':0e-9 }
  
    channelDict['BBNAPS1-12'] = {'name':'BBNAPS1-12', 'channelType':'quadratureMod', 'isLogical':False, 'isPhysical':True, 'isGenerator':False, 'AWGName':'BBNAPS1', 'IChannel':'ch1', 'QChannel':'ch2', 'channelShift':0e-9, 'correctionT':[[1,0],[0,1]], 'ampFactor':1.0, 'phaseSkew':0.0, 'carrierGen':'QPC1-1691'}
    channelDict['BBNAPS1-34'] = {'name':'BBNAPS1-34', 'channelType':'quadratureMod', 'isLogical':False, 'isPhysical':True, 'isGenerator':False, 'AWGName':'BBNAPS1', 'IChannel':'ch3', 'QChannel':'ch4', 'channelShift':0e-9, 'correctionT':[[1,0],[0,1]], 'ampFactor':1.0, 'phaseSkew':0.0, 'carrierGen':'Agilent1'}

    channelDict['TekAWG2-12'] = {'name':'TekAWG2-12', 'channelType':'quadratureMod', 'isLogical':False, 'isPhysical':True, 'isGenerator':False, 'AWGName':'TekAWG2', 'IChannel':'ch1', 'QChannel':'ch2', 'channelShift':0e-9,  'correctionT':[[1,0],[0,1]], 'ampFactor':1.0, 'phaseSkew':0.0, 'carrierGen':'QPC1-1691'}
    channelDict['TekAWG2-34'] = {'name':'TekAWG2-34', 'channelType':'quadratureMod', 'isLogical':False, 'isPhysical':True, 'isGenerator':False, 'AWGName':'TekAWG2', 'IChannel':'ch3', 'QChannel':'ch4', 'channelShift':0e-9,  'correctionT':[[1,0],[0,1]], 'ampFactor':1.0, 'phaseSkew':0.0, 'carrierGen':'Agilent1'}


    channelDict['QPC1-1691'] = {'name':'QPC1-1691', 'channelType':'generator', 'isLogical':False, 'isPhysical':False, 'isGenerator':True, 'gateChannel':'TekAWG1-ch1m1', 'gateChannelShift':-50.0e-9, 'gateBuffer':20e-9, 'gateMinWidth':100e-9, 'frequency':5}    
    channelDict['Agilent1'] = {'name':'Agilent1', 'channelType':'generator', 'isLogical':False, 'isPhysical':False, 'isGenerator':True, 'gateChannel':'TekAWG1-ch2m1', 'gateChannelShift':0.0, 'gateBuffer':20e-9, 'gateMinWidth':100e-9, 'frequency':5}    

    save_channel_info(channelDict, 'ChannelParams.json')

    #Look to see if iPython's event loop is running
    app = QtCore.QCoreApplication.instance()
    if app is None:
        app = QtGui.QApplication(sys.argv)

    channelWindow = ChannelInfoView('ChannelParams.json')
    channelWindow.show()

    try: 
        from IPython.lib.guisupport import start_event_loop_qt4
        start_event_loop_qt4(app)
    except ImportError:
        sys.exit(app.exec_())


    
    