'''
Channels is where we store information for mapping virtual (qubit) channel to real channels.

Created on Jan 19, 2012

@author: cryan


'''

import sys
import json
import PulseShapes
import config

from PySide import QtGui, QtCore
from operator import itemgetter
from math import tan,cos,pi
import config

ChannelDict = {}

class Channel(object):
    '''
    Channel superclass
    '''
    def __init__(self, name=None):
        self.name = name

    @property
    def isLogical(self):
        return False
    @property
    def isPhysical(self):
        return False
    @property
    def isGenerator(self):
        return False
    @property
    def isAWG(self):
        return False

    # add __eq__ method to make hashable
    def __eq__(self, other):
        return id(self) == id(other)

class LogicalChannel(Channel):
    '''
    The main class from which we will generate sequences. 
    At some point it needs to be assigned to a physical channel.
    '''
    def __init__(self, name=None, AWGName=None, channelType=None, samplingRate=1.2e9):
        self.name = name
        self.physicalChannel = physicalChannel

    @property
    def isLogical(self):
        return True
    
class PhysicalChannel(Channel):
    '''
    The main class for actual AWG channels.
    '''
    def __init__(self, name=None, AWG=None):
        super(PhysicalChannel, self).__init__(name=name)
        self.AWG = AWG

    @property
    def isPhysical(self):
        return True

    @property
    def samplingRate(self):
        return ChannelDict[self.AWG].samplingRate

    # add __eq__ method to make hashable
    def __eq__(self, other):
        return id(self) == id(other)


class PhysicalMarkerChannel(PhysicalChannel):
    '''
    An digital output channel on an AWG.
    '''
    def __init__(self, name=None, AWG=None, channel=None, delay=0.0, **kwargs):
        super(PhysicalMarkerChannel, self).__init__(name=name, AWG=AWG)
        self.delay = delay
        self.channel = channel
        
class PhysicalQuadratureChannel(PhysicalChannel):
    '''
    Something used to implement a standard qubit channel with two analog channels and a microwave gating channel.
    '''
    def __init__(self, name=None, AWG=None, carrierGen=None, IChannel=None, QChannel=None, delay=0.0, ampFactor=1.0, phaseSkew=0.0, **kwargs):
        super(PhysicalQuadratureChannel, self).__init__(name=name, AWG=AWG)
        self.carrierGen = carrierGen
        self.IChannel = IChannel
        self.QChannel = QChannel
        self.delay = delay
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
        super(LogicalMarkerChannel, self).__init__(name=name, physicalChannel=physicalChannel)        
    
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
    def __init__(self, name=None, physicalChannel=None, freq=None, piAmp=0.0, pi2Amp=0.0, shapeFun=PulseShapes.gaussian, pulseLength=0.0, bufferTime=0.0, dragScaling=0, cutoff=2, **kwargs):
        super(Qubit, self).__init__(name=name, physicalChannel=physicalChannel)
        if name is not None:
            # try to look up data in parameter file
            with open(config.ChannelParams) as f:
                channelParams = json.load(f)
            if name in channelParams.keys():
                qubitParams = channelParams[name]
            else:
                raise NameError("Did not find '{0}' in channel parameter file".format(name))
            self.shapeFun = getattr(PulseShapes, qubitParams['pulseType'])
            self.pulseLength = qubitParams['pulseLength']
            self.bufferTime = qubitParams['bufferTime']
            self.piAmp = qubitParams['piAmp']
            self.pi2Amp = qubitParams['pi2Amp']
            self.dragScaling = qubitParams['dragScaling']
            self.cutoff = qubitParams['cutoff']
        else:
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
    def __init__(self, name=None, gateChannel=None, gateBuffer=0.0, gateMinWidth=0.0, gateDelay=0.0):
        super(Generator, self).__init__(name=name)
        self.gateChannel = gateChannel
        self.gateBuffer = gateBuffer
        self.gateMinWidth = gateMinWidth
        self.gateDelay = gateDelay

    @property
    def isGenerator(self):
        return True

class AWG(Channel):
    '''
    Although not quite a channel, it is tightly linked to channels.
    '''
    def __init__(self, name=None, model=None):
        super(AWG, self).__init__(name=name)
        self.model = model
        
    @property
    def isAWG(self):
        return True

def save_channel_info(fileName=None):
    '''
    Helper function to save a channelInfo dictionary to a JSON file or string.
    '''
    #Convert the channel into a dictionary
    if fileName is None:
        fileName = config.channelParamsFile
    with open(fileName,'w') as FID:
        json.dump(ChannelDict, FID, sort_keys=True, indent=2, default=json_serializer)
    
def load_channel_dict(fileName=None):
    '''
    Helper function to convert back from a JSON file to a channel dictionary
    '''
    with open(fileName,'r') as FID:
        return json.load(FID)
        
def update_channel_info(fileName=None):
    '''
    Helper function to load a channel info file into channel and channelInfo objects
    '''
    global ChannelDict
    if fileName is None:
        fileName = config.channelParamsFile
    with open(fileName,'r') as FID:
        ChannelDict = json.load(FID, object_hook=json_deserializer)

def json_serializer(obj):
    '''
    Helper function to flatten the channel classes to a dictionary for json serialization.
    We just keep the class name and the properties
    '''
    jsonDict = {'__class__': obj.__class__.__name__}
    jsonDict.update(obj.__dict__)
    #Deal with shape function handles specially
    if 'shapeFun' in jsonDict:
        jsonDict['shapeFun'] = jsonDict['shapeFun'].__name__
    return jsonDict

def json_deserializer(jsonDict):
    '''
    Helper function to convert a json representation of a channel back into an object.
    '''
    #Extract the class name from the dictionary
    #If there is no class then assume top level dictionary
    if '__class__' not in jsonDict:
        return jsonDict
    else:
        className = jsonDict.pop('__class__')
        class_ = getattr(sys.modules[__name__], className)
        #Deal with shape functions
        if 'shapeFun' in jsonDict:
            jsonDict['shapeFun'] = getattr(PulseShapes, jsonDict['shapeFun'])
        return class_(**jsonDict)


'''  
*****************************************************************************
GUI Stuff.
*****************************************************************************
'''
class ChannelInfoView(QtGui.QMainWindow):
    def __init__(self):
        super(ChannelInfoView, self).__init__()
        
        global ChannelDict
        
        #Create an item view for the logical channels
        self.logicalChannelListModel = QtGui.QStringListModel([tmpKey for tmpKey, tmpChan in ChannelDict.items() if tmpChan.isLogical])
        self.logicalChannelListModel.sort(0)
        self.logicalChannelListView = QtGui.QListView()
        self.logicalChannelListView.setModel(self.logicalChannelListModel)
        self.logicalChannelListView.clicked.connect(lambda(idx): self.update_channelView(idx, self.logicalChannelListModel))
        
        #Create an item view for the physical channels
        self.physicalChannelListModel = QtGui.QStringListModel([tmpKey for tmpKey, tmpChan in ChannelDict.items() if tmpChan.isPhysical])
        self.physicalChannelListModel.sort(0)
        self.physicalChannelListView = QtGui.QListView()
        self.physicalChannelListView.setModel(self.physicalChannelListModel)
        self.physicalChannelListView.clicked.connect(lambda(idx): self.update_channelView(idx, self.physicalChannelListModel))

        #Create an item view for the physical channels
        self.generatorListModel = QtGui.QStringListModel([tmpKey for tmpKey, tmpChan in ChannelDict.items() if tmpChan.isGenerator])
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
        for chanName, chan in ChannelDict.items():
            self.channelWidgets[chanName] = ChannelView(chan, self)
            hSplitter.addWidget(self.channelWidgets[chanName])
            self.channelWidgets[chanName].hide()
        
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
        print(self.physicalChannelListModel.stringList())
        self.physicalChannelListModel.removeRows(0,1)
        print(self.physicalChannelListModel.stringList())
        msgBox = QtGui.QMessageBox()
        msgBox.setText("Not implemented yet!")
        msgBox.exec_()
        
    def delete_channel(self):
        msgBox = QtGui.QMessageBox()
        msgBox.setText("Not implemented yet!")
        msgBox.exec_()
        
        
class ChannelView(QtGui.QWidget):
    def __init__(self, channel, parent = None):
        super(ChannelView, self).__init__()
        self.channel = channel
        
        #Some hidden fields
        skipFields = ['name']

        #Create the layout as a vbox of hboxes
        form = QtGui.QFormLayout()
        self.GUIhandles = {}
        #Do the channelType on top for information purposes
        form.addRow('channelType', QtGui.QLabel(self.channel.__class__.__name__))

        #Helper function to update         
        for key,value in sorted(channel.__dict__.items(), key=itemgetter(0)):
            if key not in skipFields:
                #For the physical channel field we'll pop up a combo box
                if key == 'physicalChannel':
                    chanType = ChannelDict[value].__class__
                    tmpModel = QtGui.QStringListModel()
                    parent.physicalChannelListModel.rowsRemoved.connect(lambda a,b,c : tmpModel.removeRow(b))
                    tmpModel.setStringList( [tmpChan for tmpChan in parent.physicalChannelListModel.stringList() if ChannelDict[tmpChan].__class__ == chanType] )
                    tmpWidget = QtGui.QComboBox()
                    tmpWidget.setModel(tmpModel)
                elif isinstance(value, basestring):
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

    ChannelDict['q1'] = Qubit(name='q1', piAmp=1.0, pi2Amp=0.5, pulseType='drag', pulseLength=40e-9, bufferTime=2e-9, dragScaling=1, physicalChannel='BBNAPS1-12', carrierGen='QPC1-1691')
    ChannelDict['q2'] = Qubit(name='q2', piAmp=1.0, pi2Amp=0.5, pulseType='drag', pulseLength=40e-9, bufferTime=2e-9, dragScaling=1, physicalChannel='BBNAPS1-34', carrierGen='Agilent1')

    ChannelDict['digitizerTrig'] = LogicalMarkerChannel(name='digitizerTrig', physicalChannel='BBNAPS1-2m1')

    ChannelDict['BBNAPS1-12'] = PhysicalQuadratureChannel(name='BBNAPS1-12', AWG='BBNAPS1', IChannel='ch1', QChannel='ch2', delay=0e-9, ampFactor=1, phaseSkew=0)
    ChannelDict['BBNAPS1-34'] = PhysicalQuadratureChannel(name='BBNAPS1-34', AWG='BBNAPS1', IChannel='ch3', QChannel='ch4', delay=0e-9, ampFactor=1, phaseSkew=0)
    ChannelDict['BBNAPS1-2m1'] = PhysicalMarkerChannel(name='BBNAPS1-2m1', AWG='BBNAPS1')

    ChannelDict['QPC1-1691'] = Generator(name='QPC1-1691', gateChannel='TekAWG1-ch1m1', gateDelay=-50.0e-9, gateBuffer=20e-9, gateMinWidth=100e-9)   
    ChannelDict['Agilent1'] = Generator(name='Agilent1', gateChannel='TekAWG1-ch1m1', gateDelay=-50.0e-9, gateBuffer=20e-9, gateMinWidth=100e-9)   

    ChannelDict['TekAWG1'] = AWG(name='TekAWG1', model='Tek5000')
    ChannelDict['BBNAPS1'] = AWG(name='BBNAPS1', model='BBNAPS')

    save_channel_info()
    update_channel_info()

    #Look to see if iPython's event loop is running
    app = QtCore.QCoreApplication.instance()
    if app is None:
        app = QtGui.QApplication(sys.argv)

    channelWindow = ChannelInfoView()
    channelWindow.show()

    try: 
        from IPython.lib.guisupport import start_event_loop_qt4
        start_event_loop_qt4(app)
    except ImportError:
        sys.exit(app.exec_())


    
    