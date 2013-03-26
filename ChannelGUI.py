'''
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
from Channels import *
import PulseShapes 

import inspect
from PySide import QtGui, QtCore




class ChannelInfoView(QtGui.QMainWindow):
    def __init__(self):
        super(ChannelInfoView, self).__init__()

        #Create an item view for the logical channels
        self.logicalChannelListModel = QtGui.QStringListModel([tmpKey for tmpKey, tmpChan in ChannelDict.items() if isinstance(tmpChan, LogicalChannel)])
        self.logicalChannelListModel.sort(0)
        self.logicalChannelListView = QtGui.QListView()
        self.logicalChannelListView.setModel(self.logicalChannelListModel)
        self.logicalChannelListView.clicked.connect(lambda(idx): self.update_channelView(idx, self.logicalChannelListModel))
        
        #Create an item view for the physical channels
        self.physicalChannelListModel = QtGui.QStringListModel([tmpKey for tmpKey, tmpChan in ChannelDict.items() if isinstance(tmpChan, PhysicalChannel)])
        self.physicalChannelListModel.sort(0)
        self.physicalChannelListView = QtGui.QListView()
        self.physicalChannelListView.setModel(self.physicalChannelListModel)
        self.physicalChannelListView.clicked.connect(lambda(idx): self.update_channelView(idx, self.physicalChannelListModel))

        #Create an item view for the generators channels
        self.generatorListModel = QtGui.QStringListModel([tmpKey for tmpKey, tmpChan in ChannelDict.items() if isinstance(tmpChan, Generator)])
        self.generatorListModel.sort(0)
        self.generatorListView = QtGui.QListView()
        self.generatorListView.setModel(self.generatorListModel)
        self.generatorListView.clicked.connect(lambda(idx): self.update_channelView(idx, self.generatorListModel))

        #Create an item view for the AWG channels
        self.awgListModel = QtGui.QStringListModel([tmpKey for tmpKey, tmpChan in ChannelDict.items() if isinstance(tmpChan, AWG)])
        self.awgListModel.sort(0)
        self.awgListView = QtGui.QListView()
        self.awgListView.setModel(self.awgListModel)
        self.awgListView.clicked.connect(lambda(idx): self.update_channelView(idx, self.awgListModel))

        tmpWidget = QtGui.QWidget()
        vBox = QtGui.QVBoxLayout(tmpWidget)
        
        self.tabWidget = QtGui.QTabWidget()
        self.tabWidget.addTab(self.logicalChannelListView, 'Logical')        
        self.tabWidget.addTab(self.physicalChannelListView, 'Physical')        
        self.tabWidget.addTab(self.generatorListView, 'Generator')
        self.tabWidget.addTab(self.awgListView, 'AWG')
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

    #Make a list of all the pulse names
    pulseNames = [x[0] for x in inspect.getmembers(PulseShapes) if inspect.isfunction(x[1])]
    pulseNameModel = QtGui.QStringListModel(pulseNames)

    def __init__(self, channel, parent = None):
        super(ChannelView, self).__init__()
        self.channel = channel
        
        #Some hidden fields
        skipFields = ['name', 'pulseParams']

        #Create the layout as a vbox of hboxes
        form1 = QtGui.QFormLayout()
        self.GUIhandles = {}
        #Do the channelType on top for information purposes
        form1.addRow('channelType', QtGui.QLabel(self.channel.__class__.__name__))

        #Helper function to update         
        for key,value in sorted(channel.__dict__.items(), key=itemgetter(0)):
            if key not in skipFields:
                #For the physical channel field we'll pop up a combo box with all the other channels that match this one
                if key == '_physicalChannel':
                    chanType = type(ChannelDict[value])
                    tmpModel = QtGui.QStringListModel()
                    #Make sure that if we delete the channel then it gets removed from here too
                    parent.physicalChannelListModel.rowsRemoved.connect(lambda a,b,c : tmpModel.removeRow(b))
                    tmpModel.setStringList( [tmpChan for tmpChan in parent.physicalChannelListModel.stringList() if type(ChannelDict[tmpChan]) == chanType] )
                    tmpWidget = QtGui.QComboBox()
                    tmpWidget.setModel(tmpModel)
                elif isinstance(value, basestring):
                    tmpWidget = QtGui.QLineEdit(value)
                else:
                    tmpWidget = QtGui.QLineEdit(str(value))
                    tmpWidget.setValidator(QtGui.QDoubleValidator())
                form1.addRow(key.lstrip('_'), tmpWidget)
                self.GUIhandles[key] = tmpWidget

        vBox = QtGui.QVBoxLayout()
        vBox.addLayout(form1)

        #Now do the pulse params
        form2 = QtGui.QFormLayout()
        if hasattr(channel, 'pulseParams'):
            vBox2 = QtGui.QVBoxLayout()
            for key,value in channel.pulseParams.items():
                #List all available shape functions
                if key == 'shapeFun':
                    tmpWidget = QtGui.QComboBox()
                    tmpWidget.setModel(self.pulseNameModel)
                elif isinstance(value, basestring):
                    tmpWidget = QtGui.QLineEdit(value)
                else:
                    tmpWidget = QtGui.QLineEdit(str(value))
                    tmpWidget.setValidator(QtGui.QDoubleValidator())
                form1.addRow(key, tmpWidget)
                self.GUIhandles[key] = tmpWidget
            
            #Some add/delete buttons
            hBox = QtGui.QHBoxLayout()
            addChanButton = QtGui.QPushButton('Add')
            # addChanButton.clicked.connect(self.add_channel)
            hBox.addWidget(addChanButton)
            deleteChanButton = QtGui.QPushButton('Delete')
            # deleteChanButton.clicked.connect(self.delete_channel)
            hBox.addWidget(deleteChanButton)
            hBox.addStretch(1)
            vBox2.addLayout(hBox)                


            vBox.addLayout(vBox2)

        self.setLayout(vBox)
            
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

