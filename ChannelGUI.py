from Channels import *
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

