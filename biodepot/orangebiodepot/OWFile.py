import sys, os
from orangebiodepot.util.BwBase import OWBwBWidget, ConnectionDict
from Orange.widgets import widget, gui, settings
from PyQt5 import QtWidgets, QtGui

class OWBamFile(widget.OWWidget):
    name = "File"
    description = "Set one or more data files"
    category = "Data"
    icon = "icons/Bamfile.svg"
    priority = 2

    inputs = [("File", str, "set_file")]
    outputs = [("File", str)]

    want_main_area = False
    want_control_area = True

    filename = settings.Setting('', schema_only=True)
    inputConnections = settings.Setting({},schema_only=True)

    def __init__(self):
        super().__init__()
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid #1a8ac6; border-radius: 2px;}
        QPushButton:pressed { background-color: #158805; border-style: inset;}
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        QPushButton:hover {background-color: #1588f5; }
        '''  
        self.setStyleSheet(css)
        icon=QtGui.QIcon('/biodepot/orangebiodepot/icons/bluefile.png')
        disabledFlag=True
        
        if not self.inputConnections:
            disabledFlag=False
        self.file_edit = gui.lineEdit(None, self, "filename",disabled=disabledFlag)
        self.btn_send=gui.button(None, self, " Enter ", callback=self.sendFileEdit, autoDefault=False,disabled=disabledFlag)
        self.btn_file = gui.button(None, self, "", callback=self.get_file, autoDefault=False,disabled=disabledFlag)
        self.btn_file.setStyleSheet("border: 1px solid #1a8ac6; border-radius: 2px;")
        self.btn_file.setIcon(icon)
        self.btn_send.setStyleSheet("border: 1px solid #1a8ac6; border-radius: 2px; heighr 25")
        self.file_edit.setClearButtonEnabled(True)
        self.file_edit.setPlaceholderText("Enter file(s)")
        self.buttonsArea.layout().addWidget(self.file_edit)
        self.buttonsArea.layout().addWidget(self.btn_file)
        self.buttonsArea.layout().addWidget(self.btn_send)
        self.buttonsArea.layout().setSpacing(4)
        self.buttonsArea.layout().addSpacing(4)
        self.buttonsArea.setMinimumWidth(400)

        if self.filename is not None:
            self.send("File", self.filename)

    """
    Called when button pushed
    """
    def get_file(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.filename = ' '.join(QtWidgets.QFileDialog.getOpenFileNames(self, caption="Open Data File", directory=defaultDir, filter="Any file (*.*)")[0])
    
    def sendFileEdit(self):
        self.send("File",self.filename)
        
    def findInputValue(self):
        #checks inputConnections and returns first one for now
        if not self.inputConnections:
            return None
        for sourceId in self.inputConnections:
            return self.inputConnections[sourceId]

    def set_file (self, path, sourceId=None):
        if self.inputConnections and path is None:
            self.inputConnections.pop(sourceId,None) 
        elif path:
            self.inputConnections[sourceId]=path
        else:
            return
        inputPath=self.findInputValue()
        self.file_edit.setEnabled(inputPath is None)
        self.btn_file.setEnabled(inputPath is None)
        self.btn_send.setEnabled(inputPath is None)
    
        if inputPath:
            self.filename=inputPath
            self.send('File',inputPath)

