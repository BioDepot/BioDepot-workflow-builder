import sys, os
from Orange.widgets import widget, gui, settings
from PyQt5 import QtWidgets, QtGui

class OWDirectory(widget.OWWidget):
    name = "Directory"
    description = "Set an input as a Directory"
    category = "Data"
    icon = "icons/directory.svg"
    priority = 2

    inputs = [("Dir", str, "set_dir")]
    outputs = [("Dir", str)]

    want_main_area = False
    want_control_area = True

    # Jimmy March-28-2017, persisting path when we reload workflow, the path will still there.
    #    use schema_only to ensure the stored path only loaded from .ows files
    dirname= settings.Setting('', schema_only=True)
    inputConnections = settings.Setting({},schema_only=True)
    
    def __init__(self):
        super().__init__()
        icon=QtGui.QIcon('/biodepot/orangebiodepot/icons/bluefile.png')
        disabledFlag=True
        if not self.inputConnections:
            disabledFlag=False
        self.dir_edit = gui.lineEdit(None, self, "dirname",disabled=disabledFlag)
        self.btn_send=gui.button(None, self, " Enter ", callback=self.sendDirEdit, autoDefault=False,disabled=disabledFlag)
        self.btn_dir = gui.button(None, self, "", callback=self.get_dir, autoDefault=False,disabled=disabledFlag)
        self.btn_dir.setStyleSheet("border: 1px solid #1a8ac6; border-radius: 2px;")
        self.btn_dir.setIcon(icon)
        self.btn_send.setStyleSheet("border: 1px solid #1a8ac6; border-radius: 2px; heighr 25")
        self.dir_edit.setClearButtonEnabled(True)
        self.dir_edit.setPlaceholderText("Enter directory")
        self.buttonsArea.layout().addWidget(self.dir_edit)
        self.buttonsArea.layout().addWidget(self.btn_dir)
        self.buttonsArea.layout().addWidget(self.btn_send)
        self.buttonsArea.layout().setSpacing(4)
        self.buttonsArea.layout().addSpacing(4)
        self.buttonsArea.setMinimumWidth(400)

        if self.dirname is not None:
            self.send("Dir", self.dirname)

    """
    Called when button pushed
    """
    def get_dir(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.dirname = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate directory", directory=defaultDir)
    
    def sendDirEdit(self):
        self.send("Dir",self.dirname)
        
    def findInputValue(self):
        #checks inputConnections and returns first one for now
        if not self.inputConnections:
            return None
        for sourceId in self.inputConnections:
            return self.inputConnections[sourceId]

    def set_dir (self, path, sourceId=None):
        if self.inputConnections and path is None:
            self.inputConnections.pop(sourceId,None) 
        elif path:
            self.inputConnections[sourceId]=path
        else:
            return
        inputPath=self.findInputValue()
        self.dir_edit.setEnabled(inputPath is None)
        self.btn_dir.setEnabled(inputPath is None)
        self.btn_send.setEnabled(inputPath is None)
    
        if inputPath:
            self.dirname=inputPath
            self.send('Dir',inputPath)
