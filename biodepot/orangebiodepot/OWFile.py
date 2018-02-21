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
    inputConnections = settings.Setting(None,schema_only=True)
    checked=settings.Setting(False,schema_only=True)
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
        self.btn_send.setStyleSheet("border: 1px solid #1a8ac6; border-radius: 2px;")
        self.file_edit.setClearButtonEnabled(True)
        self.file_edit.setPlaceholderText("Enter file(s)")
        self.checkbox=gui.checkBox(None, self,'checked',label=None)
        self.checkbox.setStyleSheet("QCheckBox::indicator{width: 20; height: 20}")
        self.file_edit.setEnabled(self.checkbox.isChecked())
        self.btn_file.setEnabled(self.checkbox.isChecked())
        self.btn_send.setEnabled(self.checkbox.isChecked())
        self.file_edit.setStyleSheet(":disabled { color: #282828}")
        self.checkbox.stateChanged.connect(self.onCheckboxChange)
        self.buttonsArea.layout().addWidget(self.checkbox)
        self.buttonsArea.layout().addWidget(self.file_edit)
        self.buttonsArea.layout().addWidget(self.btn_file)
        self.buttonsArea.layout().addWidget(self.btn_send)
        self.buttonsArea.layout().setSpacing(4)
        self.buttonsArea.layout().addSpacing(4)
        self.buttonsArea.setMinimumWidth(400)
        self.oldFilename=self.filename
        if self.filename:
            self.send("File", self.filename)
            self.btn_send.setText('ReEnter')
        else:
            self.file_edit.clear()
            self.btn_send.setText('Enter')

    """
    Called when button pushed
    """
    def onCheckboxChange(self):
        self.file_edit.setEnabled(self.checkbox.isChecked())
        self.btn_file.setEnabled(self.checkbox.isChecked())
        self.btn_send.setEnabled(self.checkbox.isChecked())
        if not self.checkbox.isChecked():
            self.filename=self.oldFilename
            if self.oldFilename:
                self.file_edit.setText(self.oldFilename)
            else:
                self.file_edit.clear()
                self.btn_send.setText('Enter')
        
    def get_file(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.filename = ' '.join(QtWidgets.QFileDialog.getOpenFileNames(self, caption="Open Data File", directory=defaultDir, filter="Any file (*.*)")[0])
    
    def sendFileEdit(self):
        self.oldFilename=self.filename
        if self.filename:
            self.send("File",self.filename)
            self.checkbox.setChecked(False)
            self.btn_send.setText('ReEnter')
        else:
            self.btn_send.setText('Enter')

        
    def findInputValue(self):
        #checks inputConnections and returns first one for now
        if not self.inputConnections:
            self.checkbox.setChecked(False)
            return None
        for sourceId in self.inputConnections:
            self.checkbox.setEnabled(False)
            self.checkbox.setChecked(False)
            return self.inputConnections[sourceId]

    def set_file (self, path):
        if path is None:
            self.inputConnections is None
        elif path:
            self.inputConnections=path
            self.file_edit.setText(path)
            self.sendFileEdit()


