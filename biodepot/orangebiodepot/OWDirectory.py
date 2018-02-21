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
        self.dir_edit = gui.lineEdit(None, self, "dirname",disabled=disabledFlag)
        self.btn_send=gui.button(None, self, " Enter ", callback=self.sendDirEdit, autoDefault=False,disabled=disabledFlag)
        self.btn_dir = gui.button(None, self, "", callback=self.get_dir, autoDefault=False,disabled=disabledFlag)
        self.btn_dir.setStyleSheet("border: 1px solid #1a8ac6; border-radius: 2px;")
        self.btn_dir.setIcon(icon)
        self.btn_send.setStyleSheet("border: 1px solid #1a8ac6; border-radius: 2px;")
        self.dir_edit.setClearButtonEnabled(True)
        self.dir_edit.setPlaceholderText("Enter directory")
        self.checkbox=gui.checkBox(None, self,'checked',label=None)
        self.checkbox.setStyleSheet("QCheckBox::indicator{width: 20; height: 20}")
        self.dir_edit.setEnabled(self.checkbox.isChecked())
        self.btn_dir.setEnabled(self.checkbox.isChecked())
        self.btn_send.setEnabled(self.checkbox.isChecked())
        self.dir_edit.setStyleSheet(":disabled { color: #282828}")
        self.checkbox.stateChanged.connect(self.onCheckboxChange)
        self.buttonsArea.layout().addWidget(self.checkbox)
        self.buttonsArea.layout().addWidget(self.dir_edit)
        self.buttonsArea.layout().addWidget(self.btn_dir)
        self.buttonsArea.layout().addWidget(self.btn_send)
        self.buttonsArea.layout().setSpacing(4)
        self.buttonsArea.layout().addSpacing(4)
        self.buttonsArea.setMinimumWidth(400)
        self.oldDirname=self.dirname
        if self.dirname:
            self.send("Dir", self.dirname)
            self.btn_send.setText('ReEnter')
        else:
            self.dir_edit.clear()
            self.btn_send.setText('Enter')

    """
    Called when button pushed
    """
    def onCheckboxChange(self):
        self.dir_edit.setEnabled(self.checkbox.isChecked())
        self.btn_dir.setEnabled(self.checkbox.isChecked())
        self.btn_send.setEnabled(self.checkbox.isChecked())
        if not self.checkbox.isChecked():
            self.dirname=self.oldDirname
            if self.oldDirname:
                self.dir_edit.setText(self.oldDirname)
            else:
                self.dir_edit.clear()
                self.btn_send.setText('Enter')
        
    def get_dir(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.dirname = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate directory", directory=defaultDir)
    
    def sendDirEdit(self):
        self.oldDirname=self.dirname
        if self.dirname:
            self.send("Dir",self.dirname)
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

    def set_dir (self, path):
        if path is None:
            self.inputConnections is None
        elif path:
            self.inputConnections=path
            self.dir_edit.setText(path)
            self.sendDirEdit()

