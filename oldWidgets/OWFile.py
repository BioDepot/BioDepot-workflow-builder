import sys, os
from orangebiodepot.util.BwBase import OWBwBWidget, ConnectionDict
from Orange.widgets import widget, gui, settings
from PyQt5.QtCore import Qt
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

    filesList = settings.Setting([], schema_only=True)
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
        browseIcon=QtGui.QIcon('/biodepot/Miscellaneous/icons/bluefile.png')
        addIcon=QtGui.QIcon('/biodepot/Miscellaneous/icons/add.png')
        removeIcon=QtGui.QIcon('/biodepot/Miscellaneous/icons/remove.png')
        submitIcon=QtGui.QIcon('/biodepot/Miscellaneous/icons/submit.png')
        reloadIcon=QtGui.QIcon('/biodepot/Miscellaneous/icons/reload.png')
        disabledFlag=True
    
        if not self.inputConnections:
            disabledFlag=False
        self.filename=""
        self.file_edit = gui.lineEdit(None, self, "filename",disabled=disabledFlag)
        self.boxEdit=QtGui.QListWidget(self)
        self.boxEdit.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
        self.boxEdit.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        for f in self.filesList:
            self.boxEdit.addItem(f)
        self.btn_send=gui.button(None, self, "", callback=self.sendFileEdit, autoDefault=False,disabled=disabledFlag)
        self.btn_add=gui.button(None, self, "", callback=self.addFileEdit, autoDefault=False,disabled=disabledFlag)
        self.btn_remove=gui.button(None, self, "", callback=self.deleteFiles, autoDefault=False,disabled=disabledFlag)
        self.btn_reload=gui.button(None, self, "", callback=self.reloadAll, autoDefault=False,disabled=disabledFlag)
        self.btn_file = gui.button(None, self, "", callback=self.get_file, autoDefault=False,disabled=disabledFlag)
        if self.file_edit.isEnabled():
            self.btn_reload.setEnabled(False)
        if not self.file_edit.text():
            self.btn_add.setEnabled(False)
        if not self.boxEdit.count():
            self.btn_send.setEnabled(False)
        if not self.boxEdit.selectedItems():
            self.btn_remove.setEnabled(False)
        self.btn_file.setIcon(browseIcon)
        self.btn_add.setIcon(addIcon)
        self.btn_remove.setIcon(removeIcon)
        self.btn_reload.setIcon(reloadIcon)
        self.btn_send.setIcon(submitIcon)
        self.btn_add.setStyleSheet('background: None; border: None ; border-radius: 0;')
        self.btn_file.setStyleSheet('background: None; border: None ; border-radius: 0;')
        self.btn_remove.setStyleSheet('background: None; border: None ; border-radius: 0;')
        self.btn_send.setStyleSheet('background: None; border: None ; border-radius: 0;')
        self.btn_reload.setStyleSheet('background: None; border: None ; border-radius: 0;')
        self.file_edit.setClearButtonEnabled(True)
        self.file_edit.setPlaceholderText("Enter file(s)")
        self.file_edit.setStyleSheet(":disabled { color: #282828}")
        self.boxEdit.setStyleSheet(":disabled { color: #282828}")
        self.file_edit.textChanged.connect(lambda: self.btn_add.setEnabled(bool(self.file_edit.text())))
        self.boxEdit.itemSelectionChanged.connect(lambda: self.btn_remove.setEnabled(bool(self.boxEdit.selectedItems())))
        self.bigBox=gui.widgetBox(self.buttonsArea)
        boxLayout=QtGui.QVBoxLayout()     
        layout=QtGui.QGridLayout()
        layout.addWidget(self.file_edit,1,0)
        layout.addWidget(self.btn_file,1,1)
        layout.addWidget(self.btn_add,1,2)
        layout.addWidget(self.btn_remove,1,3)
        layout.addWidget(self.btn_send,1,4)
        layout.addWidget(self.btn_reload,1,5)
        self.buttonsArea.setMinimumWidth(400)
        self.bigBox.layout().addLayout(boxLayout)
        self.scroll_area = QtWidgets.QScrollArea(verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn)
        self.scroll_area.setWidget(self.bigBox)
        self.scroll_area.setWidgetResizable(True)
        self.buttonsArea.layout().addWidget(self.scroll_area)
        boxLayout.addWidget(self.boxEdit)
        boxLayout.addLayout(layout)
        self.buttonsArea.layout().addLayout(boxLayout)
        if self.btn_reload.isEnabled:
            self.sendFileEdit()
    """
    Called when button pushed
    """

    def deleteFiles(self):
        if self.boxEdit.selectedItems():
            for item in self.boxEdit.selectedItems():
                self.boxEdit.takeItem(self.boxEdit.row(item))
        if self.boxEdit.count():
            self.btn_send.setEnabled(True)
        else:
            self.btn_send.setEnabled(False)
            self.btn_remove.setEnabled(False)
            
    def addFileEdit(self):
        if self.filename:
            self.boxEdit.addItem(self.filename)
            self.file_edit.clear()
            self.btn_send.setEnabled(True)
            self.btn_add.setEnabled(False)
        
    def get_file(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        self.boxEdit.addItems(QtWidgets.QFileDialog.getOpenFileNames(self, caption="Open Data File", directory=defaultDir, filter="Any file (*.*)")[0])
        if self.boxEdit.count():
            self.btn_send.setEnabled(True)
            
    def sendFileEdit(self):
        if self.boxEdit.count():
            self.filesList=[]
            for index in range(self.boxEdit.count()):
                self.filesList.append(self.boxEdit.item(index).text())
            self.send("File","\n".join( self.filesList))
            self.btn_reload.setEnabled(True)
            self.btn_send.setEnabled(False)
            self.file_edit.setEnabled(False)
            self.btn_add.setEnabled(False)
            self.btn_file.setEnabled(False)
            self.boxEdit.setEnabled(False)            
            
    def reloadAll(self):
        self.btn_reload.setEnabled(False)
        self.btn_send.setEnabled(True)
        self.file_edit.setEnabled(True)
        self.file_edit.clear()
        self.btn_add.setEnabled(False)
        self.btn_file.setEnabled(True)
        self.boxEdit.setEnabled(True)
        self.send("File",None)
        
    def set_file (self, path):
        if path is None:
            self.inputConnections is None
            self.btn_reload.setEnabled(False)
            self.btn_send.setEnabled(True)
            self.file_edit.setEnabled(True)
            self.file_edit.clear()
            self.btn_add.setEnabled(False)
            self.btn_file.setEnabled(True)
        elif path:
            self.inputConnections=path
            self.boxEdit.addItems(str.splitlines(path))
            self.filesList=str.splitlines(path)
            self.btn_reload.setEnabled(False)
            self.btn_send.setEnabled(False)
            self.file_edit.setEnabled(False)
            self.file_edit.clear()
            self.btn_add.setEnabled(False)
            self.btn_file.setEnabled(False)
            self.boxEdit.setEnabled(False)  
            self.send("File",path)


