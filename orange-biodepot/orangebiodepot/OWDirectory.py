import sys
import numpy

import Orange.data
from Orange.widgets import widget, gui
from PyQt4 import QtGui

class OWDtoxsAlignment(widget.OWWidget):
    name = "Directory"
    description = "Set an input as a Directory"
    category = "Data"
    icon = "icons/directory.svg"
    priority = 2

    inputs = [("Dir", str, "set_dir")]
    outputs = [("Dir", str)]

    want_main_area = False
    want_control_area = True

    def __init__(self):
        super().__init__()
        self.dir_edit = QtGui.QLineEdit()
        self.btn_dir = gui.button(None, self, "☰", callback=self.get_dir, autoDefault=False)

        self.buttonsArea.layout().addWidget(self.btn_dir)
        self.buttonsArea.layout().addWidget(self.dir_edit)
        self.buttonsArea.layout().addSpacing(8)
        self.buttonsArea.setMinimumWidth(400)

    """
    Called when button pushed
    """
    def get_dir(self):
        dir  = QtGui.QFileDialog.getExistingDirectory(self)
        self.set_dir(dir)

    """
    Called when input set or button pushed
    """
    def set_dir(self, path):
        # Make sure path isn't empty or whitespace
        if type(path) is str and path.strip():
            path = path.strip()
            self.dir_edit.setText(path)
            self.send("Dir", path)