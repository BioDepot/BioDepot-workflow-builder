import sys
import numpy

import Orange.data
from Orange.widgets import widget, gui
from PyQt5 import QtWidgets

class OWBamFile(widget.OWWidget):
    name = "Bam File"
    description = "Set an bam file"
    category = "Data"
    icon = "icons/bamfile.svg"
    priority = 2

    inputs = [("File", str, "set_file")]
    outputs = [("File", str)]

    want_main_area = False
    want_control_area = True

    def __init__(self):
        super().__init__()
        self.file_edit = QtWidgets.QLineEdit()
        self.btn_file = gui.button(None, self, "☰", callback=self.get_file, autoDefault=False)

        self.buttonsArea.layout().addWidget(self.btn_file)
        self.buttonsArea.layout().addWidget(self.file_edit)
        self.buttonsArea.layout().addSpacing(8)
        self.buttonsArea.setMinimumWidth(400)

    """
    Called when button pushed
    """
    def get_file(self):
        file = QtWidgets.QFileDialog.getOpenFileName(self, "Open BAM File", ".", "BAM file (*.bam)")
        self.set_file(file)

    """
    Called when input set or button pushed
    """
    def set_file(self, path):
        # Make sure path isn't empty or whitespace
        if type(path[0]) is str:
            samfile = path[0].strip()
            self.file_edit.setText(samfile)
            self.send("File", samfile)


def main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv))

    ow = OWBamFile()
    ow.show()
    ow.raise_()

    ow.handleNewSignals()
    app.exec_()
    ow.handleNewSignals()
    return 0

if __name__=="__main__":
    sys.exit(main())