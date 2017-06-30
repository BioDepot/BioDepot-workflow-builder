import sys, os

from Orange.widgets import widget, gui, settings
from PyQt5 import QtWidgets, QtGui

class OWBamFile(widget.OWWidget):
    name = "File"
    description = "Set an Data file"
    category = "Data"
    icon = "icons/Bamfile.svg"
    priority = 2

    inputs = [("File", str, "set_file")]
    outputs = [("File", str)]

    want_main_area = False
    want_control_area = True

    file_name = settings.Setting('', schema_only=True)

    def __init__(self):
        super().__init__()
        self.file_edit = gui.lineEdit(None, self, "file_name")
        self.btn_file = gui.button(None, self, "☰", callback=self.get_file, autoDefault=False)

        self.buttonsArea.layout().addWidget(self.btn_file)
        self.buttonsArea.layout().addWidget(self.file_edit)
        self.buttonsArea.layout().setSpacing(4)
        self.buttonsArea.layout().addSpacing(4)
        self.buttonsArea.setMinimumWidth(400)

        if self.file_name is not "":
            self.send("File", self.file_name)

    """
    Called when button pushed
    """
    def get_file(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        file = QtWidgets.QFileDialog.getOpenFileName(self, caption="Open Data File", directory=defaultDir, filter="Any file (*.*)", options=QtGui.QFileDialog.DontUseNativeDialog)
        self.set_file(file)

    """
    Called when input set or button pushed
    """
    def set_file(self, path):
        # Make sure path isn't empty or whitespace
        if type(path[0]) is str:
            self.file_name = path[0].strip()
            self.file_edit.setText(self.file_name)
            self.send("File", self.file_name)


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