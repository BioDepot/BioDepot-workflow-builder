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
    directory_path = settings.Setting('', schema_only=True)

    def __init__(self):
        super().__init__()

        self.dir_edit = gui.lineEdit(None, self, "directory_path") #QtWidgets.QLineEdit()
        self.btn_dir = gui.button(None, self, "☰", callback=self.get_dir, autoDefault=False)

        self.buttonsArea.layout().addWidget(self.btn_dir)
        self.buttonsArea.layout().addWidget(self.dir_edit)
        self.buttonsArea.layout().setSpacing(4)
        self.buttonsArea.layout().addSpacing(4)
        self.buttonsArea.setMinimumWidth(400)

        # Jimmy March-28-2017, if we loaded settings from workflow, trigger the output channel
        if self.directory_path is not "":
            self.send("Dir", self.directory_path)

    """
    Called when button pushed
    """
    def get_dir(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        dir = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate Directory", directory=defaultDir)
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


def main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv))

    ow = OWDirectory()
    ow.show()
    ow.raise_()

    ow.handleNewSignals()
    app.exec_()
    ow.handleNewSignals()
    return 0

if __name__=="__main__":
    sys.exit(main())