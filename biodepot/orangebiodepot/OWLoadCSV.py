from Orange.widgets import widget, gui
import Orange.data
import sys, os
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import pyqtSlot, Qt, QThread, pyqtSignal
from PyQt5.QtGui import QStandardItem, QTextDocument, QFont

class OWLoadCSV(widget.OWWidget):
    name = "Load CSV"
    description = "Load a CSV file to Data Table"
    category = "General"
    icon = "icons/csvfile.svg"

    priority = 3

    inputs = [("CSV file", str, "setCSVfile")]
    outputs = [("DataTable", Orange.data.Table)]

    want_main_area = False
    want_control_area = True

    def __init__(self):
        super().__init__()

        self.csv_filename = ''

        self.edtCSVFile = gui.lineEdit(None, self, "csv_filename")
        self.btnSelectFile = gui.button(None, self, "☰", callback=self._OnChooseCSV, autoDefault=False)

        self.buttonsArea.layout().addWidget(self.edtCSVFile)
        self.buttonsArea.layout().addWidget(self.btnSelectFile)
        self.buttonsArea.layout().setSpacing(4)
        self.buttonsArea.layout().addSpacing(4)
        self.buttonsArea.setMinimumWidth(400)

    def _OnChooseCSV(self):
        defaultDir = '~'
        if os.path.exists('/data'):
            defaultDir = '/data'

        file = QtWidgets.QFileDialog.getOpenFileName(self, caption="Open CSV File", directory=defaultDir,
                                                         filter="CSV|TSV file (*.csv *.tsv)")
        if type(file[0]) is str:
            file = file[0].strip()
            self.setCSVfile(file)

    def setCSVfile(self, file):
        self.csv_filename = file

        tsvReader = Orange.data.io.FileFormat.get_reader(self.csv_filename)
        data = None
        try:
            data = tsvReader.read()
        except Exception as ex:
            print(ex)

        if data:
            self.send("DataTable", data)