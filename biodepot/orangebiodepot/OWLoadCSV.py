from Orange.widgets import widget, gui, settings
import Orange.data
import sys, os
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QComboBox

class OWLoadCSV(widget.OWWidget):
    name = "Load CSV"
    description = "Load a CSV file to Data Table"
    category = "General"
    icon = "icons/csvfile.svg"

    priority = 3

    inputs = [("CSV file", str, "setCSVfile"), ("Directory", str, "setDirectory")]
    outputs = [("DataTable", Orange.data.Table)]

    want_main_area = False
    want_control_area = True

    FROM_FOLDER, SPECIFIC_FILE = range(2)
    source = settings.Setting('', schema_only=True)
    choosePath = settings.Setting('', schema_only=True)
    chooseFile = settings.Setting('', schema_only=True)

    def __init__(self):
        super().__init__()

        self.csv_filename = ''
        self.folder_path = ''

        self.edtCSVFile = gui.lineEdit(None, self, "csv_filename")
        self.btnSelectFile = gui.button(None, self, "☰", callback=self._OnChooseCSV, autoDefault=False)
        self.edtFolder = gui.lineEdit(None, self, "choosePath")
        self.btnSelFolder =  gui.button(None, self, "☰", callback=self._OnChooseFolder, autoDefault=False)
        self.cboFileList = QComboBox(self, sizeAdjustPolicy=QComboBox.AdjustToContents)
        self.cboFileList.activated[int].connect(self.load_data)

        #info_box = gui.widgetBox(self.controlArea, "Info")

        vlayoutBase = QtWidgets.QVBoxLayout()
        vlayoutBase.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        vlayoutBase.setContentsMargins(0, 0, 0, 0)
        vlayoutBase.setSpacing(0)


        radioGroup = gui.radioButtons(None, self, "source", box = True, addSpace=False,
                                      callback=self.OnFileTypeChanged, addToLayout=False)

        # Folder section
        hbox_folderbutton = gui.hBox(None, addToLayout=False, margin=0)
        folderButton = gui.appendRadioButton(radioGroup, "Folder:", addToLayout=False)
        folderButton.setFixedSize(QtCore.QSize(85, 25))
        hbox_folderbutton.layout().addWidget(folderButton)
        hbox_folder = gui.hBox(None, addToLayout=False, margin=0)
        hbox_folder.layout().setSpacing(0)

        hbox_folderEditor = gui.hBox(None, addToLayout=False, margin=0)
        hbox_folderEditor.layout().addWidget(self.edtFolder)
        hbox_folderEditor.layout().addWidget(self.btnSelFolder)

        vbox_folderRight = gui.vBox(None, addToLayout=False, margin=0)
        vbox_folderRight.layout().setSpacing(0)
        vbox_folderRight.layout().addWidget(hbox_folderEditor)
        vbox_folderRight.layout().addWidget(self.cboFileList)

        hbox_folder.layout().addWidget(hbox_folderbutton)
        hbox_folder.layout().addWidget(vbox_folderRight)

        # File section
        hbox_filebutton = gui.hBox(None, addToLayout=False, margin=0)
        fileButton = gui.appendRadioButton(radioGroup, "File:", addToLayout=False)
        fileButton.setFixedSize(QtCore.QSize(85, 25))
        hbox_filebutton.layout().addWidget(fileButton)

        hbox_file = gui.hBox(None, addToLayout=False, margin=0)
        hbox_file.layout().setSpacing(0)

        fileArea = gui.hBox(None, addToLayout=False, margin = 0)
        fileArea.layout().addWidget(self.edtCSVFile)
        fileArea.layout().addWidget(self.btnSelectFile)

        hbox_file.layout().addWidget(hbox_filebutton)
        hbox_file.layout().addWidget(fileArea)

        vlayoutBase.addWidget(hbox_folder)
        vlayoutBase.addWidget(hbox_file)
        #vlayoutBase.addWidget(info_box)

        self.controlArea.layout().addLayout(vlayoutBase)

        # set window size
        self.controlArea.setMinimumSize(500,110)

        if not self.source:
            self.source = self.FROM_FOLDER
            self.OnFileTypeChanged()

    def _OnChooseCSV(self):
        defaultDir = '~'
        if os.path.exists('/data'):
            defaultDir = '/data'

        file = QtWidgets.QFileDialog.getOpenFileName(self, caption="Open CSV File", directory=defaultDir,
                                                         filter="CSV|TSV file (*.csv *.tsv)")
        if not file:
            return

        if type(file[0]) is str:
            file = file[0].strip()
            self.setCSVfile(file)

    def _OnChooseFolder(self):
        defaultDir = '/root'
        if os.path.exists('/data'):
            defaultDir = '/data'
        dir = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate Directory", directory=defaultDir)
        if not dir:
            return

        self.setDirectory(dir)

    def _scanDirectory(self):
        self.filelist = []
        for root, dirs, files in os.walk(self.choosePath):
            for filename in files:
                if filename.endswith(('.csv', '.tsv')):
                    self.filelist.append(filename)

        self.filelist = [os.path.join(self.choosePath, x) for x in self.filelist]

        self.cboFileList.clear()
        if not self.filelist:
            self.cboFileList.addItem("(none)")
            self.cboFileList.model().item(0).setEnabled(False)
        else:
            for i, file in enumerate(self.filelist):
                self.cboFileList.addItem(os.path.basename(file))
                self.cboFileList.model().item(i).setToolTip(file)
                self.cboFileList.setItemData(i, file)

    def setDirectory(self, path):
        if path is None:
            self.source = self.SPECIFIC_FILE
            return

        if not os.path.exists(path):
            return

        self.source = self.FROM_FOLDER
        self.choosePath = path
        self.OnFileTypeChanged()

    def setCSVfile(self, file):
        if not file:
            return

        self.source = self.SPECIFIC_FILE
        self.csv_filename = file
        self.OnFileTypeChanged()


    def OnFileTypeChanged(self):
        if self.source is self.FROM_FOLDER:
            self.edtFolder.setEnabled(True)
            self.cboFileList.setEnabled(True)
            self.edtCSVFile.setEnabled(False)
            self.btnSelFolder.setEnabled(True)
            self.btnSelectFile.setEnabled(False)
            self._scanDirectory()
        else:
            self.edtFolder.setEnabled(False)
            self.cboFileList.setEnabled(False)
            self.edtCSVFile.setEnabled(True)
            self.btnSelFolder.setEnabled(False)
            self.btnSelectFile.setEnabled(True)

        self.load_data()

    def load_data(self):
        filename = ''
        if self.source is self.FROM_FOLDER:
            filename = self.cboFileList.itemData(self.cboFileList.currentIndex())
            print(filename)
        else:
            if self.csv_filename:
                filename = self.csv_filename

        if filename and os.path.exists(filename):
            data = None
            try:
                tsvReader = Orange.data.io.FileFormat.get_reader(filename)
                data = tsvReader.read()
            except Exception as ex:
                print(ex)

            if data:
                self.send("DataTable", data)



def main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv))

    ow = OWLoadCSV()
    #ow.setDirectory('/Users/Jimmy/Downloads/STAR/')
    #ow.setHostMountedDir("/Users/Jimmy/Downloads/STAR/fasta_input", None, None)
    ow.show()
    ow.raise_()

    ow.handleNewSignals()
    app.exec_()
    return 0

if __name__=="__main__":
    sys.exit(main())