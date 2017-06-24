import sys, os
from Orange.widgets import widget
from orangebiodepot.util.DockerClient import DockerClient
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import QStandardItem, QFont
from PyQt5.QtCore import QThread

class OWGenericTask(widget.OWWidget):
    name = "Customized Container"
    description = "Run docker container"
    category = "General"
    icon = "icons/container.svg"

    priority = 2

    inputs =  [("Mounted Directory", str, "setHostMountedDir", widget.Multiple)]
    outputs = [("Output", str)]

    want_main_area = False


    dockerClient = DockerClient('unix:///var/run/docker.sock', 'local')

    def __init__(self):
        super().__init__()

        self.base_dir = os.path.dirname(os.path.abspath(__file__))

        # stylesheet
        css = '''
            QPushButton {background-color: #1588c5; color: white; height: 25px; border: 1px solid #1a8ac6; border-radius: 2px;}
            QPushButton:pressed { background-color: #158805; border-style: inset;}
            QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
            QPushButton:hover {background-color: #1588f5; }   
            QListView::item{color: black;margin: 5px;}
            QFrame#frameTitle {background: #1588c5; color: #1588c5}
            QFrame#frameSubtitle {background: #1998de; }
            QFrame#frameSearchName {background: #1998de; }
            QFrame[frameShape="4"] { color: #c6d1da; }
            QLabel#lblFormTitle{ background-color: #daeaf7;padding-left: 5px;}
            QFrame#mainContent {background-color: #1588f5;border: 1px solid #c6d1da;}
            QComboBox {border: 1px solid #bdbdbd;border-radius: 3px;padding-left: 5px;color: #181818;height: 30px;}
            QComboBox::drop-down {border: 0px;width: 20px;margin-right:2px;}
            QComboBox::down-arrow {border: 0px;image: url(icons/ui_combo_dropdown_arrow.png);background-position: center center; height: 10px; width: 10px}
            QComboBox::down-arrow:on { top: 1px;left: 1px; }
            QComboBox QAbstractItemView::item {min-height: 12px;}
        '''
        self.setStyleSheet(css)
        self.vlayoutBase = QtWidgets.QVBoxLayout(self.controlArea)
        self.vlayoutBase.setContentsMargins(0, 0, 0, 0)
        self.vlayoutBase.setSpacing(0)
        self.vlayoutBase.setObjectName("vlayoutBase")

        self.frameTitle = QtWidgets.QFrame()
        self.frameTitle.setMinimumHeight(35)
        self.frameTitle.setMaximumHeight(35)
        self.frameTitle.setObjectName("frameTitle")
        self.lblIcon = QtWidgets.QLabel()
        self.lblIcon.setFixedSize(QtCore.QSize(45, 45))
        self.lblIcon.setText("")
        pixmap = QtGui.QPixmap(os.path.join(self.base_dir, "icons/generic_task.svg"))
        pixmap.scaled(self.lblIcon.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
        self.lblIcon.setPixmap(pixmap)
        self.lblIcon.setScaledContents(True)
        self.lblIcon.setObjectName("lblIcon")
        self.lblFormTitle = QtWidgets.QLabel()
        self.lblFormTitle.setObjectName("lblFormTitle")

        self.hlayout_title = QtWidgets.QHBoxLayout(self.frameTitle)
        self.hlayout_title.setContentsMargins(0, 0, 0, 0)
        self.hlayout_title.setSpacing(0)
        self.hlayout_title.addWidget(self.lblIcon)
        self.hlayout_title.addWidget(self.lblFormTitle)

        self.frameTitle.setLayout(self.hlayout_title)
        self.vlayoutBase.addWidget(self.frameTitle)

        self.frameLine = QtWidgets.QFrame()
        self.frameLine.setFrameShape(QtWidgets.QFrame.HLine)
        self.frameLine.setFrameShadow(QtWidgets.QFrame.Plain)
        self.vlayoutBase.addWidget(self.frameLine)

        self.mainContent = QtWidgets.QFrame()
        self.mainContent.setMinimumHeight(60)
        self.mainContent.setContentsMargins(3,3,3,3)
        self.vlayout_main = QtWidgets.QVBoxLayout()
        self.vlayout_main.setContentsMargins(8, 8, 8, 8)
        self.vlayout_main.setSpacing(2)

        self.lblDockerImage = QtWidgets.QLabel()
        self.lblDockerImage.setMinimumHeight(25)
        self.lblDockerImage.setMaximumHeight(25)
        self.vlayout_main.addWidget(self.lblDockerImage)

        self.cboDockerImage = QtWidgets.QComboBox(self.mainContent)
        self.cboDockerImage.setMinimumHeight(30)
        self.cboDockerImage.setMaximumHeight(30)
        self.cboDockerImage.setObjectName("cboDockerImage")
        self.delegate4Combo = QtWidgets.QStyledItemDelegate()
        self.cboDockerImage.setItemDelegate(self.delegate4Combo)
        self.vlayout_main.addWidget(self.cboDockerImage)

        self.vlayout_main.addItem(
            QtWidgets.QSpacerItem(-1, 10, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum))

        self.hlayout_mapping = QtWidgets.QHBoxLayout()
        self.hlayout_mapping.setObjectName("hlayout_mapping")
        self.hlayout_mapping.setSpacing(5)

        self.lblVMapping = QtWidgets.QLabel()
        self.lblVMapping.setMinimumHeight(25)
        self.lblVMapping.setMaximumHeight(25)
        self.hlayout_mapping.addWidget(self.lblVMapping)

        self.hlayout_mapping.addItem(
            QtWidgets.QSpacerItem(40, -1, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum))

        self.btnAddMapping = QtWidgets.QPushButton(self.mainContent)
        self.btnAddMapping.setFixedSize(QtCore.QSize(20, 20))
        self.btnAddMapping.setObjectName("btnAddMapping")
        self.btnDeleteMapping = QtWidgets.QPushButton(self.mainContent)
        self.btnDeleteMapping.setFixedSize(QtCore.QSize(20, 20))
        self.btnDeleteMapping.setObjectName("btnDeleteMapping")
        self.hlayout_mapping.addWidget(self.btnAddMapping)
        self.hlayout_mapping.addWidget(self.btnDeleteMapping)

        self.vlayout_main.addLayout(self.hlayout_mapping)

        self.mainContent.setLayout(self.vlayout_main)
        self.vlayoutBase.addWidget(self.mainContent)

        self.lstMapping = QtWidgets.QTableView(self.mainContent)
        self.lstMapping.setObjectName("lstMapping")
        self.lstMapping.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.lstMapping.verticalHeader().setVisible(False)
        self.vlayout_main.addWidget(self.lstMapping)

        self.lblCmd = QtWidgets.QLabel()
        self.lblCmd.setMinimumHeight(25)
        self.lblCmd.setMaximumHeight(25)
        self.vlayout_main.addWidget(self.lblCmd)

        self.txtCommand = QtWidgets.QPlainTextEdit(self.mainContent)
        self.txtCommand.setObjectName("txtCommand")
        self.vlayout_main.addWidget(self.txtCommand)

        self.vlayout_main.addItem(QtWidgets.QSpacerItem(-1, 5, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum))

        self.frameLine2 = QtWidgets.QFrame()
        self.frameLine2.setFrameShape(QtWidgets.QFrame.HLine)
        self.frameLine2.setFrameShadow(QtWidgets.QFrame.Plain)
        self.vlayout_main.addWidget(self.frameLine2)

        self.vlayout_main.addItem(
            QtWidgets.QSpacerItem(-1, 5, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum))

        # buttons area
        self.hlayout_buttons = QtWidgets.QHBoxLayout()
        self.hlayout_buttons.setObjectName("hlayout_buttons")
        self.hlayout_buttons.setSpacing(5)
        self.lblMessage = QtWidgets.QLabel()
        self.lblMessage.setMinimumHeight(25)
        self.lblMessage.setMaximumHeight(25)
        self.lblMessage.setObjectName("lblMessage")
        self.hlayout_buttons.addWidget(self.lblMessage)
        self.btnRun = QtWidgets.QPushButton(self.mainContent)
        self.btnRun.setFixedSize(QtCore.QSize(90, 30))
        self.btnRun.setObjectName("btnRun")
        self.hlayout_buttons.addItem(
            QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum))
        self.hlayout_buttons.addWidget(self.btnRun)
        self.vlayout_main.addLayout(self.hlayout_buttons)

        self.controlArea.setContentsMargins(0,0,0,0)
        self.controlArea.layout().addLayout(self.vlayoutBase)
        self.controlArea.setMinimumSize(600, 500)

        # initialize UI components
        self.retranslateUi(self)

        # create UI events
        self.btnAddMapping.clicked.connect(self.OnAddVMapping)
        self.btnDeleteMapping.clicked.connect(self.OnDeleteVMapping)
        self.lstMapping.doubleClicked.connect(self.OnVMapDoubleClicked)
        self.btnRun.clicked.connect(self.OnRunContainer)


        self.InitializeUI()

    def retranslateUi(self, Widget):
        _translate = QtCore.QCoreApplication.translate
        self.lblFormTitle.setText(_translate("Widget", "Generic Task Runner"))
        self.lblDockerImage.setText(_translate("Widget", "Select docker image:"))
        self.lblVMapping.setText(_translate("Widget", "Mount Mapping: "))
        self.lblCmd.setText(_translate("Widget", "Run command:"))

        self.btnAddMapping.setText(_translate("Widget", "+"))
        self.btnDeleteMapping.setText(_translate("Widget", "-"))
        self.btnRun.setText(_translate("Widget", "Start"))

    def InitializeUI(self):
        self.loadDockerImages()

        self.model_vmap = QtGui.QStandardItemModel(self.lstMapping)
        self.model_vmap.setColumnCount(2)
        self.model_vmap.setHorizontalHeaderLabels(['From', 'To container'])
        self.lstMapping.setModel(self.model_vmap)
        width = self.lstMapping.width()
        self.lstMapping.setColumnWidth(0, 395)
        self.lstMapping.setColumnWidth(1, 170)


    def setHostMountedDir(self, directory, id):
        if os.path.exists(directory):
            items = self.model_vmap.findItems(directory)
            if not items:
                itemFrom = QStandardItem(directory)
                self.model_vmap.appendRow([itemFrom, QStandardItem()])

    def loadDockerImages(self):
        images = self.dockerClient.images()
        if not type(images) is list:
            self.cboDockerImage.addItem('No Docker Image', 'NULL')
        else:
            for img in images:
                if not img['RepoTags']: continue
                name = img['RepoTags'][0]
                if name == '<none>:<none>': continue
                id = img['Id'].split(':')[1][:12]
                self.cboDockerImage.addItem(name, id)

        self.cboDockerImage.view().setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)

    def OnAddVMapping(self):
        width = self.lstMapping.width() * 0.99
        self.lstMapping.setColumnWidth(0, width * 0.7)
        self.lstMapping.setColumnWidth(1, width * 0.3)
        self.model_vmap.appendRow([QStandardItem(), QStandardItem()])

    def OnDeleteVMapping(self):
        current = self.lstMapping.currentIndex()
        if current:
            self.model_vmap.removeRow(current.row())

    def OnVMapDoubleClicked(self, mode_index):
        column = mode_index.column()
        if column == 0:
            item = self.model_vmap.itemFromIndex(mode_index)
            if not item.text():
                dir = QtWidgets.QFileDialog.getExistingDirectory(self)
                item.setText(dir)

    def _enableUIElements(self, enabled=True):
        self.btnAddMapping.setEnabled(enabled)
        self.btnRun.setEnabled(enabled)
        self.btnDeleteMapping.setEnabled(enabled)
        self.cboDockerImage.setEnabled(enabled)
        self.lstMapping.setEnabled(enabled)
        self.txtCommand.setReadOnly(not enabled)

    def OnRunContainer(self):
        #disable UI elements
        self._enableUIElements(False)

        imageName = self.cboDockerImage.currentText()
        # get dir mapping
        volumes = {}
        for row in range(self.model_vmap.rowCount()):
            sFrom = self.model_vmap.data(self.model_vmap.index(row, 0))
            sTo = self.model_vmap.data(self.model_vmap.index(row, 1))
            if not os.path.exists(sFrom): continue
            if not sTo: continue
            volumes[sFrom] = sTo

        commands = self.txtCommand.toPlainText()
        #print(commands)

        self.lblMessage.setText('Running ' + imageName + ' ...')
        self.setStatusMessage('Running...')

        # Run the container in a new thread
        self.running_thread = GenericDockerRunner(self.dockerClient, imageName, volumes, commands)
        self.running_thread.finished.connect(self.ThreadEvent_OnRunCompleted)
        self.running_thread.start()

    def ThreadEvent_OnRunCompleted(self):
        self.lblMessage.setText("")
        self.setStatusMessage('Finished!')
        # restore controls
        self._enableUIElements()
        # notify next procedure
        #self.send("Output", self.fastqDirectory)


class GenericDockerRunner(QThread):
    def __init__(self, cli, image_name, volumes, commands):
        QThread.__init__(self)
        self.docker = cli
        self.image_name = image_name
        self.volumes = volumes
        self.commands = commands

    def __del__(self):
        self.wait()

    def run(self):
        response = self.docker.create_container(self.image_name, volumes=self.volumes, commands=self.commands)
        if response['Warnings'] == None:
            self.containerId = response['Id']
            self.docker.start_container(self.containerId)
        else:
            print(response['Warnings'])
        i = 1
        # Keep running until container is exited
        while self.docker.container_running(self.containerId):
            self.sleep(1)
        # Remove the container now that it is finished
        self.docker.remove_container(self.containerId)


def main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv))

    ow = OWGenericTask()
    ow.setHostMountedDir("/Users/Jimmy/Downloads/STAR/fasta_input", None)
    ow.show()
    ow.raise_()

    ow.handleNewSignals()
    app.exec_()
    return 0

if __name__=="__main__":
    sys.exit(main())