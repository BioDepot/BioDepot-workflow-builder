from Orange.widgets import widget, gui
import sys, os, fnmatch, tempfile, re
sys.path.append('/coreutils')
import requests, json
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import pyqtSlot, Qt, QThread, pyqtSignal
from PyQt5.QtGui import QStandardItem, QTextDocument, QFont
from DockerClient import DockerClient, DockerThread_BuildImage
from UIDockerfileEditor import DockerSyntaxHighlighter


class OWImageBuilder(widget.OWWidget):
    name = "Image Builder"
    description = "Build Custom Image"
    category = "Bwb-core"
    icon = "icons/imagebuilder.png"

    priority = 2

    inputs = []
    outputs = []#[("ImageTag", str)]

    want_main_area = True
    want_control_area = False

    def __init__(self):
        super().__init__()

        self.base_dir = os.path.dirname(os.path.abspath(__file__))
        self.FLAG_binder_string = '<BinderCompatible>'
        self.FLAG_binder_compatible = False

        # stylesheet
        css = '''
                QPushButton {background-color: #1588c5; color: white; height: 25px; border: 1px solid #1a8ac6; border-radius: 2px;}
                QPushButton:pressed { background-color: #158805; border-style: inset;}
                QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
                QPushButton:hover {background-color: #1588f5; }   
                QFrame#frameTitle {background: #1588c5;color: #1588c5}
                QFrame#frameSubtitle {background: #1998de; }
                QFrame#framePackages{border: 1px solid #1588c5;border-radius: 2px;}
                QFrame#frameSearchName {background: #1998de; }
                QLabel#lblInfoTitle{ 
                    background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 #1998de, stop:1 #1588c5);
                    padding-left: 3px;color: white;
                }
                QLineEdit#edtPackageName{
                    background-image: url(icons/search_icon.png);
                    background-repeat: no-repeat;
                    background-position: left;
                    border: 1px solid #1a8ac6; border-radius: 10px;
                    padding: 2 2 2 20;
                }
            '''
        # GUI
        self.setStyleSheet(css)
        fontTitle = QtGui.QFont()
        fontTitle.setPointSize(18) if sys.platform == 'darwin' else fontTitle.setPointSize(12)
        fontSubtitle = QtGui.QFont()
        fontSubtitle.setPointSize(14) if sys.platform == 'darwin' else fontTitle.setPointSize(11)
        fontSmallTitle = QtGui.QFont()
        fontSmallTitle.setPointSize(11) if sys.platform == 'darwin' else fontSmallTitle.setPointSize(8)
        self.fontDockerfileEditor = QtGui.QFont()
        self.fontDockerfileEditor.setFamily("Helvetica Neue" if sys.platform == "darwin" else "Courier New")
        self.fontDockerfileEditor.setPointSize(
            12) if sys.platform == 'darwin' else self.fontDockerfileEditor.setPointSize(10)

        # creat controls
        self.vlayoutBase = QtWidgets.QVBoxLayout()
        self.vlayoutBase.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.vlayoutBase.setContentsMargins(0, 0, 0, 0)
        self.vlayoutBase.setSpacing(0)
        self.vlayoutBase.setObjectName("vlayoutBase")
        # Title area
        self.frameTitle = QtWidgets.QFrame(self)
        self.frameTitle.setMinimumSize(QtCore.QSize(0, 65))
        self.frameTitle.setMaximumSize(QtCore.QSize(16777215, 65))
        self.frameTitle.setFrameShape(QtWidgets.QFrame.Box)
        self.frameTitle.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frameTitle.setObjectName("frameTitle")
        self.lblFormTitle = QtWidgets.QLabel()
        self.lblFormTitle.setFixedSize(QtCore.QSize(300, 30))
        self.lblFormTitle.setFont(fontTitle)
        self.lblFormTitle.setStyleSheet("color: white;")
        self.lblFormTitle.setObjectName("lblFormTitle")
        self.lblIcon = QtWidgets.QLabel()
        self.lblIcon.setFixedSize(QtCore.QSize(41, 31))
        self.lblIcon.setText("")
        self.lblIcon.setPixmap(QtGui.QPixmap(os.path.join(self.base_dir, "icons/builder_icon.png")))
        self.lblIcon.setScaledContents(True)
        self.lblIcon.setObjectName("lblIcon")
        self.lblDockerVersion = QtWidgets.QLabel()
        self.lblDockerVersion.setFixedSize(QtCore.QSize(220, 16))
        self.lblDockerVersion.setFont(fontSmallTitle)
        self.lblDockerVersion.setStyleSheet("color: white;")
        self.lblDockerVersion.setObjectName("lblDockerVersion")
        self.lblBuiding = QtWidgets.QLabel()
        self.lblBuiding.setFixedSize(QtCore.QSize(500, 16))
        self.lblBuiding.setFont(fontSmallTitle)
        self.lblBuiding.setStyleSheet("color: white;")
        self.lblBuiding.setObjectName("lblBuilding")
        self.lblBuiding.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.lblBuiding.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.lblBuidingStep = QtWidgets.QLabel()
        self.lblBuidingStep.setFixedSize(QtCore.QSize(300, 16))
        self.lblBuidingStep.setFont(fontSmallTitle)
        self.lblBuidingStep.setStyleSheet("color: white;")
        self.lblBuidingStep.setObjectName("lblBuidingStep")
        self.lblBuidingStep.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        # grid layout for title area
        self.glayout_title = QtWidgets.QGridLayout(self.frameTitle)
        self.glayout_title.addWidget(self.lblIcon, 0, 0, 2, 0)  # rowspan=2
        self.glayout_title.addWidget(self.lblFormTitle, 0, 1, Qt.AlignLeft | Qt.AlignTop)
        self.glayout_title.addWidget(self.lblBuidingStep, 0, 2, Qt.AlignRight | Qt.AlignTop)
        self.glayout_title.addWidget(self.lblDockerVersion, 1, 1, Qt.AlignLeft | Qt.AlignBottom)
        self.glayout_title.addWidget(self.lblBuiding, 1, 2, Qt.AlignRight | Qt.AlignBottom)
        self.glayout_title.setColumnMinimumWidth(0, 41)
        self.frameTitle.setLayout(self.glayout_title)
        self.vlayoutBase.addWidget(self.frameTitle)

        self.frameSubtitle = QtWidgets.QFrame(self)
        self.frameSubtitle.setEnabled(True)
        self.frameSubtitle.setMinimumSize(QtCore.QSize(0, 35))
        self.frameSubtitle.setMaximumSize(QtCore.QSize(16777215, 35))
        self.frameSubtitle.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frameSubtitle.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frameSubtitle.setObjectName("frameSubtitle")
        self.lblSubtitle = QtWidgets.QLabel(self.frameSubtitle)
        self.lblSubtitle.setGeometry(QtCore.QRect(15, 2, 400, 30))
        self.lblSubtitle.setFont(fontSubtitle)
        self.lblSubtitle.setStyleSheet("color: white")
        self.lblSubtitle.setObjectName("lblSubtitle")
        self.vlayoutBase.addWidget(self.frameSubtitle)
        # Left of main area
        self.mainContent = QtWidgets.QWidget(self)
        self.mainContent.setObjectName("mainContent")
        self.hlayout_mainArea = QtWidgets.QHBoxLayout(self.mainContent)
        self.hlayout_mainArea.setContentsMargins(8, 5, 8, 8)
        self.hlayout_mainArea.setSpacing(2)
        self.hlayout_mainArea.setObjectName("hlayout_mainArea")
        self.vlayout_content = QtWidgets.QVBoxLayout()
        self.vlayout_content.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.vlayout_content.setSpacing(5)
        self.vlayout_content.setContentsMargins(1, 1, 1, 1)
        self.vlayout_content.setObjectName("vlayout_content")
        self.lblName = QtWidgets.QLabel(self.mainContent)
        self.lblName.setMinimumSize(QtCore.QSize(0, 25))
        self.lblName.setMaximumSize(QtCore.QSize(16777215, 25))
        self.lblName.setObjectName("lblName")
        self.vlayout_content.addWidget(self.lblName)
        self.edtImageName = QtWidgets.QLineEdit(self.mainContent)
        self.edtImageName.setObjectName("edtImageName")
        self.vlayout_content.addWidget(self.edtImageName)
        self.lblBaseImage = QtWidgets.QLabel(self.mainContent)
        self.lblBaseImage.setMinimumSize(QtCore.QSize(0, 25))
        self.lblBaseImage.setMaximumSize(QtCore.QSize(16777215, 25))
        self.lblBaseImage.setObjectName("lblBaseImage")
        self.vlayout_content.addWidget(self.lblBaseImage)
        self.cboBaseImage = QtWidgets.QComboBox(self.mainContent)
        self.cboBaseImage.setMinimumSize(QtCore.QSize(0, 30))
        self.cboBaseImage.setMaximumSize(QtCore.QSize(16777215, 30))
        self.cboBaseImage.setObjectName("cboBaseImage")
        self.vlayout_content.addWidget(self.cboBaseImage)
        # self.lblRScript = QtWidgets.QLabel(self.mainContent)
        # self.lblRScript.setMinimumSize(QtCore.QSize(0, 25))
        # self.lblRScript.setMaximumSize(QtCore.QSize(16777215, 25))
        # self.lblRScript.setObjectName("lblRScript")
        # self.vlayout_content.addWidget(self.lblRScript)
        # self.hlayout_RScript = QtWidgets.QHBoxLayout()
        # self.hlayout_RScript.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        # self.hlayout_RScript.setContentsMargins(0, -1, -1, -1)
        # self.hlayout_RScript.setSpacing(1)
        # self.hlayout_RScript.setObjectName("horizontalLayout")
        # self.edtRScript = QtWidgets.QLineEdit(self.mainContent)
        # self.edtRScript.setObjectName("edtRScript")
        # self.edtRScript.setReadOnly(True)
        # self.hlayout_RScript.addWidget(self.edtRScript)
        # self.btnSelectScriptFile = QtWidgets.QPushButton(self.mainContent)
        # self.btnSelectScriptFile.setMinimumSize(QtCore.QSize(24, 24))
        # self.btnSelectScriptFile.setMaximumSize(QtCore.QSize(24, 24))
        # self.btnSelectScriptFile.setObjectName("btnSelectScriptFile")
        # self.hlayout_RScript.addWidget(self.btnSelectScriptFile)
        # self.vlayout_content.addLayout(self.hlayout_RScript)

        self.hlayout_dockertitle = QtWidgets.QHBoxLayout()
        self.hlayout_dockertitle.setObjectName("hlayout_dockertitle")
        self.hlayout_dockertitle.setContentsMargins(0, -1, -1, -1)
        self.hlayout_dockertitle.setSpacing(1)

        self.lblDockerfile = QtWidgets.QLabel(self.mainContent)
        self.lblDockerfile.setObjectName("lblDockerfile")
        # self.chkBinderCompatible = QtWidgets.QCheckBox()

        self.hlayout_dockertitle.addWidget(self.lblDockerfile)
        self.hlayout_dockertitle.addItem(
            QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum))
        # self.hlayout_dockertitle.addWidget(self.chkBinderCompatible)
        self.vlayout_content.addLayout(self.hlayout_dockertitle)
        # self.vlayout_content.addWidget(self.lblDockerfile)
        self.txtDockerfile = QtWidgets.QPlainTextEdit(self.mainContent)
        self.txtDockerfile.setObjectName("txtDockerfile")
        # self.txtDockerfile.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)      ### text wrap
        self.vlayout_content.addWidget(self.txtDockerfile)
        self.pbrBuildPrgoress = QtWidgets.QProgressBar(self.mainContent)
        self.pbrBuildPrgoress.setProperty("value", 0)
        self.pbrBuildPrgoress.setVisible(False)
        self.pbrBuildPrgoress.setObjectName("pbrBuildPrgoress")
        self.vlayout_content.addWidget(self.pbrBuildPrgoress)
        # buttons area
        self.hlayout_buttons = QtWidgets.QHBoxLayout()
        self.hlayout_buttons.setObjectName("hlayout_buttons")
        self.hlayout_buttons.setSpacing(5)
        self.btnOpen = QtWidgets.QPushButton(self.mainContent)
        self.btnOpen.setFixedSize(QtCore.QSize(75, 30))
        self.btnOpen.setObjectName("btnOpen")
        self.btnSave = QtWidgets.QPushButton(self.mainContent)
        self.btnSave.setFixedSize(QtCore.QSize(75, 30))
        self.btnSave.setObjectName("btnSave")
        self.btnBuild = QtWidgets.QPushButton(self.mainContent)
        self.btnBuild.setFixedSize(QtCore.QSize(90, 30))
        self.btnBuild.setObjectName("btnBuild")
        self.hlayout_buttons.addWidget(self.btnOpen)
        self.hlayout_buttons.addWidget(self.btnSave)
        self.hlayout_buttons.addItem(
            QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum))
        self.hlayout_buttons.addWidget(self.btnBuild)
        self.vlayout_content.addLayout(self.hlayout_buttons)
        # Right bioconductor packages aera
        self.vlayout_packages = QtWidgets.QVBoxLayout()
        self.vlayout_packages.setSpacing(0)
        self.vlayout_packages.setObjectName("vlayout_packages")
        self.vlayout_packages.setContentsMargins(0, 0, 0, 0)
        self.lblInfoTitle = QtWidgets.QLabel(self.mainContent)
        self.lblInfoTitle.setMinimumSize(QtCore.QSize(0, 30))
        self.lblInfoTitle.setMaximumSize(QtCore.QSize(16777215, 30))
        self.lblInfoTitle.setObjectName("lblInfoTitle")
        self.vlayout_packages.addWidget(self.lblInfoTitle)

        self.frameSearchName = QtWidgets.QFrame()
        self.frameSearchName.setContentsMargins(5, 5, 2, 5)
        self.frameSearchName.setObjectName("frameSearchName")
        self.frameSearchName.setMinimumHeight(35)
        self.frameSearchName.setMaximumHeight(35)
        self.hlayout_list_search = QtWidgets.QHBoxLayout()
        self.hlayout_list_search.setObjectName("hlayout_list_search")
        self.hlayout_list_search.setContentsMargins(1, 0, 0, 2)
        self.edtPackageName = QtWidgets.QLineEdit(self.mainContent)
        self.edtPackageName.setClearButtonEnabled(True)
        self.edtPackageName.setObjectName("edtPackageName")
        self.edtPackageName.setPlaceholderText("Search package")
        self.edtPackageName.setMinimumHeight(25)
        self.edtPackageName.setMaximumHeight(25)
        self.hlayout_list_search.addWidget(self.edtPackageName)
        self.frameSearchName.setLayout(self.hlayout_list_search)
        self.vlayout_packages.addWidget(self.frameSearchName)
        self.lstPackages = QtWidgets.QTableView(self.mainContent)
        self.lstPackages.setObjectName("lstPackages")
        self.lstPackages.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.lstPackages.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.lstPackages.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.lstPackages.verticalHeader().setVisible(False)
        self.lstPackages.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.vlayout_packages.addWidget(self.lstPackages)

        self.framePackages = QtWidgets.QFrame(self)
        self.framePackages.setContentsMargins(0, 0, 0, 0)
        self.framePackages.setFrameShape(QtWidgets.QFrame.Box)
        self.framePackages.setFrameShadow(QtWidgets.QFrame.Plain)
        self.framePackages.setObjectName("framePackages")
        # self.framePackages.setLayout(self.vlayout_packages)

        # main window horizontal splitter
        container_mainLeft = QtWidgets.QWidget()
        container_mainRight = self.framePackages
        container_mainLeft.setLayout(self.vlayout_content)
        container_mainRight.setLayout(self.vlayout_packages)
        self.splitter_main = QtWidgets.QSplitter()
        self.splitter_main.addWidget(container_mainLeft)
        self.splitter_main.addWidget(container_mainRight)
        self.splitter_main.setOrientation(Qt.Horizontal)
        self.splitter_main.setSizes([200, 150])
        self.hlayout_mainArea.addWidget(self.splitter_main)

        self.vlayoutBase.addWidget(self.mainContent)

        self.mainArea.layout().addLayout(self.vlayoutBase)
        self.mainArea.layout().setContentsMargins(0, 0, 0, 0)

        #LHH commented out to add scrollbars 
        # set window size
        #width = 900
        #height = 700
        #self.mainArea.setMinimumSize(width, height)

        # initialize UI components
        self.retranslateUi(self)
        self.model_package = QtGui.QStandardItemModel(self.lstPackages)
        self.model_package.setColumnCount(2)
        headerNames = ['Package Name', 'Detail']
        self.model_package.setHorizontalHeaderLabels(headerNames)
        self.package_list_proxy = QtCore.QSortFilterProxyModel(self)
        self.package_list_proxy.setSourceModel(self.model_package)

        # create UI events
        #self.btnSelectScriptFile.clicked.connect(self.OnChooseScriptFile)
        self.btnBuild.clicked.connect(self.OnBuildClicked)
        self.btnOpen.clicked.connect(self.OnLoadDockerfile)
        self.btnSave.clicked.connect(self.OnSaveDockerfile)
        self.cboBaseImage.currentIndexChanged.connect(self.OnBaseImageSelectChanged)
        self.edtPackageName.textChanged.connect(self.OnPackageNameChanged)
        self.model_package.itemChanged.connect(self.OnPackageListSelectedChanged)
        
        #LHH patch to add scrollbars
        self.scroll_area = QtWidgets.QScrollArea(
            verticalScrollBarPolicy=Qt.ScrollBarAlwaysOn
        )
        self.scroll_area.setWidget(self.mainContent)
        self.scroll_area.setWidgetResizable(True)
        self.mainArea.layout().addWidget(self.scroll_area)
        self.InitializeUI()
        #self.show()

    def retranslateUi(self, Widget):
        _translate = QtCore.QCoreApplication.translate
        Widget.setWindowTitle(_translate("Widget", "Bioconductor Docker Image Builder"))
        self.lblFormTitle.setText(_translate("Widget", "Bioconductor Docker Image Builder"))
        self.lblSubtitle.setText(_translate("Widget", "Customize your own bioconductor image"))
        self.lblName.setText(_translate("Widget", "Image Name: "))
        self.lblBaseImage.setText(_translate("Widget", "Select base image:"))
        self.lblDockerfile.setText(_translate("Widget", "Builder detail: "))
        self.btnBuild.setText(_translate("Form", "Build"))
        self.lblDockerVersion.setText(_translate("Form", "Docker Engine: {0}"))
        # self.lblBuiding.setText(_translate("Form", 'Ready'))
        self.lblDockerfile.setText(_translate("Form", "Docker file:"))
        self.btnOpen.setText(_translate("Form", "Open"))
        self.btnSave.setText(_translate("Form", "Save"))
        self.btnBuild.setText(_translate("Form", "Build"))
        self.lblInfoTitle.setText(_translate("Form", "Bioconductor/CRAN R packages"))
        # self.chkBinderCompatible.setText(_translate("Form", "Binder Compatible"))

    def InitializeUI(self):
        # Init Docker Engine
        self.dockerInitialized = False
        self.dockerClient = None
        strDockerInfo = ""
        try:
            self.dockerClient = DockerClient('unix:///var/run/docker.sock', 'local')
            strDockerInfo = self.lblDockerVersion.text().format(self.dockerClient.version()['Version'])
            self.dockerInitialized = True
        except:
            e = sys.exc_info()[0]
            print(e)
            strDockerInfo = "No Docker Engine installed OR missing docker-py"
            self.btnBuild.setEnabled(False)

        self.SelectedBiocPackage = []
        self.lblDockerVersion.setText(strDockerInfo)
        # self.dockerClient.images()

        # scan docker files
        self.cboBaseImage.addItem("From scratch", "-SCRATCH-")
        self.dockerfile_dir = os.path.join(self.base_dir, 'DockerFiles')
        dockerfiles = [x for x in fnmatch.filter(os.listdir(self.dockerfile_dir), '*.Dockerfile')]
        for f in dockerfiles:
            filename, _ = os.path.splitext(os.path.basename(f))
            self.cboBaseImage.addItem(filename, os.path.join(self.dockerfile_dir, f))

        self.building_log_file = os.path.join(self.base_dir, 'building.log')
        self.LoadBiocPakageList()

    @pyqtSlot()
    def OnChooseScriptFile(self):
        start_file = os.path.expanduser("~/")
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open R Script File', start_file)
        if not filename:
            return

        self.edtRScript.setText(filename)

    def OnBaseImageSelectChanged(self, index):
        dockerfile = self.cboBaseImage.itemData(index)
        self.documentFromDockerfile(dockerfile)

    def documentFromDockerfile(self, filename):
        doc = QTextDocument(self)
        doc.setDocumentLayout(QtWidgets.QPlainTextDocumentLayout(doc))
        if filename == "-SCRATCH-":
            doc.setPlainText("")
        else:
            with open(filename, 'r') as f:
                doc.setPlainText(f.read())
        doc.setDefaultFont(self.fontDockerfileEditor)
        doc.highlighter = DockerSyntaxHighlighter(doc)
        # doc.modificationChanged[bool].connect(self.onModificationChanged)
        doc.setModified(False)
        self._cachedDocument = doc
        self.txtDockerfile.setDocument(doc)

        plainText = self.txtDockerfile.toPlainText()
        if self._find_bioclite(self.FLAG_binder_string, plainText):
            self.FLAG_binder_compatible = True
        else:
            self.FLAG_binder_compatible = False

        return self._cachedDocument

    def LoadBiocPakageList(self):
        # loading
        item = QStandardItem(" >>> Loading bioconductor package list <<< ")
        self.model_package.appendRow(item)
        self.lstPackages.setModel(self.package_list_proxy)
        self.lstPackages.resizeColumnsToContents()

        self.loadpackage_thread = BiocPackageList()
        self.loadpackage_thread.load_completed.connect(self.ThreadEvent_OnLoadBiocPackageCompleted)
        self.loadpackage_thread.start()

    def ThreadEvent_OnLoadBiocPackageCompleted(self, packagelist):
        self.model_package.removeRows(0, self.model_package.rowCount())
        for pkg in packagelist:
            itemName = QStandardItem(pkg['Name'])
            itemName.setCheckable(True)
            itemTitle = QStandardItem(pkg['Title'])
            self.model_package.appendRow([itemName, itemTitle])

        # Apply the model to the list view
        self.lstPackages.setModel(self.package_list_proxy)
        # self.lstPackages.resizeColumnToContents(0)
        self.lstPackages.setColumnWidth(0, 190)
        self.lstPackages.setColumnWidth(1, 130)

    @pyqtSlot(str)
    def OnPackageNameChanged(self, text):
        search = QtCore.QRegExp(text, QtCore.Qt.CaseInsensitive, QtCore.QRegExp.RegExp)
        self.package_list_proxy.setFilterRegExp(search)

    def _move_editor_cursor(self, start, end):
        cursor = self.txtDockerfile.textCursor()
        cursor.setPosition(start)
        cursor.movePosition(QtGui.QTextCursor.Right, QtGui.QTextCursor.KeepAnchor, end - start)
        self.txtDockerfile.setTextCursor(cursor)

    def _find_bioclite(self, query, text):
        pattern = re.compile(query)
        return pattern.search(text, 0) if query != "" else None

    def _update_bioc_package_in_dockerfile(self, previous_package):
        base_bioclite = "RUN Rscript -e \"source('https://bioconductor.org/biocLite.R');biocLite(c({0}),ask=FALSE)\"\n"
        if self.FLAG_binder_compatible:
            # if self.chkBinderCompatible.isChecked():
            base_bioclite = "RUN echo \"source('http://bioconductor.org/biocLite.R'); biocLite(c({0}))\" | R --vanilla\n"

        previous = ','.join("'{0}'".format(w) for w in previous_package)
        current = ','.join("'{0}'".format(w) for w in self.SelectedBiocPackage)

        # no package selected, delete
        delete_mode = previous_package and not self.SelectedBiocPackage
        if delete_mode: current = ''

        plainText = self.txtDockerfile.toPlainText()

        matched = self._find_bioclite(previous, plainText)

        if matched:
            self._move_editor_cursor(matched.start(), matched.end())
            if delete_mode:
                cursor = self.txtDockerfile.textCursor()
                cursor.movePosition(QtGui.QTextCursor.StartOfBlock)
                cursor.movePosition(QtGui.QTextCursor.EndOfBlock, QtGui.QTextCursor.KeepAnchor)
                self.txtDockerfile.setTextCursor(cursor)
        else:
            # try locate "CMD"
            matched_cmd = self._find_bioclite("CMD", plainText)
            if not matched_cmd: matched_cmd = self._find_bioclite("WORKDIR", plainText)
            if matched_cmd:
                self._move_editor_cursor(matched_cmd.start(), matched_cmd.end())
                cursor = self.txtDockerfile.textCursor()
                cursor.movePosition(QtGui.QTextCursor.StartOfBlock)
                self.txtDockerfile.setTextCursor(cursor)
            else:
                self.txtDockerfile.moveCursor(QtGui.QTextCursor.End)
                base_bioclite = "\n" + base_bioclite

        cursor = self.txtDockerfile.textCursor()
        if matched and cursor.hasSelection():
            cursor.insertText(current)
        else:
            cursor.insertText(base_bioclite.format(current))
        self.txtDockerfile.setTextCursor(cursor)

    def OnPackageListSelectedChanged(self, item):
        package_name = item.text()
        previous_package = self.SelectedBiocPackage.copy()
        if not item.checkState():
            self.SelectedBiocPackage.remove(package_name)
        else:
            self.SelectedBiocPackage.append(package_name)

        self._update_bioc_package_in_dockerfile(previous_package)

    def OnSaveDockerfile(self):
        filename = self.cboBaseImage.currentData()
        if filename == "-SCRATCH-":
            filename = os.path.expanduser("~/")

        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save Dockerfile',
            filename,
            'Docker files (*)'
        )

        if filename:
            f = open(filename, 'w')
            f.write(self.txtDockerfile.toPlainText())
            f.close()

    def OnLoadDockerfile(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open Dockerfile',
            os.path.expanduser("~/"),
            'Docker files (*)'
        )
        if filename:
            itemName = '... ' + os.path.join(os.path.basename(os.path.dirname(filename)), os.path.basename(filename))
            index = self.cboBaseImage.findData(filename)
            if index >= 0:
                self.cboBaseImage.setCurrentIndex(index)
            else:
                self.cboBaseImage.addItem(itemName, filename)
                self.cboBaseImage.setCurrentIndex(self.cboBaseImage.findData(filename))

    def _set_building_text(self, labelCtrl, text):
        metrics = QtGui.QFontMetrics(labelCtrl.font())
        elidedText = metrics.elidedText(text, Qt.ElideRight, labelCtrl.width())
        labelCtrl.setText(elidedText)

    def _enableUIElements(self, enabled=False):
        self.edtImageName.setEnabled(enabled)
        #self.btnSelectScriptFile.setEnabled(enabled)
        self.txtDockerfile.setReadOnly(not enabled)
        self.lstPackages.setEnabled(enabled)
        self.edtPackageName.setEnabled(enabled)
        self.btnOpen.setEnabled(enabled)
        self.btnSave.setEnabled(enabled)
        self.btnBuild.setEnabled(self.dockerInitialized and enabled)
        self.cboBaseImage.setEnabled(enabled)

    def _building_message_processor(self, message):
        # write to log file
        #     with open(self.building_log_file, 'a') as f:
        #         f.write(message)

        # eliminate blank lines
        message = [s for s in message.splitlines() if s.strip()]
        if not message:
            return
        # message = ''.join(message)
        message = message[0]
        self._set_building_text(self.lblBuiding, message)

        endOfBuiding = False
        # extract step progress
        #   example: 'Step 10/10 : FROM ubuntu:16.04'
        step_pattern = re.compile(r'Step \d+\/\d+')
        if step_pattern.search(message):
            progress = message.split(' : ')[0].split(' ')[1].split('/')
            progress = [int(x) for x in progress]
            if progress[0] == 1:
                self.pbrBuildPrgoress.setMaximum(progress[1])
                self.progressBarInit()

            self.pbrBuildPrgoress.setValue(progress[0])
            self.progressBarSet(progress[0]/progress[1]*100)
            self._set_building_text(self.lblBuidingStep, message)

            #endOfBuiding = progress[0] == progress[1]

        successful = re.search(r'Successfully built ([0-9a-f]+)', message)

        # errorDetail={'code': 2, 'message': "The command '/bin/sh -c pip3 install -r orange3/requirements-core.txt' returned a non-zero code: 2"}, error="The command '/bin/sh -c pip3 install -r orange3/requirements-core.txt' returned a non-zero code: 2"
        failure = message.find("errorDetail=") >= 0 and message.find("error=") >= 0

        if endOfBuiding or successful or failure:
            print(message)
            if not failure:
                self._set_building_text(self.lblBuidingStep, '')
                self.setStatusMessage('Build completed.')
            else:
                self.error('Built with errors.')
                self.setStatusMessage('Built with errors.')
            self._enableUIElements(True)
            self.pbrBuildPrgoress.setVisible(False)
            self.progressBarFinished()
            if self.dockerfile_for_build:
                if sys.platform != 'win32':
                    os.unlink(self.dockerfile_for_build)
                else:
                    os.remove(self.dockerfile_for_build)
                self.dockerfile_for_build = ''

    @pyqtSlot()
    def OnBuildClicked(self):
        imagename = self.edtImageName.text()
        self.clear_messages()
        if not imagename:
            self.error('No Image Name')
            # msg = QtWidgets.QMessageBox()
            # msg.setText('No Image Name')
            # msg.setInformativeText("Please specify a image name")
            # msg.setIcon(QtWidgets.QMessageBox.Warning)
            # msg.exec()
            self.edtImageName.setFocus()
            return

        # verify image name
        pattern = re.compile(
            "^(?:(?=[^:\/]{1,253})(?!-)[a-zA-Z0-9-]{1,63}(?<!-)(?:\.(?!-)[a-zA-Z0-9-]{1,63}(?<!-))*(?::[0-9]{1,5})?/)?((?![._-])(?:[a-z0-9._-]*)(?<![._-])(?:/(?![._-])[a-z0-9._-]*(?<![._-]))*)(?::(?![.-])[a-zA-Z0-9_.-]{1,128})?$")
        if not pattern.search(imagename):
            self.error('Invalid image name')
            msg = QtWidgets.QMessageBox()
            msg.setText('Invalid image name')
            msg.setInformativeText("Typical image name:\n    registry/image-name[:version] \n\n"
                                   "For example: \n    biodepot/bwb:latest")
            msg.setIcon(QtWidgets.QMessageBox.Warning)
            #msg.exec()
            self.edtImageName.setFocus()
            return

        dockerfile = self.txtDockerfile.toPlainText()

        # skip empty docker file
        if not dockerfile.strip():
            self._set_building_text(self.lblBuiding, 'Error: empty Docker file')
            self.warning('Error: empty Docker file')
            return

        self._enableUIElements(False)
        self.pbrBuildPrgoress.setVisible(True)
        self.setStatusMessage('Building...')

        buildpath = self.cboBaseImage.currentData()
        if buildpath == "-SCRATCH-":
            buildpath = self.base_dir
        else:
            buildpath = os.path.dirname(buildpath)
            if buildpath == self.dockerfile_dir:
                buildpath = self.base_dir

        # create a temp shadow dockerfile for building
        shadow_dockerfile = ''
        if sys.platform == 'win32':
            shadow_dockerfile = 'Dockerfile.build.tmp'
            self.dockerfile_for_build = os.path.join(buildpath, 'Dockerfile.build.tmp')
            with open(self.dockerfile_for_build, 'wb') as fp:
                fp.write(dockerfile.encode('utf-8'))
        else:
            fp = tempfile.NamedTemporaryFile(delete=False)
            fp.write(dockerfile.encode('utf-8'))
            fp.close()
            self.dockerfile_for_build = fp.name
            shadow_dockerfile = fp.name

        self._set_building_text(self.lblBuidingStep, "Start building")
        # Call docker API to build image using thread
        self.buildimage_thread = DockerThread_BuildImage(self.dockerClient, imagename, buildpath, shadow_dockerfile)
        self.buildimage_thread.build_process.connect(self.ThreadEvent_OnImageBuilding)
        self.buildimage_thread.build_complete.connect(self.ThreadEvent_OnImageBuildComplete)
        self.buildimage_thread.start()

    def ThreadEvent_OnImageBuilding(self, message):
        # filter message, eliminate terminal colors
        ESC = r'\x1b'
        CSI = ESC + r'\['
        OSC = ESC + r'\]'
        CMD = '[@-~]'
        ST = ESC + r'\\'
        BEL = r'\x07'
        pattern = '(' + CSI + '.*?' + CMD + '|' + OSC + '.*?' + '(' + ST + '|' + BEL + ')' + ')'
        plainMessage = re.sub(pattern, '', message)
        self._building_message_processor(plainMessage)

    def ThreadEvent_OnImageBuildComplete(self):
        self._set_building_text(self.lblBuidingStep, '')
        self._enableUIElements(True)
        self.pbrBuildPrgoress.setVisible(False)
        self.progressBarFinished()


from bs4 import BeautifulSoup


class BiocPackageList(QThread):
    load_completed = pyqtSignal(object)
    # fetch_url = "http://bioconductor.org/packages/release/BiocViews.html#___Software"
    package_json_url = [
        "https://bioconductor.org/packages/json/3.5/bioc/packages.js",
        "https://bioconductor.org/packages/json/3.5/data/experiment/packages.js",
        "https://bioconductor.org/packages/json/3.5/data/annotation/packages.js"
    ]

    cran_packages_url = "https://cran.r-project.org/web/packages/available_packages_by_name.html"

    def __init__(self):
        QThread.__init__(self)

    def __del__(self):
        self.wait()

    def run(self):
        package_list = []
        # load package list from CRAN
        try:
            cranhtml = requests.get(self.cran_packages_url)
            soup = BeautifulSoup(cranhtml.text, "html.parser")
            for table in soup.findAll('table'):
                for row in table.findAll('tr'):
                    cells = row.findAll('td')
                    if len(cells) == 2:
                        links = cells[0].findAll("a")
                        if len(links) > 0:
                            name = links[0].contents[0]
                        title = cells[1].string
                        package_list.append({"Name": name, "Title": title})
        except Exception as e:
            print(str(e))

        # load package list from bioconductor
        try:
            for url in self.package_json_url:
                rawhtml = requests.get(url)
                jsonValue = '{%s}' % (rawhtml.text.split('{', 1)[1].rsplit('}', 1)[0],)
                packages = json.loads(jsonValue)
                for item in packages['content']:
                    package_list.append({"Name": item[0], "Title": item[2]})
        except Exception as e:
            print(str(e))

        package_list = sorted(package_list, key=lambda x: x['Name'].lower())
        # print(len(package_list))

        self.load_completed.emit(package_list)
