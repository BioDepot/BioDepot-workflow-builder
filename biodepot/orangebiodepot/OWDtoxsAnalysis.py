import os, fnmatch
import Orange.data
from Orange.data.io import FileFormat
from Orange.widgets import widget, gui, settings
from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5 import QtGui, QtWidgets, QtCore
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread

class OWDtoxsAnalysis(widget.OWWidget):
    name = "Dtoxs Analysis"
    description = "Step 2 of Dtoxs SOP. Uses edgeR for differential expression analysis."
    category = "RNASeq"
    icon = "icons/dtoxs-analysis2.svg"
    priority = 10

    image_name = "biodepot/dtoxs_analysis"
    image_version = "latest"

    inputs = [("Counts", str, "set_aligns", widget.Default),
              ("Configs", str, "set_configs"),
              ("Params", str, "set_params")]
    outputs = [("Results", str), ("Top 40", Orange.data.Table)]

    want_main_area = True
    want_control_area = False

    Exp_Design_file = settings.Setting('', schema_only=True)

    def __init__(self):
        super().__init__()

        self.result_folder_name = "Results"
        self.host_counts_dir = None
        self.host_config_dir = None
        self.host_param_dir = None

        # This client talks to your local docker
        self.docker = DockerClient('unix:///var/run/docker.sock', 'local')

        # GUI
        self.setStyleSheet("QPushButton{\n"
                           "    background-color: #1588c5;\n"
                           "    color: white;\n"
                           "    height: 25px;\n"
                           "    border: 1px solid #1a8ac6;\n"
                           "}")
        # creat controls
        self.vlayoutBase = QtWidgets.QVBoxLayout()
        self.vlayoutBase.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.vlayoutBase.setContentsMargins(0, 0, 0, 0)
        self.vlayoutBase.setSpacing(0)
        self.vlayoutBase.setObjectName("vlayoutBase")

        self.vlayoutBase = QtWidgets.QVBoxLayout()
        self.vlayoutBase.setContentsMargins(0, 0, 0, 0)
        self.vlayoutBase.setSpacing(0)
        self.vlayoutBase.setObjectName("vlayoutBase")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.frameTitle = QtWidgets.QFrame(self.mainArea)
        self.frameTitle.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frameTitle.sizePolicy().hasHeightForWidth())
        self.frameTitle.setSizePolicy(sizePolicy)
        self.frameTitle.setMinimumSize(QtCore.QSize(0, 45))
        self.frameTitle.setMaximumSize(QtCore.QSize(16777215, 65))
        self.frameTitle.setStyleSheet("QFrame#frameTitle{\n"
                                      "    background: #1588c5;\n"
                                      "    color: #1588c5\n"
                                      "}")
        self.frameTitle.setFrameShape(QtWidgets.QFrame.Box)
        self.frameTitle.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frameTitle.setObjectName("frameTitle")
        self.lblFormTitle = QtWidgets.QLabel(self.frameTitle)
        self.lblFormTitle.setGeometry(QtCore.QRect(10, 13, 221, 21))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setWeight(50)
        self.lblFormTitle.setFont(font)
        self.lblFormTitle.setStyleSheet("\n"
                                        "    color: white;\n"
                                        "")
        self.lblFormTitle.setObjectName("lblFormTitle")
        self.verticalLayout.addWidget(self.frameTitle)
        self.mainContent = QtWidgets.QWidget(self.mainArea)
        self.mainContent.setObjectName("mainContent")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.mainContent)
        self.verticalLayout_4.setContentsMargins(15, 15, 15, 15)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.vlayout_content = QtWidgets.QVBoxLayout()
        self.vlayout_content.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.vlayout_content.setSpacing(2)
        self.vlayout_content.setObjectName("vlayout_content")
        self.lblExpFile = QtWidgets.QLabel(self.mainContent)
        self.lblExpFile.setMinimumSize(QtCore.QSize(0, 25))
        self.lblExpFile.setMaximumSize(QtCore.QSize(16777215, 25))
        self.lblExpFile.setObjectName("lblExpFile")
        self.vlayout_content.addWidget(self.lblExpFile)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.horizontalLayout.setContentsMargins(0, -1, -1, -1)
        self.horizontalLayout.setSpacing(1)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.edtExpDesignFile = QtWidgets.QLineEdit(self.mainContent)
        self.edtExpDesignFile.setObjectName("edtExpDesignFile")
        self.horizontalLayout.addWidget(self.edtExpDesignFile)
        self.btnSelectExpFile = QtWidgets.QPushButton(self.mainContent)
        self.btnSelectExpFile.setMinimumSize(QtCore.QSize(24, 24))
        self.btnSelectExpFile.setMaximumSize(QtCore.QSize(24, 24))
        self.btnSelectExpFile.setObjectName("btnSelectExpFile")
        self.horizontalLayout.addWidget(self.btnSelectExpFile)
        self.vlayout_content.addLayout(self.horizontalLayout)
        self.line_2 = QtWidgets.QFrame(self.mainContent)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.vlayout_content.addWidget(self.line_2)
        self.lblInfoTitle = QtWidgets.QLabel(self.mainContent)
        self.lblInfoTitle.setMinimumSize(QtCore.QSize(0, 25))
        self.lblInfoTitle.setMaximumSize(QtCore.QSize(16777215, 25))
        self.lblInfoTitle.setStyleSheet(
            "background-color:qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 rgba(235, 235, 235, 255), stop:1 rgba(217, 217, 217, 255));\n"
            "padding-left: 3px;")
        self.lblInfoTitle.setObjectName("lblInfoTitle")
        self.vlayout_content.addWidget(self.lblInfoTitle)
        self.line = QtWidgets.QFrame(self.mainContent)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.vlayout_content.addWidget(self.line)
        self.infoLabel = QtWidgets.QLabel(self.mainContent)
        self.infoLabel.setMinimumSize(QtCore.QSize(0, 45))
        self.infoLabel.setStyleSheet("padding-left: 3px")
        self.infoLabel.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.infoLabel.setAlignment(QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
        self.infoLabel.setWordWrap(True)
        self.infoLabel.setObjectName("infoLabel")
        self.vlayout_content.addWidget(self.infoLabel)
        self.btn_run = QtWidgets.QPushButton(self.mainContent)
        self.btn_run.setMinimumSize(QtCore.QSize(90, 0))
        self.btn_run.setMaximumSize(QtCore.QSize(90, 16777215))
        self.btn_run.setObjectName("btn_run")
        self.vlayout_content.addWidget(self.btn_run, 0, QtCore.Qt.AlignRight)
        self.verticalLayout_4.addLayout(self.vlayout_content)
        self.verticalLayout.addWidget(self.mainContent)
        self.vlayoutBase.addLayout(self.verticalLayout)

        self.mainArea.layout().addLayout(self.vlayoutBase)
        self.mainArea.layout().setContentsMargins(0, 0, 0, 0)
        self.mainArea.setMinimumSize(568, 246)

        # events
        self.btnSelectExpFile.clicked.connect(self.OnChooseExpFile)
        self.btn_run.clicked.connect(self.start_analysis)

        self.retranslateUi(self.mainArea)


    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Widget", "Form"))
        self.lblFormTitle.setText(_translate("Widget", "DTOXS Analysis"))
        self.lblExpFile.setText(_translate("Widget", "Select experiment design file:"))
        self.btnSelectExpFile.setText(_translate("Widget", " ☰ "))
        self.lblInfoTitle.setText(_translate("Widget", "Information"))
        self.infoLabel.setText(_translate("Widget",
                                          "Connect to Dtoxs Alignment or use a Directory widget to specify location of UMI read counts data files"))
        self.btn_run.setText(_translate("Widget", "Run"))
        self.edtExpDesignFile.setText(self.Exp_Design_file)

    """
    Set input
    """
    def set_aligns(self, path):
        if not type(path) is str:
            # TODO create warning
            print('Tried to set alignment directory to None')
        elif not os.path.exists(path):
            # TODO create warning
            print('Tried to set alignment directory to none existent directory: ' + str(path))
        else:
            self.host_counts_dir = path.strip()
            # Jimmy March-29-2017, once the counts input was set, automatically create an output fold as a sibling of "Seqs"
            parent_path = os.path.abspath(os.path.join(self.host_counts_dir, '..'))
            self.host_results_dir = os.path.join(parent_path, self.result_folder_name)
            if not os.path.exists(self.host_results_dir):
                os.makedirs(self.host_results_dir)

            if os.path.exists(self.Exp_Design_file):
                import shutil
                shutil.copy2(self.Exp_Design_file, self.host_counts_dir)

            if self.host_param_dir and self.host_config_dir and self.Exp_Design_file:
                self.infoLabel.setText('All set.\nWaiting to run...')
                #run by default when all is set - TODO add a checkbox to this
                self.start_analysis()

    def set_configs(self, path):
        if not type(path) is str:
             print('Tried to set configs directory to None')
        elif not os.path.exists(path):
            print('Tried to set configs directory to none existent directory: ' + str(path))
        else:
            self.host_config_dir = path.strip()
            if self.host_counts_dir and self.host_param_dir and self.Exp_Design_file:
                self.infoLabel.setText('All set.\nWaiting to run...')

    def set_params(self, path):
        if not type(path) is str:
             print('Tried to set params directory to None')
        elif not os.path.exists(path):
            print('Tried to set params directory to none existent directory: ' + str(path))
        else:
            self.host_param_dir = path.strip()
            if self.host_counts_dir and self.host_config_dir and self.Exp_Design_file:
                self.infoLabel.setText('All set.\nWaiting to run...')

    def OnChooseExpFile(self):
        start_file = os.path.expanduser("~/")
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open Experiment Design File', start_file)
        if not filename:
            return

        self.edtExpDesignFile.setText(filename)
        self.Exp_Design_file = filename
        if self.host_counts_dir and self.host_config_dir and self.host_param_dir:
            import shutil
            shutil.copy2(self.Exp_Design_file, self.host_counts_dir)
            self.infoLabel.setText('All set.\nWaiting to run...')

    """
    Pull image
    """
    def pull_image(self):
        self.infoLabel.setText('Pulling \'' + self.image_name + ":" + self.image_version + '\' from Dockerhub...')
        self.setStatusMessage("Downloading...")
        self.progressBarInit()
        self.btn_run.setEnabled(False)
        # Pull the image in a new thread
        self.pull_image_thread = PullImageThread(self.docker, self.image_name, self.image_version)
        self.pull_image_thread.pull_progress.connect(self.pull_image_progress)
        self.pull_image_thread.finished.connect(self.pull_image_done)
        self.pull_image_thread.start()

    def pull_image_progress(self, val):
        self.progressBarSet(val)

    def pull_image_done(self):
        self.infoLabel.setText('Finished pulling \'' + self.image_name + ":" + self.image_version + '\'')
        self.progressBarFinished()
        self.run_analysis()

    """
    Analysis
    """
    def start_analysis(self):
        # Make sure the docker image is downloaded
        if not self.docker.has_image(self.image_name, self.image_version):
            self.pull_image()
        # Make sure there is an alignment directory set
        elif self.host_counts_dir and self.host_config_dir and self.host_param_dir:
            self.run_analysis()
        else:
            self.infoLabel.setText('Set counts, params and configs directories before running.')

    def run_analysis(self):
        self.infoLabel.setText('Running analysis...')
        self.setStatusMessage('Running...')
        #self.progressBarInit()
        # Run the container in a new thread
        self.run_analysis_thread = RunAnalysisThread(self.docker,
                                                     self.image_name,
                                                     self.host_counts_dir,
                                                     self.host_results_dir,
                                                     self.host_config_dir,
                                                     self.host_param_dir)

        self.run_analysis_thread.analysis_progress.connect(self.run_analysis_progress)
        self.run_analysis_thread.finished.connect(self.run_analysis_done)
        self.run_analysis_thread.start()

    def run_analysis_progress(self, val):
        self.progressBarSet(val)

    def run_analysis_done(self):
        self.infoLabel.setText("Finished running analysis!")
        self.btn_run.setEnabled(True)
        self.btn_run.setText('Run again')
        self.setStatusMessage('Finished!')
        # self.progressBarFinished()
        self.send("Results", self.host_results_dir)

        # Jimmy March 29 2017 added, create a dataset for DataTable widget to show TOP-40.tsv
        tsvFile = os.path.join(self.host_results_dir, 'FDR-0.1/TOP-40.tsv');
        tsvReader = FileFormat.get_reader(tsvFile)
        data = None
        try:
            data = tsvReader.read()

            # TODO Jimmy March 29 2017, tried to define domain for DataTable, didn't work
            '''
            domain = Orange.data.Domain(["Cell", "Plate", "Gene", "logFC", "logCPM", "PValue"], data.domain)

            domain = Orange.data.Domain(data.domain.attributes, data.domain.class_vars,
                [Orange.data.DiscreteVariable("Cell"),
                Orange.data.DiscreteVariable("Plate"),
                Orange.data.DiscreteVariable("Condition1"),
                Orange.data.DiscreteVariable("Condition2"),
                Orange.data.StringVariable("Gene"),
                Orange.data.ContinuousVariable("logFC"),
                Orange.data.ContinuousVariable("logCPM"),
                Orange.data.StringVariable("PValue")])

            data = data.from_table(domain, data)'''
        except Exception as ex:
            print(ex)
        self.send("Top 40", data)


"""
Run Container Thread
"""


class RunAnalysisThread(QThread):
    analysis_progress = pyqtSignal(int)

    container_aligns_dir = "/home/user/Counts"
    container_results_dir = "/home/user/Results"
    container_config_dir = "/home/user/Configs"
    container_param_dir = "/home/user/Params"

    def __init__(self, cli, image_name, host_aligns_dir, host_results_dir, host_config_dir, host_param_dir):
        QThread.__init__(self)
        self.docker = cli
        self.image_name = image_name
        self.host_aligns_dir = host_aligns_dir
        self.host_results_dir = host_results_dir
        self.host_config_dir = host_config_dir
        self.host_param_dir = host_param_dir

    def __del__(self):
        self.wait()

    """
    Run should first create a container and then start it
    """
    def run(self):

        configs = [os.path.join(self.container_config_dir, x) for x in
                   fnmatch.filter(os.listdir(self.host_config_dir), 'Configs.*')]
        if len(configs) > 0:
            config_file = configs[0]
        else:
            print('No config file found, exit.')
            return

        cmd_s1 = "Rscript /home/user/Programs/Compare-Molecule-Expression.R {} /home/user /home/user/Programs >& /home/user/Results/log.txt".format(config_file)

        # See Steps 6 and 7 of DTOXS SOP Identification of Differentially Expressed Genes
        commands = ["Rscript /home/user/Programs/Extract-Gene-Expression-Samples.R /home/user/Counts",
                    cmd_s1,
                    "exit"]

        print(commands)
        volumes = {self.host_aligns_dir: self.container_aligns_dir,
                   self.host_results_dir: self.container_results_dir,
                   self.host_config_dir: self.container_config_dir,
                   self.host_param_dir: self.container_param_dir}

        response = self.docker.create_container(self.image_name,
                                                volumes=volumes,
                                                commands=commands)
        if response['Warnings'] is None:
            self.containerId = response['Id']
            self.docker.start_container(self.containerId)
        else:
            print(response['Warnings'])

        # Keep running until container is exited
        while self.docker.container_running(self.containerId):
            self.sleep(1)
        # Remove the container now that it is finished
        self.docker.remove_container(self.containerId)


if __name__ == "__main__":
    import sys
    from AnyQt.QtWidgets import QApplication

    a = QApplication(sys.argv)
    ow = OWDtoxsAnalysis()
    ow.show()
    a.exec_()
    ow.saveSettings()
