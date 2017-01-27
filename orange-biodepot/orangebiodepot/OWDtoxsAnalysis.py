import os
import Orange.data
from Orange.widgets import widget, gui
from PyQt4.QtCore import QThread, SIGNAL
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread

class OWDtoxsAnalysis(widget.OWWidget):
    name = "Dtoxs Analysis"
    description = "Step 2 of Dtoxs SOP. Uses edgeR for differential expression analysis."
    category = "RNASeq"
    icon = "icons/dtoxs-analysis2.svg"
    priority = 10

    image_name = "biodepot/dtoxs_analysis"
    image_version = "latest"

    inputs = [("Counts", str, "set_aligns")]
    outputs = [("Results", str),("Top 40", Orange.data.Table)]

    want_main_area = False

    # The directory of the alignment data needed to run
    # the container will be set by the user before the
    # container can be run
    host_counts_dir = ""

    # The default write location of all widgets should be
    # ~/BioDepot/WidgetName/
    # TODO is this an issue multiple containers write to the same place?
    host_results_dir = os.path.expanduser('~') + '/BioDepot/Dtoxs_Analysis/Results'

    def __init__(self):
        super().__init__()

        if not os.path.exists(self.host_results_dir):
            os.makedirs(self.host_results_dir)
        # Docker Client
        # This client talks to your local docker
        self.docker = DockerClient('unix:///var/run/docker.sock', 'local')

        # GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.infoLabel = gui.widgetLabel(box, 'Connect to Dtoxs Alignment or use a Directory widget'
                                              ' to specify location of UMI read counts data files')
        self.infoLabel.setWordWrap(True)

        # TODO let user specify results directory
        # self.btn_res = gui.button(None, self, "Open", callback=self.set_results_dir)
        # self.resultsEdit = QtGui.QLineEdit(self.host_results_dir)
        # self.buttonsArea.layout().addWidget(self.btn_res)
        # self.buttonsArea.layout().addWidget(self.resultsEdit)

        self.autoRun = True
        gui.checkBox(self.controlArea, self, 'autoRun', 'Run automatically when input set')
        self.btn_run = gui.button(self.controlArea, self, "Run", callback=self.start_analysis)

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
            if self.autoRun:
                self.start_analysis()
            else:
                self.infoLabel.setText('Alignment directory set.\nWaiting to run...')

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
        self.connect(self.pull_image_thread, SIGNAL("pull_progress"), self.pull_image_progress)
        self.connect(self.pull_image_thread, SIGNAL("finished()"), self.pull_image_done)
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
        elif self.host_counts_dir:
            self.run_analysis()
        else:
            self.infoLabel.setText('Set alignment directory before running.')

    def run_analysis(self):
        self.infoLabel.setText('Running analysis...')
        self.setStatusMessage('Running...')
        self.progressBarInit()
        # Run the container in a new thread
        self.run_analysis_thread = RunAnalysisThread(self.docker,
                                                     self.image_name,
                                                     self.host_counts_dir,
                                                     self.host_results_dir)

        self.connect(self.run_analysis_thread, SIGNAL('analysis_progress'), self.run_analysis_progress)
        self.connect(self.run_analysis_thread, SIGNAL("finished()"), self.run_analysis_done)
        self.run_analysis_thread.start()

    def run_analysis_progress(self, val):
        self.progressBarSet(val)

    def run_analysis_done(self):
        self.infoLabel.setText("Finished running analysis!")
        self.btn_run.setEnabled(True)
        self.btn_run.setText('Run again')
        self.setStatusMessage('Finished!')
        self.progressBarFinished()
        self.send("Results", self.host_results_dir)
        # TODO create Orange.data.Table from TOP-40.tsv


"""
Run Container Thread
"""
class RunAnalysisThread(QThread):

    container_aligns_dir = "/home/user/Counts"
    container_results_dir = "/home/user/Results"

    # See Steps 6 and 7 of DTOXS SOP Identification of Differentially Expressed Genes
    commands = ["Rscript /home/user/Programs/Extract-Gene-Expression-Samples.R /home/user/Counts",
                "Rscript /home/user/Programs/Compare-Molecule-Expression.R "
                "/home/user/Configs/Configs.LINCS.Dataset.Gene.LINCS.20150409.tsv "
                "/home/user "
                "/home/user/Programs"
                ">& /home/user/Results/log.txt",
                "exit"]

    def __init__(self, cli, image_name, host_aligns_dir, host_results_dir):
        QThread.__init__(self)
        self.docker = cli
        self.image_name = image_name
        self.host_aligns_dir = host_aligns_dir
        self.host_results_dir = host_results_dir

    def __del__(self):
        self.wait()

    """
    Run should first create a container and then start it
    """
    def run(self):
        volumes = { self.host_aligns_dir: self.container_aligns_dir,
                    self.host_results_dir: self.container_results_dir }
        response = self.docker.create_container(self.image_name,
                                                volumes=volumes,
                                                commands=self.commands)
        if response['Warnings'] == None:
            self.containerId = response['Id']
            self.docker.start_container(self.containerId)
            #self.docker.remove_container(self.containerId)
        else:
            # TODO emit warning to Widget
            print(response['Warnings'])
        i = 1
        # Keep running until container is exited
        while self.docker.container_running(self.containerId):
            self.emit(SIGNAL('analysis_progress'), i * 0.55)
            self.sleep(1)
            i += 1
        # Remove the container now that it is finished
        self.docker.remove_container(self.containerId)