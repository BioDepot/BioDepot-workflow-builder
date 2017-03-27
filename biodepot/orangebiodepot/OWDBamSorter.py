import os
import Orange.data
from Orange.widgets import widget, gui
from PyQt5.QtCore import QThread, pyqtSignal
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread
import pysam

class OWDBamSorter(widget.OWWidget):
    name = "BAM sorter"
    description = "Step 2 of soritng bam file. Uses BAM file widget to read bam file."
    category = "RNASeq"
    icon = "icons/BamRank.svg"
    priority = 10

    inputs = [("Counts", str, "set_aligns")]
    outputs = [("Results", str)]

    want_main_area = False

    # The directory of the alignment data needed to run
    # the container will be set by the user before the
    # container can be run
    host_counts_dir = ""

    # The default write location of all widgets should be
    # ~/BioDepot/WidgetName/
    # TODO is this an issue if multiple containers write to the same place?
    host_results_dir = "" #os.path.expanduser('~') + '/BioDepot/Bam_Sort/Results'

    def __init__(self):
        super().__init__()

        #if not os.path.exists(self.host_results_dir):
        #    os.makedirs(self.host_results_dir)

        # GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.infoLabel = gui.widgetLabel(box, 'Connect to Bam file  widget'
                                              ' to specify location of bam file to sort.')
        self.infoLabel.setWordWrap(True)

        self.autoRun = True
        gui.checkBox(self.controlArea, self, 'autoRun', 'Run automatically when input set')
        self.btn_run = gui.button(self.controlArea, self, "Run", callback=self.start_analysis)

    """
    Set input
    """
    def set_aligns(self, path):
        if not type(path) is str:
            # TODO create warning
            print('Tried to set the bam file to None')
        elif not os.path.exists(path):
            # TODO create warning
            print('Tried to set bam file to non existant path: ' + str(path))
        else:
            self.host_counts_dir = path.strip()
            self.host_results_dir = os.path.dirname(self.host_counts_dir)
            if self.autoRun:
                self.start_analysis()
            else:
                self.infoLabel.setText('Bam file set.\nWaiting to run...')

    """
    Analysis
    """
    def start_analysis(self):
        # Make sure the docker image is downloaded
        #if not self.docker.has_image(self.image_name, self.image_version):
            #self.pull_image()
        # Make sure there is an alignment directory set
        if self.host_counts_dir:
            self.run_analysis()
        else:
            self.infoLabel.setText('Set bam file before running.')

    def run_analysis(self):
        self.infoLabel.setText('Running analysis...')
        self.setStatusMessage('Running...')
        pysam.sort("-o", str(self.host_results_dir) + "/output.sorted.bam", str(self.host_counts_dir))
        self.infoLabel.setText('Sorted. The output file is at ' + str(self.host_results_dir) + "/output.sorted.bam")
        self.setStatusMessage('Sorting completed.')
