import os
import sys
from Orange.widgets import widget, gui
import pysam

class OWBamSort(widget.OWWidget):
    name = "BAM sorting"
    description = "Step 2 of sorting a BAM file. Used the bam file widget to read a bam file"
    category = "RNASeq"
    icon = "icons/BamRank.svg"
    priority = 10

    inputs = [("File", str, "set_aligns")]
    outputs = [("Results", str)]

    want_main_area = False

    # TODO is this an issue if multiple containers write to the same place?
    host_results_dir = "" #os.path.expanduser('~') + '/BioDepot/Dtoxs_Analysis/Results'

    def __init__(self):
        super().__init__()
		
        #if not os.path.exists(self.host_results_dir):
        #    os.makedirs(self.host_results_dir)

        #GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.infoLabel = gui.widgetLabel(box, 'Connect to BAM file widget '
                                              ' to specify location of the BAM file to sort')
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
            print('Tried to set alignment directory to None')
        elif not os.path.exists(path):
            # TODO create warning
            print('Tried to set alignment directory to none existent directory: ' + str(path))
        else:
            self.host_counts_dir = path.strip()
            self.host_results_dir = os.path.dirname(self.host_counts_dir)
            if self.autoRun:
                self.start_analysis()
            else:
                self.infoLabel.setText('Alignment directory set.\nWaiting to run...')

    """
    Anaysis
    """
    def start_analysis(self):
        if self.host_counts_dir:
            self.run_analysis()
        else:
            self.infoLabel.setText('Set alignment directory before running.')

    def run_analysis(self):
        self.infoLabel.setText('Running analysis...')
        self.setStatusMessage('Running...')
		
        pysam.AlignmentFile(self.host_counts_dir, "rb")
        pysam.sort("-o", self.host_results_dir + "/output.sorted.bam" , self.host_counts_dir)
		
        self.infoLabel.setText('Sorted bam file. The output is stored at ' + self.host_results_dir + "/output.sorted.bam")

def main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv))

    ow = OWBamSort()
    ow.set_aligns("/Users/Jimmy/Downloads/Bam files/toy.bam")
    ow.show()
    ow.raise_()

    ow.handleNewSignals()
    app.exec_()
    ow.handleNewSignals()
    return 0

if __name__=="__main__":
    sys.exit(main())