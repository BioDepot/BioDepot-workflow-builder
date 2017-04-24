import sys
import os
import fnmatch
from Orange.widgets import widget, gui
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread
from PyQt5.QtCore import QThread, pyqtSignal

class OWSTARAlignment(widget.OWWidget):
    name = "Star Alignment"
    description = "Performs alignment of fastqs to bam via STAR"
    category = "RNASeq"
    icon = "icons/staralignment.svg"

    priority = 10

    inputs = [("FastQ Files", str, "setFastqInput", widget.Default),
              ("Genome Dir", str, "setGenomeDir")]
    outputs = [("Directory", str)]

    want_main_area = False

    dockerClient = DockerClient('unix:///var/run/docker.sock', 'local')

    image_name = "biodepot/ubuntu-star"
    image_version = "latest"

    def __init__(self):
        super().__init__()

        self.bFastqDirSet = False
        self.bStarIndexDirSet = False

        # GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.info_fastq = gui.widgetLabel(box, 'Please specify a directory with fastq files.')
        self.info_starindex = gui.widgetLabel(box, 'Please specify a genome directory.')
        self.info_fastq.setWordWrap(True)
        self.info_starindex.setWordWrap(True)

        self.AutoStarAlignment = False
        gui.checkBox(self.controlArea, self, 'AutoStarAlignment', 'Run automatically when directory was set.')

        self.btnStarAlignment = gui.button(self.controlArea, self, "Star Alignment", callback=self.StartStarAlignment)
        self.btnStarAlignment.setEnabled(False)

    def setFastqInput(self, path):
        # When a user removes a connected Directory widget,
        # it sends a signal with path=None
        if path is None:
            self.bFastqDirSet = False
        else:
            if not os.path.exists(path):
                self.bFastqDirSet = False
                print('References set to invalid directory')
            else:
                self.fastqDirectory = path
                self.bFastqDirSet = True
                self.info_fastq.setText("Fastq files: {0}".format(self.fastqDirectory))
        self.btnStarAlignment.setEnabled(self.bFastqDirSet and self.bStarIndexDirSet)

        if self.bFastqDirSet and self.bStarIndexDirSet and self.AutoStarAlignment:
            self.StartStarAlignment()

    def setGenomeDir(self, path):
        if path is None:
            self.bStarIndexDirSet = False
        else:
            if not os.path.exists(path):
                self.bStarIndexDirSet = False
                print('References set to invalid directory')
            else:
                self.starindexDirectory = path
                self.bStarIndexDirSet = True
                self.info_starindex.setText("Star index: {0}".format(self.starindexDirectory))
        self.btnStarAlignment.setEnabled(self.bFastqDirSet and self.bStarIndexDirSet)

        if self.bFastqDirSet and self.bStarIndexDirSet and self.AutoStarAlignment:
            self.StartStarAlignment()

    def StartStarAlignment(self):
        if not self.bFastqDirSet or not self.bStarIndexDirSet:
            return

        if not self.dockerClient.has_image(self.image_name, self.image_version):
            self.pull_image()
        elif self.bFastqDirSet and self.bStarIndexDirSet:
            self.run_staralignment()


    """
        Pull image
    """
    def pull_image(self):
        self.info_fastq.setText('Pulling \'' + self.image_name + ":" + self.image_version)
        self.info_starindex.setText('')
        self.setStatusMessage("Downloading...")
        self.progressBarInit()
        self.btnStarAlignment.setEnabled(False)
        # Pull the image in a new thread
        self.pull_image_thread = PullImageThread(self.dockerClient, self.image_name, self.image_version)
        self.pull_image_thread.pull_progress.connect(self.pull_image_progress)
        self.pull_image_thread.finished.connect(self.pull_image_done)
        self.pull_image_thread.start()

    def pull_image_progress(self, val=0):
        self.progressBarSet(val)

    def pull_image_done(self):
        self.info_fastq.setText('Finished pulling \'' + self.image_name + ":" + self.image_version + '\'')
        self.progressBarFinished()
        self.run_staralignment()

    def run_staralignment(self):
        self.btnStarAlignment.setEnabled(False)
        self.info_fastq.setText('Running star alignment...')
        self.info_starindex.setText('')
        self.setStatusMessage('Running...')
        #self.progressBarInit()
        # Run the container in a new thread
        self.run_staralignment_thread = StarAlignmentThread(self.dockerClient,
                                                     self.image_name,
                                                     self.image_version,
                                                     self.fastqDirectory,
                                                     self.starindexDirectory,
                                                     self.fastqDirectory)

        self.run_staralignment_thread.analysis_progress.connect(self.run_staralignment_progress)
        self.run_staralignment_thread.finished.connect(self.run_staralignment_done)
        self.run_staralignment_thread.start()

    def run_staralignment_progress(self, val):
        self.progressBarSet(val)

    def run_staralignment_done(self):
        self.info_fastq.setText("Finished running analysis!")
        self.btnStarAlignment.setEnabled(True)
        self.btnStarAlignment.setText('Run again')
        self.setStatusMessage('Finished!')
        #self.progressBarFinished()
        self.send("Directory", self.fastqDirectory)
        self.btnStarAlignment.setEnabled(True)


class StarAlignmentThread(QThread):
    analysis_progress = pyqtSignal(int)

    container_fastq_dir = '/data'
    container_genome_dir = '/genome'

    # docker pull biodepot/ubuntu-star
    # quay.io/ucsc_cgl/star:2.4.2a--bcbd5122b69ff6ac4ef61958e47bde94001cfe80

    parameters = ['STAR',
                  '--runThreadN', '8',
                  '--genomeDir', container_genome_dir]
                  # '--outFileNamePrefix', 'rna',
                  # '--outSAMtype', 'BAM', 'SortedByCoordinate',
                  # '--outSAMunmapped', 'Within',
                  # '--quantMode', 'TranscriptomeSAM',
                  # '--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
                  # '--outFilterType', 'BySJout',
                  # '--outFilterMultimapNmax', '20',
                  # '--outFilterMismatchNmax', '999',
                  # '--outFilterMismatchNoverReadLmax', '0.04',
                  # '--alignIntronMin', '20',
                  # '--alignIntronMax', '1000000',
                  # '--alignMatesGapMax', '1000000',
                  # '--alignSJoverhangMin', '8',
                  # '--alignSJDBoverhangMin', '1',
                  # '--sjdbScore', '1',
                  # '--limitBAMsortRAM', '49268954168']



    def __init__(self, cli, image_name, image_version, host_fastq_dir, host_starindex_dir, host_results_dir):
        QThread.__init__(self)
        self.docker = cli
        self.image_name = image_name
        self.image_version = image_version
        self.host_fastq_dir = host_fastq_dir
        self.host_starindex_dir = host_starindex_dir
        self.host_results_dir = host_results_dir

    def __del__(self):
        self.wait()

    """
    Run should first create a container and then start it
    """

    def run(self):
        self.parameters.extend(['--readFilesIn'])

        self.parameters.extend([os.path.join(self.container_fastq_dir, x) for x in fnmatch.filter(os.listdir(self.host_fastq_dir), '*.fastq')])
        commands = ' '.join((str(w) for w in self.parameters))
        print(commands)

        volumes = {self.host_fastq_dir: self.container_fastq_dir,
                   self.host_starindex_dir: self.container_genome_dir}

        print (volumes)
        return
        response = self.docker.create_container(self.image_name+":"+self.image_version,
                                                volumes=volumes,
                                                commands=commands)
        if response['Warnings'] == None:
            self.containerId = response['Id']
            self.docker.start_container(self.containerId)
        else:
            print(response['Warnings'])
        i = 1
        # Keep running until container is exited
        while self.docker.container_running(self.containerId):
            #self.analysis_progress.emit(i)
            self.sleep(1)
            #i += 1
        # Remove the container now that it is finished
        self.docker.remove_container(self.containerId)

def main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv))

    ow = OWSTARAlignment()
    ow.setFastqInput("/Users/Jimmy/Downloads/Seqs/")
    ow.setGenomeDir("/Users/Jimmy/Downloads/Seqs/")
    ow.show()
    ow.raise_()

    ow.handleNewSignals()
    app.exec_()
    ow.setFastQInput(None)
    ow.handleNewSignals()
    return 0

if __name__=="__main__":
    sys.exit(main())