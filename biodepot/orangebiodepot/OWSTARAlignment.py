import sys, os, fnmatch
import re
from Orange.widgets import widget, gui
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread
from PyQt5.QtCore import QThread, pyqtSignal

class OWSTARAlignment(widget.OWWidget):
    RunMode_GenerateGenome = 0
    RunMode_STAR = 1

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
        self.result_folder_name = "/Aligned"

        # GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.info_fastq = gui.widgetLabel(box, 'Please specify a directory with fastq or fasta files.')
        self.info_starindex = gui.widgetLabel(box, 'Please specify a genome directory.')
        self.info_fastq.setWordWrap(True)
        self.info_starindex.setWordWrap(True)

        self.RunMode = OWSTARAlignment.RunMode_GenerateGenome
        self.cboRunMode = gui.comboBox(self.controlArea, self, 'RunMode', label='Running Mode:', items=('Generate Genome', 'STAR Alignment'), valueType=int)
        self.cboRunMode.setCurrentIndex(0)

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
                self.info_fastq.setText("Fasta/q files: {0}".format(self.fastqDirectory))

            host_fastq_dir = path.strip()
            # Jimmy May-03-2017, once the fastq input was set, automatically create an output fold as a sibling of "fastq"
            parent_path = os.path.abspath(os.path.join(host_fastq_dir, '..'))
            self.host_results_dir = os.path.join(parent_path + self.result_folder_name)
            if not os.path.exists(self.host_results_dir):
                os.makedirs(self.host_results_dir)

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
                                                     self.host_results_dir,
                                                     self.RunMode)

        self.run_staralignment_thread.analysis_progress.connect(self.run_staralignment_progress)
        self.run_staralignment_thread.finished.connect(self.run_staralignment_done)
        self.run_staralignment_thread.analysis_message.connect(self.run_star_message)
        self.run_staralignment_thread.start()

    def run_staralignment_progress(self, val):
        self.progressBarSet(val)

    def run_staralignment_done(self):
        self.info_fastq.setText("Finished running analysis!")
        self.btnStarAlignment.setEnabled(True)
        self.btnStarAlignment.setText('Run again')
        self.setStatusMessage('Finished!')
        output_channel = self.host_results_dir
        if self.RunMode == OWSTARAlignment.RunMode_GenerateGenome:
            output_channel = self.starindexDirectory
        self.send("Directory", output_channel)
        self.btnStarAlignment.setEnabled(True)

    def run_star_message(self, msgType, message):
        if msgType == 'warning':
            self.warning(message)
        elif msgType == 'error':
            self.error(message)


class StarAlignmentThread(QThread):
    analysis_progress = pyqtSignal(int)
    analysis_message = pyqtSignal(str, str)     # warning/error    message

    container_fastaq_dir = '/data/fastaq'
    container_genome_dir = '/data/genome'
    container_aligned_dir = '/data/aligned'

    # biodepot/ubuntu-star
    def __init__(self, cli, image_name, image_version, host_fastq_dir, host_starindex_dir, host_results_dir, running_mode):
        QThread.__init__(self)
        self.docker = cli
        self.image_name = image_name
        self.image_version = image_version
        self.host_fastq_dir = host_fastq_dir
        self.host_starindex_dir = host_starindex_dir
        self.host_results_dir = host_results_dir
        self.running_mode = running_mode

        self.parameters = ['STAR',
                      '--runThreadN', '8',
                      '--genomeDir', StarAlignmentThread.container_genome_dir]
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

    def __del__(self):
        self.wait()

    """
    Run should first create a container and then start it
    """
    def run(self):

        commands = '';

        if self.running_mode == OWSTARAlignment.RunMode_GenerateGenome:
            self.parameters.extend(['--runMode', 'genomeGenerate'])
            # search fasta files
            fastafiles = [os.path.join(self.container_fastaq_dir, x) for x in fnmatch.filter(os.listdir(self.host_fastq_dir), '*.fa')]
            if not fastafiles:
                self.analysis_message.emit('error', 'No fasta files found')
                return
            self.parameters.extend(['--genomeFastaFiles'])
            self.parameters.extend(fastafiles)
            self.parameters.extend(['--outFileNamePrefix', StarAlignmentThread.container_genome_dir + '/'])

            commands = ' '.join((str(w) for w in self.parameters))
        else:
            # search fastq files
            fastq_files = [os.path.join(self.container_fastaq_dir, x) for x in os.listdir(self.host_fastq_dir)]

            r1, r2 = [], []
            # Pattern convention: Look for "R1" / "R2" in the filename, or "_1" / "_2" before the extension
            pattern = re.compile('(?:^|[._-])(R[12]|[12]\.f)')
            for fastq in sorted(fastq_files):
                match = pattern.search(os.path.basename(fastq))
                if not match:
                    print('Invalid FASTQ file: {0}, ignore it.'.format(fastq))
                    continue
                elif '1' in match.group():
                    r1.append(fastq)
                elif '2' in match.group():
                    r2.append(fastq)

            if len(r1) != len(r2):
                self.analysis_message.emit('error', 'Check fastq names, uneven number of pairs found.')
                print('Check fastq names, uneven number of pairs found.\nr1: {}\nr2: {}'.format(r1, r2))
                return

            if len(r1) == 0:
                self.analysis_message.emit('warning', 'No fastq files found.')
                return

            # generate batch commands
            batch_commands=[]
            for i in range(len(r1)):
                readin_files = ' --readFilesIn {0} {1}'.format(r1[i], r2[i])
                readin_command = ''
                _, file_extension = os.path.splitext(r1[i])
                if file_extension == '.gz':
                    readin_command = ' --readFilesCommand zcat'

                basename = os.path.basename(r1[i]).split('_')[0]
                output_prefix = ' --outFileNamePrefix {0}/{1}.'.format(StarAlignmentThread.container_aligned_dir, basename)
                batch_commands.append(' '.join((str(w) for w in self.parameters)) + readin_files + output_prefix + readin_command)

            commands = 'bash -c "' + ';'.join(c for c in batch_commands) + '"'

        print (commands)

        volumes = {self.host_fastq_dir: self.container_fastaq_dir,
                   self.host_starindex_dir: self.container_genome_dir,
                   self.host_results_dir: self.container_aligned_dir}

        print (volumes)
        response = self.docker.create_container(self.image_name+":"+self.image_version,
                                                volumes=volumes,
                                                commands=commands)
        if response['Warnings'] == None:
            self.containerId = response['Id']
            self.docker.start_container(self.containerId)
        else:
            print(response['Warnings'])

        # Keep running until container is exited
        while self.docker.container_running(self.containerId):
            self.sleep(1)
        # Remove the container now that it is finished
        self.docker.remove_container(self.containerId)

def main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv))

    ow = OWSTARAlignment()
    ow.setFastqInput("/Users/Jimmy/Downloads/STAR/fasta_input")
    ow.setGenomeDir("/Users/Jimmy/Downloads/STAR/genome")
    ow.show()
    ow.raise_()

    ow.handleNewSignals()
    app.exec_()
    ow.setFastqInput(None)
    ow.handleNewSignals()
    return 0

if __name__=="__main__":
    sys.exit(main())