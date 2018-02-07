import os, time, fnmatch, re
from PyQt5.QtWidgets import *

from Orange.widgets import widget, gui
from Orange.widgets.settings import Setting
from AnyQt.QtCore import QObject, QThread, pyqtSignal, pyqtSlot

from orangebiodepot.util.DockerClient import DockerClient, PullImageThread

class OWDtoxsAlignment(widget.OWWidget):
    name = "Dtoxs Alignment"
    description = "Step 1 of the Dtoxs SOP. Uses a Burrows-Wheeler Aligner (BWA)."
    category = "RNASeq"
    icon = "icons/dtoxs-alignment2.svg"
    priority = 10

    image_name = "biodepot/dtoxs_alignment"
    image_version = "latest"

    inputs = [("References", str, "set_refs"),
              ("Seqs", str, "set_seqs"),
              ("Configs", str, "set_configs"),
              ("Barcodes", str, "set_barcodes_file"),
              ("Trigger",str,"onTrigger")]
    outputs = [("Counts", str)]
   
    hostDirectories = Setting({'refs': None, 'seqs': None, 'conf': None, 'counts': None, 'aligns': None, 'barcodes': None},schema_only=True)
    auto_run=Setting(False, schema_only=True)
    waitOnTrigger=Setting(False, schema_only=True)
    want_main_area = False

    def __init__(self):
        super().__init__()

        # Docker Client
        # This client talks to your local docker
        self.docker = DockerClient('unix:///var/run/docker.sock', 'local')
        
        self.labelText = {'refs': 'References: {0}', 'seqs': 'Sequences: {0}', 'conf': 'Configs: {0}', 'barcodes': 'Barcode file: {0}'}

        # folders that will created automatically
        self.result_folder_name = "Counts"
        self.aligns_folder_name = "Aligns"
        
        #add this to make sure the the checkboxes don't disable each other
        if self.waitOnTrigger and self.auto_run:
            self.auto_run=False

        # GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.lblRefs = gui.widgetLabel(box, 'References: Set') if self.hostDirectories['refs'] else gui.widgetLabel(box, 'References: Not Set')
        self.lblSeqs = gui.widgetLabel(box, 'Sequences: Set') if self.hostDirectories['seqs'] else gui.widgetLabel(box, 'Sequences: Not Set')
        self.lblConfigs = gui.widgetLabel(box, 'Configs: Set') if self.hostDirectories['conf'] else gui.widgetLabel(box, 'Configs: Not Set')
        self.lblBarcodes = gui.widgetLabel(box, "Barcode file: Set") if self.hostDirectories['barcodes'] else gui.widgetLabel(box, 'Barcode file: Not Set')

        self.infoLabel = gui.widgetLabel(box, 'Waiting...')
        self.infoLabel.setWordWrap(True)
        self.triggerReceived=False

        self.auto_runBox=gui.checkBox(self.controlArea, self, 'auto_run', 'Run automatically when input set',callback=self.setRunType)
        self.waitOnTriggerBox=gui.checkBox(self.controlArea, self, 'waitOnTrigger','Run when trigger signal received',callback=self.setRunType)
        self.btn_run = gui.button(self.controlArea, self, "Run", callback=self.btn_run_pushed)
        self.is_running = False

        self.setMinimumSize(500, 200)

    def setRunType(self):
        if self.waitOnTrigger and self.auto_run:
            self.auto_run=False
    """
    Called when the user pushes 'Run'
    """
    def btn_run_pushed(self):
        self.start_container(run_btn_pushed=True)

    """
    Set references
    """
    def set_refs(self, path):
        self.__set_directory__('refs', self.lblRefs, path)

    """
    Set seqs
    """
    def onTrigger(self, value, sourceId=None):
        if value and self.waitOnTrigger:
            self.triggerReceived=True
            self.start_container()    
    
    def set_seqs(self, path):
        self.__set_directory__('seqs', self.lblSeqs, path, startContainer=False)

        if self.hostDirectories['seqs'] is not None:
            # Jimmy March-28-2017, once the seq input was set, automatically create a result and aligns folder as a sibling of "Seqs"
            parent_path = os.path.abspath(os.path.join(self.hostDirectories['seqs'], '..'))
            self.hostDirectories['counts'] = os.path.join(parent_path, self.result_folder_name)
            self.hostDirectories['aligns'] = os.path.join(parent_path, self.aligns_folder_name)

            if not os.path.exists(self.hostDirectories['counts']):
                os.makedirs(self.hostDirectories['counts'])
            if not os.path.exists(self.hostDirectories['aligns']):
                os.makedirs(self.hostDirectories['aligns'])

        self.start_container(run_btn_pushed=False)

    def set_configs(self, path):
        self.__set_directory__('conf', self.lblConfigs, path)

    def __set_directory__(self, key, ctrlLabel, path, startContainer = True):
        # When a user removes a connected Directory widget,
        # it sends a signal with path=None
        # do not check for non-existence if waiting for signal (download)
        if path is None:
            self.hostDirectories[key] = None
            ctrlLabel.setText(self.labelText[key].format('Not Set'))
        elif not os.path.exists(path) and not self.waitOnTrigger:
            self.hostDirectories[key] = None
            ctrlLabel.setText(self.labelText[key].format('None existent directory: {}'.format(str(path))))
        else:
            self.hostDirectories[key] = path.strip()
            if key == 'barcodes':
                path = os.path.basename(path)
                self.hostDirectories['barcodes_basename'] = path
            ctrlLabel.setText(self.labelText[key].format((str(path))))

        if startContainer:
            self.start_container(run_btn_pushed=False)

    def set_barcodes_file(self, filename):
        self.__set_directory__('barcodes', self.lblBarcodes, filename)
    """
    Pull image
    """
    def pull_image(self):
        self.infoLabel.setText('Pulling \'' + self.image_name + ":" + self.image_version + '\' from Dockerhub...')
        self.setStatusMessage("Downloading...")
        self.progressBarInit()
        self.is_running = True
        self.btn_run.setEnabled(False)
        # Pull the image in a new thread
        self.pull_image_worker = PullImageThread(self.docker, self.image_name, self.image_version)
        self.pull_image_worker.pull_progress.connect(self.pull_image_progress)
        self.pull_image_worker.finished.connect(self.pull_image_finished)
        self.pull_image_worker.start()

    def pull_image_progress(self, val):
        self.progressBarSet(val)

    def pull_image_finished(self):
        self.infoLabel.setText('Finished pulling \'' + self.image_name + ":" + self.image_version + '\'')
        self.progressBarFinished()
        self.run_container()

    """
    Alignment
    """
    def start_container(self, run_btn_pushed=True):
        # Make sure both inputs are set
        all_set = all(value is not None for value in self.hostDirectories.values())

        if all_set and not self.is_running:
            if run_btn_pushed or (self.auto_run and not self.waitOnTrigger) or (self.waitOnTrigger and self.triggerReceived):
                # Make sure the docker image is downloaded
                if not self.docker.has_image(self.image_name, self.image_version):
                    self.pull_image()
                else:
                    self.run_container()


    def run_container(self):
        self.is_running = True
        self.btn_run.setEnabled(False)
        self.infoLabel.setText('Running alignment...')
        self.setStatusMessage('Running...')
        #self.progressBarInit()
        # Run the container in a new thread
        self.run_container_thread = RunAlignmentThread(self.docker,
                                                       self.image_name,
                                                       self.hostDirectories)
        self.run_container_thread.progress.connect(self.run_container_progress)
        self.run_container_thread.finished.connect(self.run_container_finished)
        self.run_container_thread.start()

    def run_container_progress(self, val):
        self.progressBarSet(val)

    def run_container_finished(self):
        self.infoLabel.setText("Finished running alignment!")
        self.btn_run.setEnabled(True)
        self.is_running = False
        self.btn_run.setText('Run again')
        self.setStatusMessage('Finished!')
        #self.progressBarFinished()
        self.send("Counts", self.hostDirectories['counts'])


"""
Run Alignment Worker
"""
class RunAlignmentThread(QThread):
    progress = pyqtSignal(int)

    container_top_dir = "/root/LINCS"
    container_ref_dir = "/root/LINCS/References"
    container_seq_dir = "/root/LINCS/Seqs"
    container_counts_dir = "/root/LINCS/Counts"
    container_aligns_dir = "/root/LINCS/Aligns"
    container_config_dir = "/root/LINCS/Configs"

    commands = ["/root/LINCS/Programs/Broad-DGE/run-alignment-analysis.sh >& /root/LINCS/Counts/run-alignment-analysis.log; "
                "exit"]

    def __init__(self, cli, image_name, hostDirectories):
        QThread.__init__(self)

        self.docker = cli
        self.image_name = image_name
        self.hostDirectories = hostDirectories
        self.containerId = ""

    def __del__(self):
        self.wait()

    def run(self):
        target_barcode_file = os.path.join(self.container_top_dir, self.hostDirectories['barcodes_basename'])

        # scan config to fetch series name
        configs = [x for x in fnmatch.filter(os.listdir(self.hostDirectories['conf']), 'Configs.*')]
        if len(configs) > 0:
            config_file = configs[0]
        else:
            print('No config file found, exit')
            return

        series_name = ''
        matched = re.compile(r"[0-9]{8}").search(config_file)
        if matched:
            series_name = matched.group(0)
        else:
            print('No series name matched, exit')
            return

        # May-31-2017 Jimmy, added for scanning lanes from all fastq files
        lanes = 6
        try:
            fastq_files = [x for x in fnmatch.filter(os.listdir(self.hostDirectories['seqs']), '*.fastq*')]
            all_lanes = re.compile(r"Lane[0-9]").findall(','.join(fastq_files))
            all_lanes = sorted(list(set(all_lanes)))
            all_lanes = re.compile(r"\d+").findall(','.join(all_lanes))
            all_lanes = list(map(int, all_lanes))
            lanes = max(all_lanes)
        except:
            print('Error on match LANES, exit')
            return

        print('LANES=', lanes)

        environment = {'ENV_SERIES_NAME': series_name,
                       'ENV_BARCODE_FILE': self.hostDirectories['barcodes_basename'],
                       'ENV_LANES': lanes}
        #print(environment)

        volumes = { self.hostDirectories['refs']: self.container_ref_dir,
                    self.hostDirectories['seqs']: self.container_seq_dir,
                    self.hostDirectories['conf']: self.container_config_dir,
                    #self.hostDirectories['param']: self.container_param_dir,
                    self.hostDirectories['counts']: self.container_counts_dir,
                    self.hostDirectories['aligns']: self.container_aligns_dir,
                    self.hostDirectories['barcodes']: target_barcode_file}
        #print(volumes)

        response = self.docker.create_container(self.image_name,
                                                volumes=volumes,
                                                environment=environment,
                                                commands=self.commands)
        if response['Warnings'] == None:
            self.containerId = response['Id']
            self.docker.start_container(self.containerId)
        else:
            print(response['Warnings'])

        #i = 1
        # Keep running until container is exited
        while self.docker.container_running(self.containerId):
            time.sleep(1)
            pass
        # Remove the container now that it is finished
        self.docker.remove_container(self.containerId)


if __name__ == "__main__":
    import sys
    from AnyQt.QtWidgets import QApplication

    a = QApplication(sys.argv)
    ow = OWDtoxsAlignment()
    ow.set_seqs('/Users/Jimmy/Developer/bioinfomatics/sample_data/dtoxs/seqs')
    ow.set_refs('/Users/Jimmy/Developer/bioinfomatics/sample_data/dtoxs/refs')
    ow.set_barcodes_file('/Users/Jimmy/Developer/bioinfomatics/sample_data/dtoxs/barcodes_trugrade_96_set2.dat')
    ow.set_configs('/Users/Jimmy/Developer/bioinfomatics/sample_data/dtoxs/configs')
    ow.show()
    a.exec_()
    ow.saveSettings()
