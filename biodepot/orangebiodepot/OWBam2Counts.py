import os
from orangebiodepot.util.BwBase import OWBwBWidget
from Orange.widgets import gui

class OWSTARAlignment(OWBwBWidget):
    name = "BAM to Counts"
    description = "Counts BAM files"
    category = "General"
    icon = "icons/bam2counts.svg"

    priority = 10

    inputs = [("GTF file", str, "setGtfFile"),
              ("Sample Table", str, "setSampleTable"),
              ("BAM files", str, "setBamFiles")]

    outputs = [("Counts", str)]

    want_main_area = False

    docker_image_name = 'biodepot/bam2counts'
    docker_image_tag = 'latest'

    # worker numbers
    worker_numbers = 4

    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)

        # GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.lblGtf = gui.widgetLabel(box, 'GTF File: Not Set')
        self.lblSampleTable = gui.widgetLabel(box, 'Sample table: Not Set')
        self.lblBamfiles = gui.widgetLabel(box, 'BAM files: Not Set')

        self.infoLabel = gui.widgetLabel(box, 'Waiting...')
        self.infoLabel.setWordWrap(True)

        self.auto_run = False
        gui.checkBox(self.controlArea, self, 'auto_run', 'Run automatically when input set')
        self.btnRun = gui.button(self.controlArea, self, "Run", callback=self.OnRunClicked)
        self.Flag_isRunning = False

        self.Initialize_InputKeys(['gtf', 'sampletable', 'bamfiles'])

        self.setMinimumSize(500, 200)

    def OnRunClicked(self):
        self.startJob(triggerdByButton=True)

    def setGtfFile(self, path):
        self.setDirectories('gtf', path, self.lblGtf)
        self.startJob()

    def setSampleTable(self, path):
        self.setDirectories('sampletable', path, self.lblSampleTable)
        self.startJob()

    def setBamFiles(self, path):
        self.setDirectories('bamfiles', path, self.lblBamfiles)
        self.startJob()

    def startJob(self, triggerdByButton = False):
        all_set = all(value is not None for value in [self.getDirectory('gtf'), self.getDirectory('sampletable'), self.getDirectory('bamfiles')])

        if all_set and (self.auto_run or triggerdByButton):
            self.btnRun.setEnabled(False)
            top_dir = '/data'
            gtf = os.path.join(top_dir, os.path.basename(self.getDirectory('gtf')))
            sample_table = os.path.join(top_dir, os.path.basename(self.getDirectory('sampletable')))
            bamfiles = os.path.join(top_dir, 'bamfiles')
            logfile = os.path.join(bamfiles, 'log.txt')

            volumes = {self.getDirectory('gtf'): gtf,
                       self.getDirectory('sampletable'): sample_table,
                       self.getDirectory('bamfiles'): bamfiles}

            cmd = 'Rscript /root/bam2counts.R {} {} {} {} >& {}'.format(gtf, sample_table, bamfiles, self.worker_numbers, logfile)
            commands = [cmd,"exit"]

            self.dockerRun(volumes, commands)


    def Event_OnRunFinished(self):
        self.btnRun.setEnabled(True)
        self.send('Counts', self.getDirectory('bamfiles'))

    def Event_OnRunMessage(self, message):
        self.infoLabel.setText(message)
