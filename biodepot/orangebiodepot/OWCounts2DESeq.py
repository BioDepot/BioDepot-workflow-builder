import os
import Orange.data
from Orange.data.io import FileFormat
from orangebiodepot.util.BwBase import OWBwBWidget
from Orange.widgets import gui


class OWCounts2Deseq(OWBwBWidget):
    name = "Counts to DESeq"
    description = "DESeq counts table"
    category = "General"
    icon = "icons/counts2deseq.svg"

    priority = 10

    inputs = [("Counts table", str, "setCountsTable"),
              ("Sample Table", str, "setSampleTable")]

    outputs = [("Results", str), ("DataTable", Orange.data.Table)]

    want_main_area = False

    docker_image_name = 'biodepot/bam2deseq'
    docker_image_tag = 'latest'

    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        # GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.lblCounts = gui.widgetLabel(box, 'Counts Table: Not Set')
        self.lblSampleTable = gui.widgetLabel(box, 'Sample table: Not Set')

        self.infoLabel = gui.widgetLabel(box, 'Waiting...')
        self.infoLabel.setWordWrap(True)

        self.auto_run = False
        gui.checkBox(self.controlArea, self, 'auto_run', 'Run automatically when input set')
        self.btnRun = gui.button(self.controlArea, self, "Run", callback=self.OnRunClicked)
        self.Flag_isRunning = False

        self.Initialize_InputKeys(['countstable', 'sampletable'])

        self.setMinimumSize(500, 180)

    def OnRunClicked(self):
        self.startJob(triggerdByButton=True)

    def setCountsTable(self, path):
        self.setDirectories('countstable', path, self.lblCounts)
        self.startJob()

    def setSampleTable(self, path):
        self.setDirectories('sampletable', path, self.lblSampleTable)
        self.startJob()

    def startJob(self, triggerdByButton=False):
        all_set = all(value is not None for value in [self.getDirectory('countstable'), self.getDirectory('sampletable')])

        if all_set and (self.auto_run or triggerdByButton):
            self.btnRun.setEnabled(False)
            top_dir = '/data'
            counts_table = os.path.join(top_dir, os.path.basename(self.getDirectory('countstable')))
            sample_table = os.path.join(top_dir, os.path.basename(self.getDirectory('sampletable')))
            self.output_dir_host = os.path.dirname(os.path.abspath(self.getDirectory('countstable')))
            output_dir = os.path.join(top_dir, 'output')
            logfile = os.path.join(output_dir, 'log.txt')

            volumes = {self.getDirectory('countstable'): counts_table,
                       self.getDirectory('sampletable'): sample_table,
                       self.output_dir_host: output_dir}

            cmd = 'Rscript /home/root/counts2DESeq.R {} {} {} >& {}'.format(counts_table, sample_table, output_dir, logfile)
            commands = [cmd, "exit"]

            self.dockerRun(volumes, commands)

    def Event_OnRunFinished(self):
        self.btnRun.setEnabled(True)
        self.send('Results', self.output_dir_host)
        tsvFile = os.path.join(self.output_dir_host, 'deseq_results.csv');
        tsvReader = FileFormat.get_reader(tsvFile)
        data = None
        try:
            data = tsvReader.read()
        except Exception as ex:
            print(ex)
        self.send("DataTable", data)

    def Event_OnRunMessage(self, message):
        self.infoLabel.setText(message)
