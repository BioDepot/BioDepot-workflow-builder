import sys
import numpy

import Orange.data
from Orange.widgets import widget, gui
from orangebiodepot.util.DockerClient import DockerClient
from PyQt4 import QtGui, QtCore
import datetime

class OWDocker(widget.OWWidget):
    name = "Docker Info"
    description = "See the status of your docker engine"
    category = "Dev"
    icon = "icons/docker-info.svg"
    priority = 1

    inputs = []
    outputs = []

    want_main_area = False
    # TODO finish supporting multiple Docker Clients for Containers, Images, Volumns tables
    # see Version table for an example
    docker_clients = [DockerClient('unix:///var/run/docker.sock', 'local')]

    font = QtGui.QFont()
    font.setFamily("Courier New")
    font.setPointSize(12)

    def __init__(self):
        super().__init__()
        self.context = 'none'

        # GUI
        self.btn_cli = gui.button(None, self, "Clients",    callback=self.get_clients,    autoDefault=False)
        self.btn_ver = gui.button(None, self, "Version",    callback=self.get_version,    autoDefault=False)
        self.btn_inf = gui.button(None, self, "Info",       callback=self.get_info,       autoDefault=False)
        self.btn_img = gui.button(None, self, "Images",     callback=self.get_images,     autoDefault=False)
        self.btn_con = gui.button(None, self, "Containers", callback=self.get_containers, autoDefault=False)
        self.btn_vol = gui.button(None, self, "Volumes",    callback=self.get_volumes,    autoDefault=False)
        self.btn_cli.setFixedWidth(100)
        self.btn_ver.setFixedWidth(100)
        self.btn_inf.setFixedWidth(100)
        self.btn_img.setFixedWidth(100)
        self.btn_con.setFixedWidth(100)
        self.btn_vol.setFixedWidth(100)

        self.hboxTop = QtGui.QHBoxLayout()
        self.hboxTop.setAlignment(QtCore.Qt.AlignLeft)
        self.hboxTop.addWidget(self.btn_con)
        self.hboxTop.addWidget(self.btn_img)
        self.hboxTop.addWidget(self.btn_vol)
        self.hboxTop.addStretch()
        self.hboxTop.addWidget(self.btn_cli)
        self.hboxTop.addWidget(self.btn_ver)
        self.hboxTop.addWidget(self.btn_inf)

        self.infoTable = QtGui.QTableWidget()
        self.infoTable.horizontalHeader().setStretchLastSection(True)
        self.infoTable.setFont(self.font)

        self.btn_rmv = gui.button(None, self, "Remove",  callback=self.remove_selected,  autoDefault=False)
        self.btn_stp = gui.button(None, self, "Stop",    callback=self.stop_selected,    autoDefault=False)
        self.btn_pau = gui.button(None, self, "Pause",   callback=self.pause_selected,   autoDefault=False)
        self.btn_unp = gui.button(None, self, "Unpause", callback=self.unpause_selected, autoDefault=False)
        self.btn_rmv.setHidden(True)
        self.btn_stp.setHidden(True)
        self.btn_pau.setHidden(True)
        self.btn_unp.setHidden(True)
        self.btn_rmv.setFixedWidth(100)
        self.btn_stp.setFixedWidth(100)
        self.btn_pau.setFixedWidth(100)
        self.btn_unp.setFixedWidth(100)

        self.rm_force = QtGui.QCheckBox("Force")
        self.rm_force.setHidden(True)

        self.hbox = QtGui.QHBoxLayout()
        self.hbox.addWidget(self.btn_rmv)
        self.hbox.addSpacing(5)
        self.hbox.addWidget(self.rm_force)
        self.hbox.addSpacing(20)
        self.hbox.addWidget(self.btn_stp)
        self.hbox.addWidget(self.btn_pau)
        self.hbox.addWidget(self.btn_unp)
        self.hbox.setAlignment(QtCore.Qt.AlignLeft)

        self.controlArea.layout().addLayout(self.hboxTop)
        self.controlArea.layout().addWidget(self.infoTable)
        self.controlArea.layout().addLayout(self.hbox)

        self.controlArea.setMinimumSize(900,298)
        self.get_containers()

    """
    Button Callbacks
    """
    def get_clients(self):
        self.context = 'clients'
        self.btn_cli.setDefault(True)
        self.btn_rmv.setHidden(True)
        self.rm_force.setHidden(True)
        self.btn_stp.setHidden(True)
        self.btn_pau.setHidden(True)
        self.btn_unp.setHidden(True)
        self.draw_clients_table()

    def get_version(self):
        self.context = 'version'
        self.btn_ver.setDefault(True)
        self.btn_rmv.setHidden(True)
        self.rm_force.setHidden(True)
        self.btn_stp.setHidden(True)
        self.btn_pau.setHidden(True)
        self.btn_unp.setHidden(True)
        self.draw_version_table()

    def get_info(self):
        self.context = 'info'
        self.btn_inf.setDefault(True)
        self.btn_rmv.setHidden(True)
        self.rm_force.setHidden(True)
        self.btn_stp.setHidden(True)
        self.btn_pau.setHidden(True)
        self.btn_unp.setHidden(True)
        self.draw_info_table()

    def get_images(self):
        self.context = 'images'
        self.btn_img.setDefault(True)
        self.btn_rmv.setHidden(False)
        self.rm_force.setHidden(False)
        self.btn_stp.setHidden(True)
        self.btn_pau.setHidden(True)
        self.btn_unp.setHidden(True)
        self.draw_images_table()

    def get_containers(self):
        self.context = 'containers'
        self.btn_con.setDefault(True)
        self.btn_rmv.setHidden(False)
        self.rm_force.setHidden(False)
        self.btn_stp.setHidden(False)
        self.btn_pau.setHidden(False)
        self.btn_unp.setHidden(False)
        self.draw_containers_table()

    def get_volumes(self):
        self.context = 'volumes'
        self.btn_vol.setDefault(True)
        self.btn_rmv.setHidden(False)
        self.rm_force.setHidden(True)
        self.btn_stp.setHidden(True)
        self.btn_pau.setHidden(True)
        self.btn_unp.setHidden(True)
        self.draw_volumes_table()

    """
    Draw Tables
    """
    def draw_clients_table(self):
        self.infoTable.clearSpans()
        self.infoTable.setColumnCount(3)
        self.infoTable.setHorizontalHeaderLabels(['', 'NAME', 'URL'])
        self.infoTable.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        self.infoTable.horizontalHeader().setVisible(True)
        self.infoTable.setRowCount(len(self.docker_clients))
        for i, client in enumerate(self.docker_clients):
            self.infoTable.setCellWidget(i, 0, QtGui.QCheckBox())
            self.infoTable.setItem(i, 0, QtGui.QTableWidgetItem(''))
            self.infoTable.setItem(i, 1, QtGui.QTableWidgetItem(client.getName()))
            self.infoTable.setItem(i, 2, QtGui.QTableWidgetItem(client.getUrl()))

    def draw_version_table(self):
        self.infoTable.clearSpans()
        self.infoTable.setColumnCount(2)
        self.infoTable.setHorizontalHeaderLabels(['', ''])
        self.infoTable.horizontalHeader().setVisible(False)
        self.infoTable.verticalHeader().setVisible(False)
        totalRows = len(self.docker_clients) + sum(len(client.version().items()) for client in self.docker_clients)
        self.infoTable.setRowCount(totalRows)
        for i, client in enumerate(self.docker_clients):
            clientSep = QtGui.QTableWidgetItem(client.getName() + ": " + client.getUrl())
            clientSep.setBackgroundColor(QtGui.QColor('lightGray'))
            self.infoTable.setItem(i, 0, clientSep)
            self.infoTable.setSpan(i, 0, 1, 2)
            self.infoTable.setCellWidget(i, 0, None)
            table = [(str(k), str(v)) for k, v in sorted(client.version().items())]
            for x, row in enumerate(table):
                self.infoTable.setItem(i+x+1, 0, QtGui.QTableWidgetItem(row[0]))
                self.infoTable.setItem(i+x+1, 1, QtGui.QTableWidgetItem(row[1]))
                self.infoTable.setCellWidget(i+x+1, 0, None)
        self.infoTable.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)

    def draw_info_table(self):
        self.infoTable.clearSpans()
        table = [(str(k), str(v)) for k, v in sorted(self.docker_clients[0].info().items())]
        self.infoTable.setColumnCount(2)
        self.infoTable.setHorizontalHeaderLabels(['', ''])
        self.infoTable.setRowCount(len(table))
        self.infoTable.horizontalHeader().setResizeMode(QtGui.QHeaderView.Interactive)
        self.infoTable.horizontalHeader().resizeSection(0, 150)
        self.infoTable.horizontalHeader().setVisible(False)
        for i, row in enumerate(table):
            self.infoTable.setItem(i, 0, QtGui.QTableWidgetItem(row[0]))
            self.infoTable.setItem(i, 1, QtGui.QTableWidgetItem(row[1]))
            self.infoTable.setCellWidget(i, 0, None)

    def draw_images_table(self):
        self.infoTable.clearSpans()
        tableHeader = ['', 'REPOSITORY', 'TAG', 'ID', 'CREATED', 'SIZE']
        self.infoTable.setColumnCount(len(tableHeader))
        self.infoTable.setHorizontalHeaderLabels(tableHeader)
        self.infoTable.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        self.infoTable.horizontalHeader().setVisible(True)
        table = []
        imgs = self.docker_clients[0].images()
        if not type(imgs) is list:
            self.infoTable.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
            self.infoTable.horizontalHeader().setResizeMode(0, QtGui.QHeaderView.ResizeToContents)
        else:
            for img in imgs:
                id = img['Id'].split(':')[1][:12]
                created = millis_to_datetime(img['Created'])
                size = sizeof_fmt(img['Size'])
                repo_info = [id, created, size]
                for repo in img['RepoTags']:
                    table.append(repo.split(':') + repo_info)

        self.infoTable.setRowCount(len(table))
        for i, row in enumerate(table):
            self.infoTable.setCellWidget(i, 0, QtGui.QCheckBox())
            self.infoTable.setItem(i, 0, QtGui.QTableWidgetItem(''))
            for x, cell in enumerate(row):
                self.infoTable.setItem(i, 1+x, QtGui.QTableWidgetItem(cell))

    def draw_containers_table(self):
        self.infoTable.clearSpans()
        tableHeader = ['', 'CONTAINER ID', 'IMAGE', 'COMMAND', 'CREATED', 'STATUS', 'PORTS', 'NAMES']
        self.infoTable.setColumnCount(len(tableHeader))
        self.infoTable.setHorizontalHeaderLabels(tableHeader)
        self.infoTable.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        self.infoTable.horizontalHeader().setVisible(True)
        self.infoTable.horizontalHeader().setResizeMode(3, QtGui.QHeaderView.Interactive)
        self.infoTable.horizontalHeader().resizeSection(3, 120)
        conts = self.docker_clients[0].containers()
        if not type(conts) is list:
            self.infoTable.setRowCount(0)
            self.infoTable.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
            self.infoTable.horizontalHeader().setResizeMode(0, QtGui.QHeaderView.ResizeToContents)
        else:
            self.infoTable.setRowCount(len(conts))
            for i, cont in enumerate(conts):
                self.infoTable.setCellWidget(i, 0, QtGui.QCheckBox())
                self.infoTable.setItem(i, 0, QtGui.QTableWidgetItem(''))
                self.infoTable.setItem(i, 1, QtGui.QTableWidgetItem(cont['Id'][:12]))
                img_split = cont['Image'].split(':')
                img_str = img_split[0] if len(img_split) < 2 else img_split[1][:12]
                self.infoTable.setItem(i, 2, QtGui.QTableWidgetItem(img_str))
                self.infoTable.setItem(i, 3, QtGui.QTableWidgetItem(cont['Command']))
                self.infoTable.setItem(i, 4, QtGui.QTableWidgetItem(millis_to_datetime(cont['Created'])))
                self.infoTable.setItem(i, 5, QtGui.QTableWidgetItem(cont['Status']))
                self.infoTable.setItem(i, 6, QtGui.QTableWidgetItem('\n'.join([str(_.get('PublicPort', ' ')) for _ in cont['Ports']])))
                self.infoTable.setItem(i, 7, QtGui.QTableWidgetItem('\n'.join(cont['Names'])))

    def draw_volumes_table(self):
        self.infoTable.clearSpans()
        tableHeader = ['', 'DRIVER', 'VOLUME NAME', 'MOUNT POINT']
        self.infoTable.setColumnCount(len(tableHeader))
        self.infoTable.setHorizontalHeaderLabels(tableHeader)
        self.infoTable.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        self.infoTable.horizontalHeader().setVisible(True)
        vols = self.docker_clients[0].volumes()
        if not type(vols) is list:
            self.infoTable.setRowCount(0)
            self.infoTable.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
            self.infoTable.horizontalHeader().setResizeMode(0, QtGui.QHeaderView.ResizeToContents)
        else:
            self.infoTable.setRowCount(len(vols))
            for i, vol in enumerate(vols):
                self.infoTable.setCellWidget(i, 0, QtGui.QCheckBox())
                self.infoTable.setItem(i, 0, QtGui.QTableWidgetItem(''))
                self.infoTable.setItem(i, 1, QtGui.QTableWidgetItem(vol['Driver']))
                self.infoTable.setItem(i, 2, QtGui.QTableWidgetItem(vol['Name']))
                self.infoTable.setItem(i, 3, QtGui.QTableWidgetItem(vol['Mountpoint']))

    """
    Remove functions
    """
    def remove_image(self, id, force):
        print('Removing image: ' + id)
        self.docker_clients[0].remove_image(id, force)

    def remove_container(self, id, force):
        print('Removing container: ' + id)
        self.docker_clients[0].remove_container(id, force)

    def remove_volume(self, id, force):
        print('Removing volume: ' + id)
        self.docker_clients[0].remove_volume(id)

    def remove_selected(self):
        remove_func = None
        draw_func = None
        id_column = 0
        if (self.context == 'images'):
            remove_func = self.remove_image
            draw_func = self.draw_images_table
            id_column = 3
        elif (self.context == 'containers'):
            remove_func = self.remove_container
            draw_func = self.draw_containers_table
            id_column = 1
        elif (self.context == 'volumes'):
            remove_func = self.remove_volume
            draw_func = self.draw_volumes_table
            id_column = 2
        for i in range(0,self.infoTable.rowCount()):
            if (self.infoTable.cellWidget(i,0).isChecked()):
                remove_func(self.infoTable.item(i,id_column).text(), self.rm_force.isChecked())
        draw_func()

    """
    Stop, Pause, Unpause only apply to containers
    """
    def stop_selected(self):
        for i in range(0,self.infoTable.rowCount()):
            if (self.infoTable.cellWidget(i,0).isChecked()):
                status = self.infoTable.item(i,5).text()
                if (status.startswith('Up') and not status.endswith('(Paused)')):
                    self.docker_clients[0].stop_container(self.infoTable.item(i,1).text())
        self.draw_containers_table()

    def pause_selected(self):
        for i in range(0,self.infoTable.rowCount()):
            if (self.infoTable.cellWidget(i,0).isChecked()):
                status = self.infoTable.item(i,5).text()
                if (status.startswith('Up') and not status.endswith('(Paused)')):
                    self.docker_clients[0].pause_container(self.infoTable.item(i,1).text())
        self.draw_containers_table()

    def unpause_selected(self):
        for i in range(0,self.infoTable.rowCount()):
            if (self.infoTable.cellWidget(i,0).isChecked()):
                status = self.infoTable.item(i,5).text()
                if (status.startswith('Up') and status.endswith('(Paused)')):
                    self.docker_clients[0].unpause_container(self.infoTable.item(i,1).text())
        self.draw_containers_table()


def table_str(table):
    str_table = ''
    col_width = [max(len(x) for x in col) for col in zip(*table)]
    for line in table:
        str_table += "   ".join("{:{}}".format(x, col_width[i]) for i, x in enumerate(line))
        str_table += '\n'

    return str_table

def sizeof_fmt(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%6.2f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)

def millis_to_datetime(millis):
    return datetime.datetime.fromtimestamp(millis).strftime('%Y-%m-%d %H:%M:%S')