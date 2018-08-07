import os
import re
import sys
import json
import jsonpickle
import pickle
import csv
from pathlib import Path
from orangebiodepot.util.createWidget import mergeWidget, createWidget, findDirectory, findIconFile
from copy import deepcopy
from collections import OrderedDict
from functools import partial
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient, PullImageThread, ConsoleProcess
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import QInputDialog, QLineEdit, QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget, QVBoxLayout

from AnyQt.QtWidgets import (
    QWidget, QButtonGroup, QGroupBox, QRadioButton, QSlider,
    QDoubleSpinBox, QComboBox, QSpinBox, QListView, QLabel,
    QScrollArea, QVBoxLayout, QHBoxLayout, QFormLayout,
    QSizePolicy, QApplication, QCheckBox
)

from collections import defaultdict
from pathlib import Path
nested_dict = lambda: defaultdict(nested_dict)




class OWWidgetBuilder(widget.OWWidget):
    name = "Widget Builder"
    description = "Build a new widget from a set of bash commands and a Docker container"
    category = "Bwb-core"
    icon = "icons/build.png"
    priority = 2

    inputs = []
    outputs = []
    pset=partial(settings.Setting,schema_only=True)
    
    want_main_area = False
    want_control_area = True
    allStates=pset({})
    allAttrs=pset({})
    
    def place_in_nested_dict(self,nested, place, value):
        if isinstance(place, str):
            place = [place]
        place = list(place)
        last = place.pop()
        for node in reversed(place):
            nested = nested[node]
        nested[last] = value
    
    def find_files(self, root=Path('.'), pattern='**/*.ows'):
        tree = nested_dict()
        for file in root.glob(pattern):
            if file.is_dir():
                continue
            parts = file.relative_to(root).parts
            name = file.name
            self.place_in_nested_dict(tree, parts, file.relative_to(root))
        return tree
    
    def populate(self, tree: dict, root: Path):
        tree_item = QTreeWidgetItem()
        tree_item.setText(0, str(root))
        for key, value in tree:
            if isinstance(value, dict):
                tree_item.addChild(populate(tree, key))
            else:
                tree_item.addChild(key)
        return tree_item
    
    def __init__(self):
        super().__init__()
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid #1a8ac6; border-radius: 2px;}
        QPushButton:pressed { background-color: #158805; border-style: inset;}       QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        QPushButton:hover {background-color: #1588f5; }
        '''  
        self.setStyleSheet(css)
        tree=self.find_files
        self.controlArea.layout().addWidget(tree)

