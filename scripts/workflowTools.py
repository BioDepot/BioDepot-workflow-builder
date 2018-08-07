import os
import re
import sys
import json
import jsonpickle
import pickle
import csv
from pathlib import Path
from copy import deepcopy
from collections import OrderedDict
from functools import partial
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import *

from AnyQt.QtWidgets import (
    QWidget, QButtonGroup, QGroupBox, QRadioButton, QSlider,
    QDoubleSpinBox, QComboBox, QSpinBox, QListView, QLabel,
    QScrollArea, QVBoxLayout, QHBoxLayout, QFormLayout,
    QSizePolicy, QApplication, QCheckBox
)

from collections import defaultdict
from pathlib import Path
nested_dict = lambda: defaultdict(nested_dict)

class WorkflowList():
    name = "Widget Editor"
    icon = "icons/build.png"
    
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
            print (file)
        return tree
    
    def __init__(self):
#        super().__init__()
        css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid #1a8ac6; border-radius: 2px;}
        QPushButton:pressed { background-color: #158805; border-style: inset;}       QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        QPushButton:hover {background-color: #1588f5; }
        '''  


wlist=WorkflowList()
tree=wlist.find_files(root=Path('../workflows'))
