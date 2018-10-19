import os
import re
import sys
import json
import jsonpickle
import pickle
import csv
import tempfile,shutil
import OWImageBuilder
import workflowTools
from xml.dom import minidom
from glob import glob
from pathlib import Path
from shutil import copyfile
from createWidget import mergeWidget, createWidget,findIconFile
from copy import deepcopy
from collections import OrderedDict
from functools import partial
from AnyQt.QtCore import QThread, pyqtSignal, Qt
from Orange.widgets import widget, gui, settings
from DockerClient import DockerClient, PullImageThread, ConsoleProcess
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWidgets import QInputDialog, QLineEdit, QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget, QVBoxLayout, QMessageBox


from AnyQt.QtWidgets import (
    QWidget, QButtonGroup, QGroupBox, QRadioButton, QSlider,
    QDoubleSpinBox, QComboBox, QSpinBox, QListView, QLabel,
    QScrollArea, QVBoxLayout, QHBoxLayout, QFormLayout,
    QSizePolicy, QApplication, QCheckBox
)
from makeToolDockCategories import *

defaultIconFile='/icons/default.png'

def categoryWidgetPath(categoryDir):
    pyFiles=glob('/biodepot/{}/OW*.py'.format(categoryDir))
    if not pyFiles:
        return None
    widgetPath=os.path.dirname((os.readlink(pyFiles[0])))
   #check if absolute path
    if widgetPath[0]=='/':
        absWidgetPath=widgetPath
    else:
        absWidgetPath=os.path.abspath('/biodepot/{}/{}'.format(categoryDir,widgetPath))
    return os.path.dirname(absWidgetPath)
    

def registerDirectory(baseToolPath):
    os.system('cd {} && pip install -e .'.format(baseToolPath))

class ToolDockEdit(widget.OWWidget):
    name = "ToolDockEditor"
    description = "Edits contents of tool dock"
    category = "Utility"
    icon = "icons/build.png"
    priority = 2
    inputs = []
    outputs = []
    pset=partial(settings.Setting,schema_only=True)
    
    want_main_area = False
    want_control_area = True
    allStates=pset({})
    allAttrs=pset({})
    data={}
    
    def __init__(self,renameData=None,canvasMainWindow=None):
        super().__init__()
        self.canvas=canvasMainWindow
        if renameData:
            #run the rename routine and then finish
            self.updateCategories()
            self.widgetRename(widgetName=renameData['widgetName'],category=renameData['category'],newName=renameData['newName'])
            return  
        self.css = '''
        QPushButton {background-color: #1588c5; color: white; height: 20px; border: 1px solid black; border-radius: 2px;}
        QPushButton:hover {background-color: #1555f5; }
        QPushButton:hover:pressed { background-color: #1588c5; color: black; border-style: inset; border: 1px solid white} 
        QPushButton:disabled { background-color: lightGray; border: 1px solid gray; }
        '''
        self.baseToolPath='/biodepot'
        self.setStyleSheet(self.css)
        self.browseIcon=QtGui.QIcon('/icons/bluefile.png')
        self.addIcon=QtGui.QIcon('/icons/add.png')
        self.removeIcon=QtGui.QIcon('/icons/remove.png')
        self.submitIcon=QtGui.QIcon('/icons/submit.png')
        self.reloadIcon=QtGui.QIcon('/icons/reload.png')
        self.controlArea.setMinimumWidth(500)
        self.controlArea.setMinimumHeight(120)
        self.startWidget()
    
    def clearLayout(self,layout):
        if layout != None:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget() is not None:
                    child.widget().deleteLater()
                elif child.layout() is not None:
                    self.clearLayout(child.layout())
                    
    def addWidgetToCategory (self,category,directory,widgetPath,confirmation=True):
        qm = QtGui.QMessageBox
        widgetName=os.path.basename(widgetPath)
        destLink='/biodepot/{}/OW{}.py'.format(directory,widgetName)
        pythonFile=glob('{}/*.py'.format(widgetPath))[0]
        destDir=categoryWidgetPath(directory)
        if not destDir:
            return
        destPath=destDir+'/'+ widgetName
        if os.path.exists(destPath) or os.path.exists(destLink):
            ret=qm.question(self,'', "symlink {} or destination widget {} exists - OverWrite ?".format(destLink,destPath), qm.Yes | qm.No)
            if ret == qm.No:
                return
            if os.path.exists(destPath):
                shutil.rmtree(destPath)
        #try:
        shutil.copytree(widgetPath,destPath)
        os.system('ln -sf {} {}'.format(pythonFile,destLink))
        if confirmation:
            widgetName=os.path.basename(widgetPath)
            title='Add {}'.format(widgetName)
            message='Added {} to {} in ToolDock'.format(widgetName,category)
            qm.information(self,title,message,QtGui.QMessageBox.Ok)
            self.canvas.reload_current()
        #except:
            #pass
    
    def startWidget(self):
        self.setWindowTitle('WidgetToolDock Editor')
        self.grid=QtGui.QGridLayout()
        self.controlArea.layout().addLayout(self.grid)
        self.initCategories()
        self.drawAddWidget()
        self.drawRemoveWidget()
        self.drawRenameWidget()
        self.drawRemoveCategory()
        
        #self.setStyleSheet(":disabled { color: #282828}")
    
    def updateCategories(self):
        self.categories=(str(os.popen('''grep -oP 'name="\K[^"]+' /biodepot/setup.py''').read())).split()
        #directories are not same as categories because Python/Linux unfriendly characters are changed
        self.directoryList=(str(os.popen('''grep -oP 'packages=\["\K[^"]+' /biodepot/setup.py''').read())).split()
        self.categoryToDirectory={}
        self.categoryToPath={}
        #check for icon link for categories where there are no widgets
        for index, category in enumerate(self.categories):
            self.categoryToDirectory[category]=self.directoryList[index]
            iconLink='/biodepot/{}/icon'.format(self.categoryToDirectory[category])
            if os.path.islink(iconLink):
                iconPath=os.readlink(iconLink)
                self.categoryToPath[category]=os.path.dirname(os.path.dirname(os.path.dirname(iconPath)))    
            else:
                #this is one of biodepots native categories if the icon directory is not a link
                self.categoryToPath[category]='/widgets/{}'.format(self.categoryToDirectory[category])
                
    def updateWidgetList(self):
        self.updateCategories()
        self.widgetList={}
        for category in self.categories:
            directory = self.categoryToDirectory[category]
            self.widgetList[category]=[os.path.basename(x)[2:-3] for x in glob('/biodepot/{}/OW*.py'.format(directory))]
        
    def initCategories(self):
        self.updateCategories()
        self.clearLayout(self.controlArea.layout())
        
    def drawAddWidget(self):
        ledit=self.makeLedit(self.grid,'Enter widget name','Add widget ',startRow=1,startColumn=1,browse=True)
        cbox=self.makeComboBox(self.grid,' to drawer:',self.categories,startRow=1,startColumn=4)
        widgetAddBtn = gui.button(None, self, "Add", callback= lambda: self.widgetAdd(ledit,cbox))
        widgetAddBtn.setFixedSize(30,20)
        widgetAddBtn.setStyleSheet(self.css)
        self.grid.addWidget(widgetAddBtn,1,6)        
        self.controlArea.layout().addLayout(self.grid)

    def widgetAdd(self,ledit,cbox,widgetsDir='/widgets'):
        category=self.getComboValue(cbox)
        inputDir=self.getLeditValue(ledit)
        if not inputDir:
            return
        directory=self.categoryToDirectory[category] 
        self.addWidgetToCategory(category,directory,inputDir)
    
    def drawRenameWidget(self):
        self.updateWidgetList()        
        self.RNWwbox=self.makeComboBox(self.grid,'Rename widget ',self.widgetList[self.categories[0]],startRow=3)
        self.RNWcbox=self.makeComboBox(self.grid,' in drawer ',self.categories,startRow=3,startColumn=4,callback=lambda: self.__onRWCategoryChoice(self.RNWcbox,self.RNWwbox))
        self.RNWledit=self.makeLedit(self.grid,'New name','to new name',startRow=3,startColumn=6,browse=False)
        widgetRenameBtn = gui.button(None, self, "Rename", callback=self.widgetRename)
        widgetRenameBtn.setFixedSize(60,20)
        widgetRenameBtn.setStyleSheet(self.css)
        self.grid.addWidget(widgetRenameBtn,3,8)
            
    def drawRemoveWidget(self):
        self.updateWidgetList()        
        self.RWwbox=self.makeComboBox(self.grid,'Remove widget ',self.widgetList[self.categories[0]],startRow=2)
        self.RWcbox=self.makeComboBox(self.grid,' from drawer ',self.categories,startRow=2,startColumn=4,callback=lambda: self.__onRWCategoryChoice(self.RWcbox,self.RWwbox))
        widgetRemoveBtn = gui.button(None, self, "Remove", callback=self.widgetRemove)
        widgetRemoveBtn.setFixedSize(60,20)
        widgetRemoveBtn.setStyleSheet(self.css)
        self.grid.addWidget(widgetRemoveBtn,2,6)
        
    def widgetRemove(self):
        qm = QtGui.QMessageBox
        widgetName=self.getComboValue(self.RWwbox)
        category=self.getComboValue(self.RWcbox)
        categoryDir=self.categoryToDirectory[category]
        symlink='/biodepot/{}/OW{}.py'.format(categoryDir,widgetName)
        if not os.path.exists(symlink):
            qm.information(self,'Non-existent symlink','Symlink {} not found'.format(symlink),QtGui.QMessageBox.Ok)
            return
        widgetPath=os.path.dirname(os.readlink(symlink))
        #check if relative path and convert to absolute path if necessary
        if widgetPath[0] != '/':
            widgetPath=os.path.abspath('/biodepot/{}/{}'.format(categoryDir,widgetPath))
        if not os.path.exists(widgetPath):
            qm.information(self,'Non-existent widget','Widget {} not found'.format(widgetPath),QtGui.QMessageBox.Ok)
            return
        #check if any ows file that is in toolbox uses the widget and issue a warning if so
        workflowPaths=self.findWorkflowsWithWidget(widgetName,category)
        print('workflow paths are')
        print(workflowPaths)
        if workflowPaths:
            workflowStr=','.join(workflowPaths)
            ret=qm.question(self,'', "Workflows {} use this widget - Do you still want to remove it?".format(workflowStr), qm.Yes | qm.No)
            if ret == qm.No:
                return
            for workflowPath in workflowPaths:
                workflowTools.removeWidgetfromWorkflow(workflowPath,workflowPath,category,widgetName)
        #finally remove the widget and symlink
        shutil.rmtree(widgetPath)
        os.unlink(symlink)
        qm.information(self,'Successfully removed','removed widget {} from {}'.format(widgetPath,category),QtGui.QMessageBox.Ok)
        self.canvas.reload_current()
    
    def widgetRename(self,widgetName=None,category=None,newName=None):
        qm = QtGui.QMessageBox
        if not widgetName:
            widgetName=self.getComboValue(self.RNWwbox)
        if not category:
            category=self.getComboValue(self.RNWcbox)
        categoryDir=self.categoryToDirectory[category]
        symlink='/biodepot/{}/OW{}.py'.format(categoryDir,widgetName)
        if not newName:
            newName=self.getLeditValue(self.RNWledit)
        if not newName or newName == widgetName:
            return
        #check if name already exists
        if not os.path.exists(symlink):
            qm.information(self,'','Non-existent symlink','Symlink {} not found'.format(symlink),QtGui.QMessageBox.Ok)
            return
        widgetPath=os.path.dirname(os.path.realpath(symlink))
        newWidgetPath=os.path.dirname(widgetPath)+'/'+newName
        if os.path.exists(newWidgetPath):
            qm.information(self,'','widget {} already exists'.format(newWidgetPath),QtGui.QMessageBox.Ok)
            return 
        if not os.path.exists(widgetPath):
            qm.information(self,'Non-existent widget','Widget {} not found'.format(widgetPath),QtGui.QMessageBox.Ok)
            return
        #check if any ows file that is in toolbox uses the widget and issue a warning if so
        print ('widgetName is '+ widgetName)
        print ('category  is' + category)
        workflowPaths=self.findWorkflowsWithWidget(widgetName,category)
        print ('workflowPaths are ')
        print (workflowPaths)
        if workflowPaths:
            myWorkflowNames=[]
            for path in workflowPaths:
                myWorkflowNames.append(os.path.basename(os.path.dirname(path)))
            workflowStr=','.join(myWorkflowNames)
            ret=qm.question(self,'', "Workflows {} use this widget - Do you still want to rename it?".format(workflowStr), qm.Yes | qm.No)
            if ret == qm.No:
                return
        for workflowPath in workflowPaths:
            workflowTools.renameWidgetInWorkflow(workflowPath,workflowPath,category,widgetName,newName)
        workflowTools.renameWidget(widgetPath,widgetName,newName)
        #finally remove the widget and symlink
        qm.information(self,'Successfully renamed','renamed widget {} to {}'.format(widgetName,newName),QtGui.QMessageBox.Ok)
        self.canvas.reload_current()
                
    def findWorkflowsWithWidget(self,widgetName,widgetCategory):
        workflowPaths=[]
        for category in self.categories:
            print(category + ' ' +self.categoryToPath[category])
            owsFiles=glob('{}/*.ows'.format(self.categoryToPath[category]))
            for owsFile in owsFiles:
                print ('checking ' + owsFile + ' for '+ widgetCategory +' '+widgetName) 
                doc = minidom.parse(owsFile)
                nodes = doc.getElementsByTagName("node")
                for node in nodes:
                    print('node pn {} node name {} category {} wname {}'.format(node.getAttribute('project_name'),node.getAttribute('name'),category,widgetName))
                    if node.getAttribute('project_name') == widgetCategory and node.getAttribute('name') == widgetName:
                        workflowPaths.append(owsFile)
                        continue
        return workflowPaths
    
    def drawRemoveCategory(self):
        cbox=self.makeComboBox(self.grid,'Remove drawer:',self.categories,startRow=4,startColumn=1)
        catRemoveBtn = gui.button(None, self, "Remove", callback= lambda: self.categoryRemove(cbox))
        catRemoveBtn.setFixedSize(60,20)
        catRemoveBtn.setStyleSheet(self.css)
        self.grid.addWidget(catRemoveBtn,4,3)
                     
    def categoryAdd(self,nameLedit,iconLedit):
        qm = QtGui.QMessageBox
        category=self.getLeditValue(nameLedit)
        if not category:
            return
        if category in self.categories:
            qm.information(self,'Add drawer','Drawer {} already in ToolDock'.format(category),QtGui.QMessageBox.Ok)
            return
        iconFile=self.getLeditValue(iconLedit)
        directory=addCategoryToToolBox(self.baseToolPath,category,iconFile=iconFile)
        registerDirectory(self.baseToolPath)
        qm.information(self,'Add drawer','Added drawer {} to directory {} in ToolDock'.format(category,directory),QtGui.QMessageBox.Ok)
        self.updateCategories()
        self.canvas.reload_current()
           
    def categoryRemove(self,nameCombo):
        qm = QtGui.QMessageBox
        category=self.getComboValue(nameCombo)
        if not category:
            return
        self.updateCategories()
        if(category not in self.categories):
            qm.information(self,'Remove drawer', 'Drawer {} not in ToolDock'.format(category),QtGui.QMessageBox.Ok)
            return
        ret=qm.question(self,'', "Remove {} from ToolDock ?".format(category), qm.Yes | qm.No)
        if ret == qm.No:
            return
        directory=self.categoryToDirectory[category]
        removeCategoryFromToolDock(self.baseToolPath,category,self.categoryToDirectory[category])
        registerDirectory(self.baseToolPath)
        self.updateCategories()
        if self.canvas.recent_schemes:
            currentDirectory=os.path.basename(os.path.dirname(self.canvas.recent_schemes[0][1]))
            if currentDirectory != directory:
                self.canvas.reload_current()
        self.canvas.reload_settings()
      
    def getCategoryList(self,widgetName):
        #categories may not be directories
        if not widgetName:
            return self.categories
        returnList=[]
        for category in self.categories:
            directory=self.categoryToDirectory[category]
            if os.path.exists('{}/{}/OW{}.py'.format(self.baseToolPath,directory,widgetName)):
                returnList.append(category)
        return returnList

    def makeLedit(self,layout,text=None,label=None,startRow=1,startColumn=1,browse=False,browseFileFlag=False):
        leditLabel=None
        if(label):
            leditLabel=QtGui.QLabel(label)
        ledit = QtGui.QLineEdit(self)
        ledit.setClearButtonEnabled(True)
        ledit.setPlaceholderText(text)
        ledit.setStyleSheet(":disabled { color: #282828}")           
        layout.addWidget(leditLabel,startRow,startColumn)
        layout.addWidget(ledit,startRow,startColumn+1)
        if browse:
            button=gui.button(None, self, "",callback= lambda: self.browseWidget(ledit,browseFileFlag),autoDefault=True, width=19, height=19)
            button.setIcon(self.browseIcon)
            layout.addWidget(button,startRow,startColumn+2)
        return ledit
        
    def browseWidget(self, ledit, browseFileFlag=False):
        if browseFileFlag:
            ledit.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, "Locate file", '/widgets')[0]))
            return
        myFileDir=QtWidgets.QFileDialog.getExistingDirectory(self, caption="Locate widget", directory='/widgets')
        ledit.setText(str(myFileDir))
            
    def makeComboBox (self,layout,label, elements,startRow=1,startColumn=1,callback=None):
        comboBoxLabel=QtGui.QLabel(label)
        comboBox=QtGui.QComboBox()
        if elements:
            comboBox.addItems(elements)
        comboBox.currentIndex=0
        if callback:
            comboBox.currentIndexChanged.connect(callback) 
        layout.addWidget(comboBoxLabel,startRow,startColumn)
        layout.addWidget(comboBox,startRow,startColumn+1)
        return comboBox
    
    def __onRWCategoryChoice(self,cbox,wbox):
        wbox.clear()
        wbox.addItems(self.widgetList[cbox.currentText()])
        
    def getComboValue(self,comboBox):
        if comboBox.isEnabled():
            return comboBox.currentText()  
        return None 
        
    def getLeditValue(self,ledit):
        if ledit.isEnabled():
            return ledit.text()
        return None
