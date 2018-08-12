#!/usr/bin/env python3
from subprocess import call
from shutil import copyfile
import sys, re, os, getopt
from os import listdir
from os.path import isfile, join
from collections import OrderedDict
import jsonpickle
import pprint
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import QInputDialog, QLineEdit

def deClass(string):
    #removes the <class 'id'> and returns id
    if string[:6] == '<class':
        m=re.match( r"(\<class ')(\w+)(')",string)
        try:
            return m.group(2)
        except:
            return string
    return string
    
def findIconFile(widgetFile):
    with open(widgetFile,'r') as f:
        data=f.read()
    f.close()
    return os.path.basename(re.findall(r'icon = "(.*?)"',data)[0])

WIDGET_HEADING ='''import os
import glob
import sys
import functools
import jsonpickle
from collections import OrderedDict
from Orange.widgets import widget, gui, settings
import Orange.data
from Orange.data.io import FileFormat
from DockerClient import DockerClient
from BwBase import OWBwBWidget, ConnectionDict, BwbGuiElements
from PyQt5 import QtWidgets, QtGui

'''

def findDirectory(inputJson):
    with open(inputJson) as f:
        data=jsonpickle.decode(f.read())
    return checkCategory(data['category'])
      
def checkCategory(category):
    categories=(str(os.popen('''grep -oP 'name="\K[^"]+' /biodepot/setup.py''').read())).split()
    if category in categories:
        #has/setattr does not like -
        directory=category.replace('-','_')
        return directory
    else:
        sys.stderr.write('*WARNING* {} not a recognized category - will place widget in User directory\n'.format(category))
        return 'User'
        
def createWidget(inputJson,outputWidget,widgetName,inputData=None):
    defaultIconFile='/biodepot/Bwb_core/icons/default.png'
    widgetPath = os.path.dirname(os.path.realpath(outputWidget)) 
    data={}
    directory='User'
    if inputJson:
        with open(inputJson) as f:
            data=jsonpickle.decode(f.read())
    elif inputData:
        data=inputData
        dataJ=jsonpickle.encode(data)
        inputJson=os.path.splitext(outputWidget)[0]+'.json'
        with open(inputJson,"w") as f:
            f.write(dataJ)
    inputPath=os.path.dirname(os.path.realpath(inputJson))
    directory=checkCategory(data['category'])        
    #write preInit
    with open(outputWidget,'w') as f:
        f.write(WIDGET_HEADING)
        className='class OW{}(OWBwBWidget):\n'.format(data['name'].replace(' ',''))
        f.write(className)
        f.write('    name = "{}"\n'.format(data['name']))
        f.write('    description = "{}"\n'.format(data['description']))
        f.write('    category = "{}"\n'.format(data['category']))
        priority=10
        if 'priority' in data and (data['priority'] == 0 or data['priority']):
            priority=data['priority']
        f.write('    priority = {}\n'.format(priority))
        iconFile=data['icon']
        os.system("mkdir -p {}/icon".format(widgetPath))
        if iconFile and os.path.exists(iconFile):
            os.system("rm {}/icon/* ".format(widgetPath))
            os.system("cp {} {}/icon/".format(iconFile,inputPath))
        else:
            icons=os.listdir(inputPath+'/icon')
            iconFile=os.path.basename(defaultIconFile)
            if not icons:
                os.system("cp {} {}/icon/".format(defaultIconFile,inputPath))
            else:
                iconFile=icons[0]
        finalIconFile = '/widgets/' + widgetName + '/icon/' + os.path.basename(iconFile)
        f.write('    icon = "{}"\n'.format(finalIconFile))
        f.write('    want_main_area = False\n')
        if not 'docker_image_name' in data:
            data['docker_image_name']='biodepot/alpine-bash'
        f.write('    docker_image_name = "{}"\n'.format(data['docker_image_name']))
        if not'docker_image_tag' in data:
            data['docker_image_tag'] = 'latest'
        f.write('    docker_image_tag = "{}"\n'.format(data['docker_image_tag']))
        #inputs and outputs
        if 'inputs' in data and data['inputs']:
            inputStr='['
            for attr, values  in data['inputs'].items():
                if 'callback' in  values and values['callback']:
                    inputStr= inputStr+ '("{}",{},"{}"),'.format(attr,deClass(str(values['type'])),values['callback'])
                else:
                    inputStr= inputStr+ '("{}",{},"handleInputs{}"),'.format(attr,deClass(str(values['type'])),attr)
            inputStr=inputStr[:-1]+']'
            f.write('    inputs = {}\n'.format(inputStr))
        
        if 'outputs' in data and data['outputs']:
            outputStr='['
            for attr, value  in data['outputs'].items():
                outputStr= outputStr+ '("{}",{}),'.format(attr,deClass(str(value['type'])))
            outputStr=outputStr[:-1]+']'
            f.write('    outputs = {}\n'.format(outputStr))
        #permanent settings
        f.write('    pset=functools.partial(settings.Setting,schema_only=True)\n')
        f.write('    runMode=pset(0)\n')
        f.write('    exportGraphics=pset(False)\n')
        f.write('    runTriggers=pset([])\n')
        f.write('    triggerReady=pset({})\n')
        f.write('    inputConnectionsStore=pset({})\n')
        f.write('    optionsChecked=pset({})\n')
        if 'parameters' in data and data['parameters']:
            for pname,pvalue in data['parameters'].items():
                if 'default' in pvalue and pvalue['default'] is not None:
                    #if it is not a number or dict type then keep quotes
                    if pvalue['type'] == 'int' or  pvalue['type'] == 'float' or pvalue['type'] == 'double' or pvalue['type'] == 'bool' or pvalue['type'][-4:] == 'dict' or pvalue['type'][-4:] == 'Dict' or pvalue['type'][-4:] == 'list':
                        f.write('    {}=pset({})\n'.format(pname,pvalue['default']))
                    else:
                        f.write('    {}=pset("{}")\n'.format(pname,pvalue['default']))
                else:               
                    if pvalue['type'] == type('str') :
                        f.write('    {}=pset("")\n'.format(pname),None)
                    if pvalue['type'] == 'bool' :
                        f.write('    {}=pset(False)\n'.format(pname))
                    elif pvalue['type'][-4:] == 'list' or pvalue['type'][-4:] == 'List':
                        f.write('    {}=pset([])\n'.format(pname))
                    elif pvalue['type'][-4:] == 'dict' or pvalue['type'][-4:] == 'Dict':
                        f.write('    {}=pset({{}})\n'.format(pname))
                    else:
                        f.write('    {}=pset(None)\n'.format(pname))
        f.write('    def __init__(self):\n')        
        f.write('        super().__init__(self.docker_image_name, self.docker_image_tag)\n')
        joinedName=data['name'].replace(' ','')
        f.write('        with open("/widgets/{}/{}") as f:\n'.format(widgetName,os.path.basename(inputJson)))
        f.write('            self.data=jsonpickle.decode(f.read())\n')
        f.write('            f.close()\n')
        #init
        f.write('        self.initVolumes()\n')
        f.write('        self.inputConnections = ConnectionDict(self.inputConnectionsStore)\n')
        f.write('        self.drawGUI()\n')
        #input callbacks
        if 'inputs' in data and data['inputs']:
            for attr, values in data['inputs'].items():
                f.write('    def handleInputs{}(self, value, *args):\n'.format(attr))
                f.write('        if args and len(args) > 0: \n')
                f.write('            self.handleInputs("{}", value, args[0][0])\n'.format(attr))
                f.write('        else:\n')
                f.write('            self.handleInputs("inputFile", value, None)\n'.format(attr))               
        if 'outputs' in data and data['outputs']:
            #generic omnibus output handler - change if you want to customize outputs
            f.write('    def handleOutputs(self):\n'.format(attr))
            for attr,values in data['outputs'].items():
                if 'default' in values and values['default'] is not None:
                    f.write('        outputValue="{}"\n'.format(values['default']))
                else:
                    f.write('        outputValue=None\n')
                f.write('        if hasattr(self,"{}"):\n'.format(attr))
                f.write('            outputValue=getattr(self,"{}")\n'.format(attr))
                f.write('        self.send("{}", outputValue)\n'.format(attr))
        f.close()

def fwrite(lines, match, replacement):
    myPattern=re.compile(match)
    for  i,line in enumerate(lines):
        if myPattern.search(line):
            lines[i]=replacement
            return
    lines.append(replacement)
        
def mergeWidget(inputJson,outputWidget,widgetName,inputData=None):
    defaultIconFile='/biodepot/Bwb_core/icons/default.png'
    widgetPath = os.path.dirname(os.path.realpath(outputWidget)) 
    data={}
    directory='User'
    if inputJson:
        with open(inputJson) as f:
            data=jsonpickle.decode(f.read())
    elif inputData:
        data=inputData
        dataJ=jsonpickle.encode(data)
        inputJson=os.path.splitext(outputWidget)[0]+'.json'
        with open(inputJson,"w") as f:
            f.write(dataJ)
    directory=checkCategory(data['category'])
    inputPath=os.path.dirname(os.path.realpath(inputJson))
    #updates the python file up to the first def command
    #split the python file into 3
    beforeClass=[]
    afterClass=[]
    afterDef=[]
    
    fpattern = re.compile("def .+:")
    with open(outputWidget,'r') as f:
        for line in f:
            if beforeClass == [] and line.startswith('class OW'):
                beforeClass=afterDef
                afterDef=[]
            elif afterClass == [] and fpattern.search(line):
                afterClass=afterDef
                afterDef=[]
            afterDef.append(line)
    
    #don't replace the header - start with the afterClass
    className='class OW{}(OWBwBWidget):\n'.format(data['name'].replace(' ',''))
    afterClass[0]=className
    fwrite(afterClass,"\s+name = ",'    name = "{}"\n'.format(data['name']))
    fwrite(afterClass,"\s+description =",'    description = "{}"\n'.format(data['description']))
    fwrite(afterClass,"\s+category =",'    category = "{}"\n'.format(data['category']))
    priority=10
    if 'priority' in data and (data['priority'] == 0 or data['priority']):
        priority=data['priority']
    fwrite(afterClass,"\s+priority =",'    priority = {}\n'.format(priority))
    iconFile=data['icon']
    os.system("mkdir -p {}/icon".format(widgetPath))
    if iconFile and os.path.exists(iconFile):
        os.system("rm {}/icon/* ".format(widgetPath))
        os.system("cp {} {}/icon/".format(iconFile,inputPath))
    else:
        icons=os.listdir(inputPath+'/icon')
        iconFile=os.path.basename(defaultIconFile)
        if not icons:
            os.system("cp {} {}/icon/".format(defaultIconFile,inputPath))
        else:
            iconFile=icons[0]
    finalIconFile = '/widgets/' + widgetName + '/icon/' + os.path.basename(iconFile)
    fwrite(afterClass,"\s+icon =",'    icon = "{}"\n'.format(finalIconFile))
    fwrite(afterClass,"\s+want_main_area =",'    want_main_area = False\n')
    if not 'docker_image_name' in data:
        data['docker_image_name']='biodepot/alpine-bash'
    if not'docker_image_tag' in data:
        data['docker_image_tag'] = 'latest'
    fwrite(afterClass,"\s+docker_image_name =",'    docker_image_name = "{}"\n'.format(data['docker_image_name']))
    fwrite(afterClass,"\s+docker_image_tag =",'    docker_image_tag = "{}"\n'.format(data['docker_image_tag']))
    #inputs and outputs
    if 'inputs' in data and data['inputs']:
        inputStr='['
        for attr, values  in data['inputs'].items():
            if 'callback' in  values and values['callback']:
                inputStr= inputStr+ '("{}",{},"{}"),'.format(attr,deClass(str(values['type'])),values['callback'])
            else:
                inputStr= inputStr+ '("{}",{},"handleInputs{}"),'.format(attr,deClass(str(values['type'])),attr)
        inputStr=inputStr[:-1]+']'
        fwrite(afterClass,"\s+inputs =", '    inputs = {}\n'.format(inputStr))
        
    if 'outputs' in data and data['outputs']:
        outputStr='['
        for attr, value  in data['outputs'].items():
            outputStr= outputStr+ '("{}",{}),'.format(attr,deClass(str(value['type'])))
        outputStr=outputStr[:-1]+']'
        fwrite(afterClass,"\s+outputs =",'    outputs = {}\n'.format(outputStr))
        #permanent settings are not rewritten - we assume that if you changed them you wanted it that way
    if 'parameters' in data and data['parameters']:
        for pname,pvalue in data['parameters'].items():
            searchStr="\s+{}=pset".format(pname)
            print (searchStr)
            if 'default' in pvalue and pvalue['default'] is not None:
                #if it is not a number or dict type then keep quotes
                if pvalue['type'] == 'int' or  pvalue['type'] == 'float' or pvalue['type'] == 'double' or pvalue['type'] == 'bool' or pvalue['type'][-4:] == 'dict' or pvalue['type'][-4:] == 'Dict' or pvalue['type'][-4:] == 'list':
                    fwrite(afterClass,searchStr,'    {}=pset({})\n'.format(pname,pvalue['default']))
                else:
                    fwrite(afterClass,searchStr,'    {}=pset("{}")\n'.format(pname,pvalue['default']))
            else:               
                if pvalue['type'] == type('str') :
                    fwrite(afterClass,searchStr,'    {}=pset("")\n'.format(pname),None)
                if pvalue['type'] == 'bool' :
                    fwrite(afterClass,searchStr,'    {}=pset(False)\n'.format(pname))
                elif pvalue['type'][-4:] == 'list' or pvalue['type'][-4:] == 'List':
                    fwrite(afterClass,searchStr,'    {}=pset([])\n'.format(pname))
                elif pvalue['type'][-4:] == 'dict' or pvalue['type'][-4:] == 'Dict':
                    fwrite(afterClass,searchStr,'    {}=pset({{}})\n'.format(pname))
                else:                   
                    fwrite(afterClass,searchStr,'    {}=pset(None)\n'.format(pname))
    with open(outputWidget,'w') as f:           
        for line in beforeClass:
            f.write(line)
        for line in afterClass:
            f.write(line)
        for line in afterDef:
            f.write(line)
