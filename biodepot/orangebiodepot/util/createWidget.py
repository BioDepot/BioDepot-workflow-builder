#!/usr/bin/env python3
from subprocess import call
from shutil import copyfile
import sys, re, os, getopt
from collections import OrderedDict
import jsonpickle
import pprint
def deClass(string):
    #removes the <class 'id'> and returns id
    if string[:6] == '<class':
        m=re.match( r"(\<class ')(\w+)(')",string)
        try:
            return m.group(2)
        except:
            return string
    return string
    

WIDGET_HEADING ='''import os
import glob
import sys
import functools
import jsonpickle
from collections import OrderedDict
from Orange.widgets import widget, gui, settings
from orangebiodepot.util.DockerClient import DockerClient
from orangebiodepot.util.BwBase import OWBwBWidget, ConnectionDict, ContainerPaths, BwbGuiElements
from PyQt5 import QtWidgets, QtGui

'''

    
def createWidget(inputJson,outputWidget, registerFlag=False): 
    data={}
    with open(inputJson) as f:
        data=jsonpickle.decode(f.read())
    f.close()

            
    #write preInit
    with open(outputWidget,'w') as f:
        f.write(WIDGET_HEADING)
        className='class OW{}(OWBwBWidget):\n'.format(data['name'].replace(' ',''))
        f.write(className)
        f.write('    name = "{}"\n'.format(data['name']))
        f.write('    description = "{}"\n'.format(data['description']))
        f.write('    category = "{}"\n'.format(data['category']))
        f.write('    priority = 10\n')
        iconFile=data['icon']
        if os.path.dirname(iconFile) != '/biodepot/orangebiodepot/icons':
            iconFile = '/biodepot/orangebiodepot/icons/' + os.path.basename(iconFile)
        f.write('    icon = "{}"\n'.format(iconFile))
        f.write('    want_main_area = False\n')
        f.write('    docker_image_name = "{}"\n'.format(data['docker_image_name']))
        f.write('    docker_image_tag = "{}"\n'.format(data['docker_image_tag']))
        #inputs and outputs
        inputStr='['
        for attr, value  in data['inputs'].items():
                inputStr= inputStr+ '("{}",{},"handleInputs{}"),'.format(attr,deClass(str(value['type'])),attr)
        inputStr=inputStr[:-1]+']'
        f.write('    inputs = {}\n'.format(inputStr))
        outputStr='['
        for attr, value  in data['outputs'].items():
                outputStr= outputStr+ '("{}",{}),'.format(attr,deClass(str(value['type'])))
        outputStr=outputStr[:-1]+']'
        f.write('    outputs = {}\n'.format(outputStr))
        #permanent settings
        f.write('    pset=functools.partial(settings.Setting,schema_only=True)\n')
        f.write('    runMode=pset(0)\n')
        f.write('    runTriggers=pset([])\n')
        f.write('    inputConnectionsStore=pset({})\n')    
        for pname,pvalue in data['parameters'].items():
            if 'default' in pvalue and pvalue['default'] is not None:
                f.write('    {}=pset({})\n'.format(pname,str(pvalue['default'])))
            else:
                f.write('    {}=pset(None)\n'.format(pname))
        f.write('    def __init__(self):\n')        
        f.write('        super().__init__(self.docker_image_name, self.docker_image_tag)\n')
        joinedName=data['name'].replace(' ','')
        f.write('        with open("/biodepot/orangebiodepot/json/{}") as f:\n'.format(os.path.basename(inputJson)))
        f.write('            self.data=jsonpickle.decode(f.read())\n')
        f.write('            f.close()\n')
        #init
        f.write('        self.initVolumes()\n')
        f.write('        self.inputConnections = ConnectionDict(self.inputConnectionsStore)\n')
        f.write('        self.drawGUI()\n')
        #input callbacks
        for attr, value in data['inputs'].items():
            f.write('    def handleInputs{}(self, value, sourceId=None):\n'.format(attr))
            f.write('        self.handleInputs(value, "{}", sourceId=None)\n'.format(attr))
        #output callback
        f.write('    def handleOutputs(self):\n'.format(attr))
        for attr in data['outputs']:
            if 'outputValue' in data['outputs']:
                outputValue=data['outputs']['outputValue']
                f.write('        self.send("{}", {}\n'.format(attr,outputValue))
            else:
                f.write('        outputValue=None\n')
                f.write('        if hasattr(self,"{}"):\n'.format(attr))
                f.write('            outputValue=getattr(self,"{}")\n'.format(attr))
                f.write('        self.send("{}", outputValue)\n'.format(attr))
        if registerFlag:
            register(inputJson,outputWidget,data['icon'])

def register(jsonFile,widgetFile,iconFile):
    sys.stderr.write('moving files to correct biodepot locations\n')
    if os.path.dirname(jsonFile) != '/biodepot/orangebiodepot/json':
        newjsonFile = '/biodepot/orangebiodepot/json/' + os.path.basename(jsonFile)
        copyfile(jsonFile, newjsonFile)
    if os.path.dirname(iconFile) != '/biodepot/orangebiodepot/icons':
        newiconFile = '/biodepot/orangebiodepot/icons/' + os.path.basename(iconFile)
        copyfile(iconFile, newiconFile)   
    if os.path.dirname(widgetFile) != '/biodepot/orangebiodepot':
        basename=os.path.basename(widgetFile)
        if basename[:2] != 'OW':
            basename='OW'+basename
        newwidgetFile = '/biodepot/orangebiodepot/' + basename
        copyfile(widgetFile, newwidgetFile)

def usage():
    sys.stderr.write('createWidgetTest -i <parms.json> -o <widget.py>\n')
    sys.stderr.write('add -r or --register flag to save widget, json and icon files to correct places in biodepot container\n')    

def main(args):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hri:o:", ["help", "output=" ,"input=","register"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    outputWidget = None
    verbose = False
    inputJson= None
    register=False
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-o", "--output"):
            outputWidget = a
        elif o in ("-i", "--input"):
            inputJson = a            
        elif o in ("-r", "--register"):
            register = True   
        else:
            assert False, "unhandled option"
    if not outputWidget and not inputJson:
        usage()
        sys.exit(2)
    sys.stderr.write('From json file {} creating widget file {}\n'.format(inputJson,outputWidget))
    createWidget(inputJson,outputWidget,registerFlag=register)
    return 0
if __name__ == "__main__":
   main(sys.argv[1:])
