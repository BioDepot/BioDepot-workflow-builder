#!/usr/bin/env python3
import sys, re, os
from collections import OrderedDict
import jsonpickle
import pprint
def deClass(string):
    #removes the <class 'id'> and returns id
    if string[:6] == '<class':
        m=re.match( r"(\<class ')(\w+)(')",string)
        return m[2]
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

    
def createWidget(inputJson,outputWidget): 
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
        f.write('    icon = "{}"\n'.format(data['icon']))
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
        f.write('        with open("/biodepot/orangebiodepot/{}.json/kallistoQuant.json") as f:\n'.format(os.path.basename(inputJson)))
        f.write('            self.data=f.read()\n')
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
    
def main(args):
    if not sys.argv[1] or not sys.argv[2]:
        sys.stderr.write('Please enter json file and output widget file\n')
        return 1
    sys.stderr.write('From json file {} creating widget file {}\n'.format(sys.argv[1],sys.argv[2]))
    createWidget(sys.argv[1],sys.argv[2])
    return 0
if __name__ == "__main__":
   main(sys.argv[1:])
