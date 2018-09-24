#!/usr/bin/env python3
import sys, os, re, subprocess
from xml.dom import minidom
from itertools import chain
    
def copyWidgets(basePath='/media/data/home/lhhung/bwb/BioDepot-workflow-builder/biodepot/'):
    dirNames=subprocess.check_output(['ls', 'newWidgets']).decode('utf-8').strip().split('\n')
    for dirName in dirNames:
        pylinks=subprocess.check_output(['find', 'biodepot/'+dirName,'-type', 'l','-name','*.py']).decode('utf-8').strip().split('\n')
        for pylink in pylinks:
            widgetPath=os.path.dirname((os.readlink(pylink)))
            #check if absolute path
            if widgetPath[0]=='/':
                absWidgetPath='/media/data/home/lhhung/bwb/BioDepot-workflow-builder/'+widgetPath
            else:
                absWidgetPath=os.path.abspath('/media/data/home/lhhung/bwb/BioDepot-workflow-builder/biodepot/{}/{}'.format(dirName,widgetPath))
            #sys.stderr.write('widgetPath is {} abs path is {}\n'.format(widgetPath,absWidgetPath))
            widgetName=os.path.basename(widgetPath)
            cmd='cp -r {} newWidgets/{}'.format(absWidgetPath,dirName)
            print (cmd)
            os.system(cmd)

def makeLinks (basePath='/media/data/home/lhhung/bwb/BioDepot-workflow-builder/biodepot/'):
    dirNames=subprocess.check_output(['ls', 'widgets']).decode('utf-8').strip().split('\n')
    for dirName in dirNames:
        if dirName[0:4] == 'Demo':
            continue
        pylinks=subprocess.check_output(['find', 'biodepot/'+dirName,'-type', 'l','-name','*.py']).decode('utf-8').strip().split('\n')
        for pylink in pylinks:
            widgetPath=os.readlink(pylink)      
            if widgetPath[0]=='/':
                widgetPath='../..' + widgetPath
            pathDirs=list(filter(None, widgetPath.split('/')))
            pathDirs.insert(-2,dirName)
            linkName='/'.join(pathDirs)
            cmd='ln -fs {} {}'.format(linkName,pylink)
            print(cmd)
            os.system(cmd)

makeLinks()
