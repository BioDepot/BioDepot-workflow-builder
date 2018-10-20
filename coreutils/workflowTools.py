#!/usr/bin/env python3
import sys, os, re 
from makeToolDockCategories import *

def findWidgetsInOWS(owsFile):
    with open(owsFile,'r') as f:
        ows=f.read()
        matches=re.findall(r'qualified_name="(.*?)"',ows)
        return set(matches)

def copyWidgets(widgets,workflowPath,basePath=''):
    #this routine assumes that all widget names are unique which is only true within a workflow or the base widget group
    for widget in widgets:
        #find path to link
        #follow link to find where the actual python script resides
        #copy the parent directory
        widgetParts=widget.split('.')
        pylink='{}/biodepot/{}/{}.py'.format(basePath,widgetParts[0],widgetParts[1])
        target=os.readlink(pylink)
        targetParts=target.split('/')
        widgetName=targetParts[-2]
        destPath='{}/widgets/{}'.format(workflowPath,widgetName)
        if not os.path.exists(destPath):
            sys.stderr.write('cp -r {}/widgets/{} {}\n'.format(basePath,widgetName,destPath))
            os.system('cp -r {}/widgets/{} {}'.format(basePath,widgetName,destPath))


def reformatOWS(workflowName,owsFile):
    sys.stderr.write("workflow is {}\n".format(workflowName))
    ows=""
    with open(owsFile,'r') as f:
        ows=f.read()
        ows=ows.replace('project_name="','project_name="{}'.format(workflowName))
        ows=ows.replace('qualified_name="','qualified_name="{}'.format(workflowName))
    with open(owsFile,'w') as f:
        f.write(ows)
    
def makeWorkflow(workflowName,directory,owsFile,basePath='',copyOWS=False,reformat=False):
    workflowPath= '{}/{}'.format(directory,workflowName)
    os.system('mkdir -p {}'.format(workflowPath))
    os.system('mkdir -p {}/widgets'.format(workflowPath))
    widgets=findWidgetsInOWS(owsFile)
    copyWidgets(widgets,workflowPath,basePath=basePath)
    reformatOWS(workflowName,owsFile)
    os.system('cp {} {}/{}.ows'.format(owsFile,workflowPath,workflowName))

makeWorkflow('dToxSDemo','./','dtoxSWorkflow.ows',basePath='/media/data/home/lhhung/bwb/BioDepot-workflow-builder',reformat=True)
