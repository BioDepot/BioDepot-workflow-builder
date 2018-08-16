#!/usr/bin/env python
import sys, re, os
import glob
#from toolDockEdit.py import addCategoryToToolBox, addWidgetsFromWorkflow
from makeToolDockCategories import *

def addCategoryToToolBox(baseToolPath,category,iconFile=None,background='light-purple'):
    directory=niceForm(category,allowDash=False)
    makeNewDirectory(baseToolPath,directory,iconficonFile,background)
    with open('{}/setup.py'.format(baseToolPath),'a+') as f:
        f.write(entryString(category,directory))
    return directory
    
def addWidgetsFromWorkflow(baseToolPath,workflowPath,destDirectory):
    workflowName=os.basename(workflowPath)
    widgetsDir='{}/widgets'.format(workflowPath)
    widgetNames=os.listdir(widgetsDir)
    for widgetName in widgetNames:
        os.system ("ln -sf  {}/{}/{}.py {}/{}/OW{}.py".format(widgetsDir,widgetName,widgetName,baseToolPath,destDirectory,widgetName))

def loadWorkflow(baseToolPath,workflowPath,background='light-yellow'):
    workflowName=os.basename(workflowPath)
    iconFile=os.listdir('{}/icon'.format(workflowPath))[0]
    directory=addCategoryToToolBox(baseToolPath,workflowName,iconFile='{}/icon/{}'.format(workflowPath,iconFile),background=background)
    addWidgetsFromWorkflow(baseToolPath,workflowPath,directory)    

def main(args):
    workflowPath=args[1]
    background=args[2]
    sys.stderr.write('workflowPath {} background {}\n'.format(workflowPath,background))
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
