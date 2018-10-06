import sys, os, re, shutil, glob, tempfile 
from makeToolDockCategories import *
from xml.dom import minidom
from collections import Counter, OrderedDict
from itertools import groupby
import datetime
        
def findMode(lst):
    maxFreq =  max(map(lst.count, lst))
    modes = [i for i in lst if lst.count(i) == maxFreq]
    return modes
    
def reformatOWS(workflowTitle,inputFile,outputFile,uniqueNames):
    workflowPath=niceForm(workflowTitle,allowDash=False)
    doc = minidom.parse(inputFile)
    scheme = doc.getElementsByTagName("scheme")[0]
    scheme.attributes['title'].value=workflowTitle
    nodes = doc.getElementsByTagName("node")
    for node in nodes:
        projectName=node.getAttribute('project_name')
        projectPath=niceForm(projectName,allowDash=False)
        widgetName=niceForm(node.getAttribute('name'),allowDash=False)
        if uniqueNames and projectName in uniqueNames and widgetName in uniqueNames[projectName]:
            uniqueName=uniqueNames[projectName][widgetName]
            
            #link will be at biodepot/workflowPath/OWuniqueName -> uniqueName/widgetName.py
            #class is still OWwidgetName
        else:
            uniqueName=widgetName
        node.attributes['name'].value=uniqueName
        node.attributes['qualified_name']=workflowPath+'.OW'+uniqueName+'.OW'+widgetName
        node.attributes['project_name'].value=workflowTitle
    with open(outputFile,'w') as f:
        f.write(doc.toxml())
            
def changeNameInOWS(oldName,newName,filename):
    doc = minidom.parse(filename)
    nodes = doc.getElementsByTagName("node")
    for node in nodes:
        node.attributes['project_name'].value=newName
        qnameParts=node.getAttribute('qualified_name').split('.')
        if qnameParts[0] == oldName:
            qnameParts[0]=newName;
            node.attributes['qualified_name'].value='.'.join(qnameParts)
        node.attributes['project_name'].value=newName
     
def copyWorkflow(inputWorkflow,outputWorkflow):
    #to copy workflow copy all the files and then change the titles in the ows file
    #could use global regex but then we could run into troubles when the workflow names are part of other names
    os.system('cp -r inputWorkflow outputWorkflow')
    oldName=os.path.basename(os.path.normpath(inputWorkflow))
    newName=os.path.basename(os.path.normpath(outputWorkflow))
    if newName != oldName :
        oldOWS=outputWorkflow+'/'+oldName+'.ows'
        newOWS=outputWorkflow+'/'+newName+'.ows'
        os.system('mv oldOWS newOWS')
        changeNameInOWS(oldName,newName,newOWS)
        
def findWidgetPathFromLink(qualifiedName,groupName,basePath=''):
    parts=qualifiedName.split('.')
    link=basePath+'/biodepot/'+'/'.join(parts[0:-1])+'.py'
    widgetPath=os.path.dirname((os.readlink(link)))
   #check if absolute path
    if widgetPath[0]=='/':
        absWidgetPath=widgetPath
    else:
        absWidgetPath=os.path.abspath('{}/biodepot/{}/{}'.format(basePath,groupName,widgetPath))
    return absWidgetPath

def widgetNameSeen(name,projectPaths):
    #see if the first part of the name is in one of the projectPaths
    namelen=len(name)
    for projectPath in projectPaths:
        projectlen=len(projectPath)
        if projectlen <= namelen and name[0:projectlen] == projectPath:
            return True
    return False
    
def modifyWidgetName(widgetName,projectName, projectNames):
    #tries to get a unique name when there
    #first try to add prefix if there is no prefix
    #we check against all project names to keep multiple chaining of names to a minimum
    projectPaths=[]
    for projectName in projectNames:
        projectPaths.append(niceForm(projectName,allowDash=False))
    if widgetNameSeen(widgetName,projectPaths):
        parts=widgetName.split('_')
        if parts and parts[-1].isdigit():
            number=int(parts[-1])+1
        else:
            number=1
        return '{}_{}'.format(widgetName,number)
    else:
        return niceForm(projectName,allowDash=False)+'_'+widgetName

def removeWidgetfromWorkflow(inputOWS,outputOWS,projectName,widgetName):
    #remove nodes that match projectName and widgetName and record them
    removedNodes=set()
    doc = minidom.parse(inputOWS)
    nodeParent=doc.getElementsByTagName("nodes")[0]
    nodes = doc.getElementsByTagName("node")
    for node in nodes:
        if node.getAttribute('project_name')==projectName and node.getAttribute('name') == widgetName:
            removedNodes.add(node.getAttribute('id'))
            nodeParent.removeChild(node)
    linkParent=doc.getElementsByTagName("links")[0]
    links=doc.getElementsByTagName("link")
    for link in links:
        if link.getAttribute('source_node_id') in removedNodes or link.getAttribute('sink_node_id') in removedNodes:
            linkParent.removeChild(link)
    propertiesParent=doc.getElementsByTagName("node_properties")[0]
    properties=doc.getElementsByTagName("properties")
    for nodeProperty in properties:
        if nodeProperty.getAttribute('node_id') in removedNodes:
            propertiesParent.removeChild(nodeProperty)
#get rid of empty with one space lines that result from removal and write    
    with open(outputOWS,'w') as f:
        f.write("".join([s for s in doc.toxml().splitlines(True) if s.strip()]))

    
def exportWorkflow (bwbOWS,outputWorkflow,projectTitle,merge=False,color=None,iconFile=None,basePath=""):
    tempDir = tempfile.mkdtemp()
    os.mkdir(tempDir+'/widgets')
    projectTitlePath=niceForm(projectTitle,allowDash=False)
    os.mkdir(tempDir+'/widgets/{}'.format(projectTitlePath))
    tempOWS='{}/{}'.format(tempDir,os.path.basename(bwbOWS))
    shutil.copyfile(bwbOWS,tempOWS)
    doc = minidom.parse(bwbOWS)
    nodes = doc.getElementsByTagName("node")
    if not nodes:
        shutil.rmtree(tempDir)
        return False 
    #get rid of anything like \n in the color string    
    if color:
        color=color.strip()
    if merge:
        projectNames=[]
        widgetPaths={}
        #gather paths and order the projects names by ProjectName and then other projects by number of occurences
        #this is the order to copy/convert the names
        for node in nodes:
            projectName=node.getAttribute('project_name')
            projectPath=niceForm(projectName,allowDash=False)
            widgetName=niceForm(node.getAttribute('name'),allowDash=False)
            qname=node.getAttribute('qualified_name')
            widgetPath=findWidgetPathFromLink(qname,projectPath,basePath)
        #keep unique list of widgetPaths for each project
            if projectName not in widgetPaths:
                widgetPaths[projectName]=[]
            if widgetPath not in widgetPaths[projectName]:
                widgetPaths[projectName].append(widgetPath)
        
        #keep track of projectNames - do not unique - want to count their occurrence and sort    
#            if projectTitle != projectName:
            projectNames.append(projectName)
            
        if projectNames:
            #sort by frequency of occurence
            projectNames.sort(key=Counter(projectNames).get, reverse=True)
            #unique the list
            projectNames=list(OrderedDict.fromkeys(projectNames))
    
#        if projectName in widgetPaths:
#            projectNames=[projectName]+projectNames
            
        #copy the widgets and make a map for renamed widgets
        nameSeen=set()
        uniqueNames={}
        for projectName in projectNames:
            for widgetPath in widgetPaths[projectName]:
                widgetName=os.path.basename(widgetPath)
                uniqueName=os.path.basename(widgetPath)
                nameChanged=False
                while uniqueName in nameSeen:
                    uniqueName=modifyWidgetName(uniqueName,projectName,projectNames)
                    nameChanged=True
                if nameChanged:
                    if projectName not in uniqueNames:
                        uniqueNames[projectName]={}
                    uniqueNames[projectName][widgetName]=uniqueName
                nameSeen.add(uniqueName)
                os.system('cp -r {} {}/widgets/{}/{}'.format(widgetPath,tempDir,projectTitlePath,uniqueName))                
        reformatOWS(projectTitle,bwbOWS,tempOWS,uniqueNames)
    else:
        widgetPaths=set()
        for node in nodes:
            projectName=node.getAttribute('project_name')
            projectPath=niceForm(projectName,allowDash=False)
            if projectPath==projectTitlePath:
                widgetName=niceForm(node.getAttribute('name'),allowDash=False)
                qname=node.getAttribute('qualified_name')
                widgetPath=findWidgetPathFromLink(qname,projectPath,basePath)
                if widgetPath not in widgetPaths:
                    widgetPaths.add(widgetPath)
                    os.system('cp -r {} {}/widgets/{}/{}'.format(widgetPath,tempDir,projectTitlePath,widgetName)) 
          
    #check if the original exists for icon and __init__.py - otherwise take them from the User drawer 
    if os.path.exists('{}/widgets/{}/icon'.format(outputWorkflow,projectTitlePath)) and os.path.exists('{}/widgets/{}/__init__.py'.format(outputWorkflow,projectTitlePath)):
        src='{}/widgets/{}/icon'.format(outputWorkflow,projectTitlePath)
        dst='{}/widgets/{}/icon'.format(tempDir,projectTitlePath)
        shutil.copytree(src, dst, symlinks=True)
        srcFile='{}/widgets/{}/__init__.py'.format(outputWorkflow,projectTitlePath)
        dstFile='{}/widgets/{}/__init__.py'.format(tempDir,projectTitlePath)
        shutil.copyfile(srcFile, dstFile)
    else:
        src='/biodepot/User/icon'
        dst='{}/widgets/{}/icon'.format(tempDir,projectTitlePath)
        shutil.copytree(src, dst,symlinks=True)
        srcFile='/biodepot/User/__init__.py'
        dstFile='{}/widgets/{}/__init__.py'.format(tempDir,projectTitlePath)
        shutil.copyfile(srcFile, dstFile)
        if not color:
            color='light-blue'
    
    changedIcon=False
    
    #check that iconFile exists and resolve if symlink
    if iconFile:
        if not os.path.exists(iconFile):
            iconFile=""
        else:
            iconFile=os.path.realpath(iconFile)
            if not os.path.isfile(iconFile):
                iconFile=""
                
    #if there is an icon
    if iconFile:
        iconDir="{}/widgets/{}/icon".format(tempDir,projectTitlePath)
        os.system('rm -f {}/* && cp {} {}/. '.format(iconDir,iconFile,iconDir))
        changedIcon=True

    if changedIcon or color:
        initFile='{}/widgets/{}/__init__.py'.format(tempDir,projectTitlePath)
        with open(initFile,'r') as f:
            lines = f.readlines()
        with open(initFile,'w') as f:
            for line in lines:
                if changedIcon and line[0:4] == 'ICON':
                    line='ICON = "icon/{}"\n'.format(os.path.basename(iconFile))
                elif color and line[0:10] == 'BACKGROUND':
                    line='BACKGROUND = "{}"\n'.format(color)
                f.write(line)
    if outputWorkflow and os.path.normpath(outputWorkflow) != '/' and os.path.isdir(outputWorkflow):
        shutil.rmtree(outputWorkflow)
        try:
            os.system('mv {} {}'.format(tempDir,outputWorkflow))
        except Exception as e:
            shutil.rmtree(tempDir)
    else:
        shutil.rmtree(tempDir)

def importWorkflow(owsFile):
    changedSetup=False
    workflowDir=os.path.dirname(owsFile)
    projectTitlePath=os.path.basename(workflowDir)
    os.system('mkdir -p /biodepot/{}'.format(projectTitlePath))

    with open('/biodepot/setup.py','r') as f:
        setupData=f.read()
    projectList=re.findall(r'setup\(name="([^"]+)"',setupData)
    changedSetup=True
    for projectName in projectList:
        projectPath=niceForm(projectName,allowDash=False)
        if projectPath == projectTitlePath:
            changedSetup=False
            break
    if changedSetup:
        setupData+=entryString(projectTitlePath,projectTitlePath)
        
    doc = minidom.parse(owsFile)
    nodes = doc.getElementsByTagName("node")        
    if not nodes:
        return
    #make links to widgets
    for node in nodes:
        #differs from export in that we want to preserve the dashes in the names for changing setup.py later
        projectPath=niceForm(node.getAttribute('project_name'),allowDash=False)
        #only care about the ones specific to this workflow
        if projectPath == projectTitlePath:
            widgetName=niceForm(node.getAttribute('name'),allowDash=False)
            qname=node.getAttribute('qualified_name')
            parts=qname.split('.')
            destLink='/biodepot/'+'/'.join(parts[0:-1])+'.py'
            pythonFile='{}/widgets/{}/{}/{}.py'.format(workflowDir,projectPath,widgetName,parts[-1][2:])
            print ('ln -sf {} {}'.format(pythonFile,destLink))
            os.system('ln -sf {} {}'.format(pythonFile,destLink))
    #make link to icons and __init__.py
    print ('ln -sf {}/widgets/{}/icon /biodepot/{}/icon'.format(workflowDir,projectTitlePath,projectTitlePath))
    if os.path.exists('/biodepot/{}/icon'.format(projectTitlePath)):
        os.unlink('/biodepot/{}/icon'.format(projectTitlePath))
    os.system('ln -sf {}/widgets/{}/icon /biodepot/{}/icon'.format(workflowDir,projectTitlePath,projectTitlePath))
    os.system('ln -sf {}/widgets/{}/__init__.py /biodepot/{}/__init__.py'.format(workflowDir,projectTitlePath,projectTitlePath))
    if changedSetup:
        with open('/biodepot/setup.py','w') as f:
            f.write(setupData)
        os.system('cd /biodepot && pip install -e .')

