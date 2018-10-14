import sys, os, re, shutil, glob, tempfile, jsonpickle
from glob import glob 
from xml.dom import minidom
from collections import Counter, OrderedDict
from itertools import groupby
import datetime
        
def findMode(lst):
    maxFreq =  max(map(lst.count, lst))
    modes = [i for i in lst if lst.count(i) == maxFreq]
    return modes
    
def reformatOWS(workflowTitle,inputFile,outputFile,uniqueNames):
    doc = minidom.parse(inputFile)
    scheme = doc.getElementsByTagName("scheme")[0]
    scheme.attributes['title'].value=workflowTitle
    nodes = doc.getElementsByTagName("node")
    for node in nodes:
        projectName=node.getAttribute('project_name')
        widgetName=node.getAttribute('name')
        if uniqueNames and projectName in uniqueNames and widgetName in uniqueNames[projectName]:
            uniqueName=uniqueNames[projectName][widgetName]
            
            #link will be at biodepot/workflowPath/OWuniqueName -> uniqueName/widgetName.py
            #class is still OWwidgetName
        else:
            uniqueName=widgetName
        node.attributes['name'].value=uniqueName
        node.attributes['qualified_name']=workflowName+'.OW'+uniqueName+'.OW'+uniqueName
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
        

def replaceName(jsonFile,key,newValue):
    with open(jsonFile) as f:
        data=jsonpickle.decode(f.read())
    data['name']=newValue
    myJdata=jsonpickle.encode(data)
    with open(jsonFile,"w") as f:
        f.write(myJdata)
        
def replaceNamePy(pyFile,oldName,newName):
    with open(pyFile) as f:
        data=f.read()
    with open(pyFile,'w') as f:
        for line in data.splitlines(True):
            #look for
            #class OWFile(OWBwBWidget):
            if line[0:8] == 'class OW':
                line='class OW{}(OWBwBWidget):\n'.format(newName)
            elif line.strip() == 'name = "{}"'.format(oldName):
                line='    name = "{}"\n'.format(newName)
            elif line.strip() == 'with open(getJsonName(__file__,"{}")) as f:'.format(oldName):
                line='        with open(getJsonName(__file__,"{}")) as f:\n'.format(newName)
            f.write(line)
    
def renameWidgetInToolDock(oldPyPath,newPyPath):
    directoryList=(str(os.popen('''grep -oP 'packages=\["\K[^"]+' /biodepot/setup.py''').read())).split()
    widgetPy=os.path.basename(newPyPath)
    for directory in directoryList:
        pyLinks=glob('/biodepot/{}/OW*.py'.format(directory))
        for pyLink in pyLinks:
            if os.path.realpath(pyLink) == os.path.realpath(oldPyPath):        
                os.unlink(pyLink)
                sys.stderr.write('ln -s {} /biodepot/{}/OW{}'.format(newPyPath,directory,widgetPy)) 
                os.system('ln -s {} /biodepot/{}/OW{}'.format(newPyPath,directory,widgetPy)) 
                
def renameWidget(srcWidget,oldName,newName):
    if oldName == newName:
        return
    #first rename parent if possible
    newPath=os.path.dirname(srcWidget)+'/'+newName
    os.system('mv {} {}'.format(srcWidget,newPath))

    #now we need to rename each of the files
    #and after renaming replace the instances of oldName in the data structures and python script with newName
    oldFiles=glob('{}/{}.*'.format(newPath,oldName))
    for oldFile in oldFiles:
        basename,extension=os.path.splitext(oldFile)
        newFile=newPath+'/'+newName+extension
        os.rename(oldFile,newFile)
        if extension == '.py':
            replaceNamePy(newFile,oldName,newName)
        else:
            replaceName(newFile,'name',newName)
    oldPyPath='{}/{}.py'.format(srcWidget,oldName)
    newPyPath='{}/{}.py'.format(newPath,newName)
    renameWidgetInToolDock(oldPyPath,newPyPath)
    

        
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
    if widgetNameSeen(widgetName,projectNames):
        parts=widgetName.split('_')
        if parts and parts[-1].isdigit():
            number=int(parts[-1])+1
        else:
            number=1
        return '{}_{}'.format(widgetName,number)
    else:
        return projectName +'_'+widgetName

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

def renameWidgetInWorkflow(inputOWS,outputOWS,projectName,widgetName,newName):
    #remove nodes that match projectName and widgetName and record them
    doc = minidom.parse(inputOWS)
    nodeParent=doc.getElementsByTagName("nodes")[0]
    nodes = doc.getElementsByTagName("node")
    for node in nodes:
        print (node.getAttribute('project_name'))
        print (node.getAttribute('name'))
        if node.getAttribute('project_name')==projectName and node.getAttribute('name') == widgetName:
            node.attributes['name']=newName
            node.attributes['qualified_name'].value='{}.OW{}.OW{}'.format(projectName,newName,newName)
#get rid of empty with one space lines that result from removal and write    
    with open(outputOWS,'w') as f:
        f.write("".join([s for s in doc.toxml().splitlines(True) if s.strip()]))
    
def exportWorkflow (bwbOWS,outputWorkflow,projectTitle,merge=False,color=None,iconFile=None,basePath=""):
    tempDir = tempfile.mkdtemp()
    os.mkdir(tempDir+'/widgets')
    os.mkdir(tempDir+'/widgets/{}'.format(projectTitle))
    tempOWS='{}/{}'.format(tempDir,os.path.basename(bwbOWS))
    shutil.copyfile(bwbOWS,tempOWS)
    doc = minidom.parse(bwbOWS)
    nodes = doc.getElementsByTagName("node")
    if not nodes:
        shutil.rmtree(tempDir)
        return False
    #find the old title
    myScheme = doc.getElementsByTagName("scheme")[0]
    oldProjectTitle = myScheme.getAttribute('title')
    
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
            widgetName=node.getAttribute('name')
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
                uniquePath='{}/widgets/{}/{}'.format(tempDir,projectTitle,uniqueName)
                shutil.copytree(widgetPath,uniquePath)
                if nameChanged:
                    renameWidget(uniquePath,widgetName,uniqueName)
        
        reformatOWS(projectTitle,bwbOWS,tempOWS,uniqueNames)
        
    else:
        widgetPaths=set()
        for node in nodes:
            projectName=node.getAttribute('project_name')
            if projectName==projectTitle:
                widgetName=node.getAttribute('name')
                qname=node.getAttribute('qualified_name')
                widgetPath=findWidgetPathFromLink(qname,projectPath,basePath)
                if widgetPath not in widgetPaths:
                    widgetPaths.add(widgetPath)
                    os.system('cp -r {} {}/widgets/{}/{}'.format(widgetPath,tempDir,projectTitle,widgetName)) 
          
    #check if the original exists for icon and __init__.py - otherwise take them from the title otherwise take it from the 
    if os.path.exists('{}/widgets/{}/icon'.format(outputWorkflow,projectTitle)) and os.path.exists('{}/widgets/{}/__init__.py'.format(outputWorkflow,projectTitle)):
        iconPath='{}/widgets/{}/icon'.format(outputWorkflow,projectTitle)
        initPath='{}/widgets/{}/__init__.py'.format(outputWorkflow,projectTitle)
    elif os.path.exists('/biodepot/{}/icon'.format(oldProjectTitle)):
        iconPath='/biodepot/{}/icon'.format(oldProjectTitle)
        initPath='/biodepot/{}/__init__.py'.format(oldProjectTitle)
    elif os.path.exists('/biodepot/{}/icon'.format(projectTitle)):
        iconPath='/biodepot/{}/icon'.format(projectTitle)
        initPath='/biodepot/{}/__init__.py'.format(projectTitle)
    else:
        iconPath='/biodepot/User/icon'
        initPath='/biodepot/User/__init__.py'
        if not color:
            color='light-blue'
    shutil.copytree(iconPath,'{}/widgets/{}/icon'.format(tempDir,projectTitle))
    shutil.copyfile(initPath, '{}/widgets/{}/__init__.py'.format(tempDir,projectTitle))
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
        iconDir="{}/widgets/{}/icon".format(tempDir,projectTitle)
        os.system('rm -f {}/* && cp {} {}/. '.format(iconDir,iconFile,iconDir))
        changedIcon=True

    if changedIcon or color:
        initFile='{}/widgets/{}/__init__.py'.format(tempDir,projectTitle)
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
    projectTitle=os.path.basename(workflowDir)
    os.system('mkdir -p /biodepot/{}'.format(projectTitle))
    with open('/biodepot/setup.py','r') as f:
        setupData=f.read()
    projectList=re.findall(r'setup\(name="([^"]+)"',setupData)
    changedSetup=True
    for projectName in projectList:
        if projectName == projectTitle:
            changedSetup=False
            break
    if changedSetup:
        setupData+=entryString(projectTitle,projectTitle)
    pythonFiles=glob('{}/widgets/*/*/*.py'.format(workflowDir))
    for pythonFile in pythonFiles:
        basePythonFile=os.path.basename(pythonFile)
        destLink='/biodepot/{}/OW{}'.format(projectTitle,basePythonFile)
        os.system('ln -sf {} {}'.format(pythonFile,destLink))
    #make link to icons and __init__.py
    print ('ln -sf {}/widgets/{}/icon /biodepot/{}/icon'.format(workflowDir,projectTitle,projectTitle))
    if os.path.exists('/biodepot/{}/icon'.format(projectTitle)):
        os.unlink('/biodepot/{}/icon'.format(projectTitle))
    os.system('ln -sf {}/widgets/{}/icon /biodepot/{}/icon'.format(workflowDir,projectTitle,projectTitle))
    os.system('ln -sf {}/widgets/{}/__init__.py /biodepot/{}/__init__.py'.format(workflowDir,projectTitle,projectTitle))
    if changedSetup:
        with open('/biodepot/setup.py','w') as f:
            f.write(setupData)
        os.system('cd /biodepot && pip install -e .')

