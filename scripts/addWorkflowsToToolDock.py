#!/usr/bin/env python3
import sys, os, re 
from glob import glob
def niceForm(badString, useDash="False"):
    # dashes cause problems with hasattr
    # strip and sub spaces
    ret = re.sub("\s+", "_", badString.strip())
    # strip other bad characters
    ret = re.sub("[^a-zA-Z0-9\_\-\.]", "", ret)
    # strip
    if useDash:
        ret = re.sub("_", "-", ret)
    else:
        ret = re.sub("-", "_", ret)
    return ret

def entryString(category,directory):
    ret='setup(\n'
    ret+='      name="{}",\n'.format(category)
    ret+='      packages=["{}"],\n'.format(directory)
    ret+='      package_data={{"{}": ["icons/*.svg"]}},\n'.format(directory)
    ret+='      entry_points={{"orange.widgets": "{} = {}"}},\n)\n'.format(category,directory)
    return ret
    
def importWorkflow(workflowDir,basePath="/biodepot"):
    projectTitlePath = os.path.basename(workflowDir)
    projectTitle = niceForm(projectTitlePath, useDash=True)
    with open("{}/setup.py".format(basePath), "r") as f:
        setupData = f.read()
    projectList = re.findall(r'setup\(\s+name="([^"]+)"', setupData)
    if projectTitle in projectList:
        return
    symDir = "{}/{}".format(basePath,projectTitlePath)
    if not os.path.exists(symDir):
        os.makedirs(symDir)
    setupData += entryString(projectTitle, projectTitlePath)
    #copy widgets
    os.system("cp -r {}/widgets/{} {}/../widgets/.".format(workflowDir,projectTitlePath,basePath))
    output=os.popen("cd {} && ls ../../widgets/{}/*/*.py".format(symDir,projectTitlePath))
    os.system("mv {}/../widgets/{}/icon {}/icon".format(basePath,projectTitlePath,symDir))
    os.system("mv {}/../widgets/{}/__init__.py {}/__init__.py".format(basePath,projectTitlePath,symDir))
    pythonFiles=output.read().splitlines()
    for pythonFile in pythonFiles:
        basePythonFile = os.path.basename(pythonFile)
        destLink = "{}/OW{}".format(symDir,basePythonFile)
        os.system ("ln -sf {} {}".format(pythonFile, destLink))
    with open("{}/setup.py".format(basePath), "w") as f:
        f.write(setupData)
            
def importWorkflows(owsFileList,basePath):
    lines=[line.rstrip('\n') for line in open(owsFileList)]
    if lines:
        for line in lines:
            importWorkflow(line,basePath)
    os.system("cd {} && pip install -e .".format(basePath))
        
if __name__ == "__main__":
    if len(sys.argv) == 2:
        basePath="/biodepot"
        importWorkflows(sys.argv[1],basePath)
