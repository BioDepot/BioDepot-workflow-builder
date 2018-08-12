#!/usr/bin/env python3
import sys, os, re

def niceForm(badString,allowDash='True'):
    #dashes cause problems with hasattr
    #strip and sub spaces
    ret=re.sub("\s+","_",badString.strip())
    #strip other bad characters
    ret=re.sub('[^a-zA-Z0-9\_\-\.]', '', ret)
    #strip 
    if not allowDash:
        ret=re.sub('-', '_', ret)
    return ret
        
def entryString(category,directory):
    ret='setup(name="{}",\n'.format(category)
    ret+='      packages=["{}"],\n'.format(directory)
    ret+='      package_data={"{}": ["icons/*.svg"]},\n'.format(directory)
    ret+='      entry_points={"orange.widgets": "{} = {}"},)\n'.format(category)
    return ret
    
def initPyString(icon,background):
    pass
    
def makeNewDirectory(basePath,directory,iconFile):
    if os.path.exists('{}/{}'.format(basePath,directory)):
        return
    os.system('mkdir -p {}/{}'.format(basePath,directory))
    os.system('mkdir -p {}/{}/icon'.format(basePath,directory))
    os.system('ln -s {}/{}/icons ../orangebiodepot/icons'.format(basePath,directory))
    os.system('touch {}/{}/__init__.py'.format(basePath,directory))
    with open('{}/{}/__init__.py'.format(basePath,directory),'w') as f:
        fwrite('import sysconfig\n')
        if iconFile and os.path.exists (iconFile):
            os.system(cp {} {}/
        else:
        
def writeSetup(categories,setupFile):
    directories=[]
    #make backup if it exists
    if os.path.exists(setupFile):
        os.system('cp {} {}.bak'.setupFile)
    with open(setupFile,'w') as f:
        for index,category in enumerate(categories):
            category=niceForm(category)
            directory=niceForm(category,allowDash=False)
            f.write(entryString(category,directory))
            directories.append(directory)




    
    




