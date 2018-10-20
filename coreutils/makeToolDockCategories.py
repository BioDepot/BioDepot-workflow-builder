#!/usr/bin/env python3
import sys, os, re 


defaultCategoryIcon='/icons/misc.png'
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
    ret+='      package_data={{"{}": ["icons/*.svg"]}},\n'.format(directory)
    ret+='      entry_points={{"orange.widgets": "{} = {}"}},)\n'.format(category,directory)
    return ret


def makeNewDirectory(basePath,directory,iconFile,background='light-purple'):
    if os.path.exists('{}/{}'.format(basePath,directory)):
        return
    os.system('mkdir -p {}/{}/icon'.format(basePath,directory))
    os.system('touch {}/{}/__init__.py'.format(basePath,directory))
    with open('{}/{}/__init__.py'.format(basePath,directory),'w') as f:
        f.write('import sysconfig\n')
        if not iconFile or not os.path.exists (iconFile):
            iconFile = defaultCategoryIcon
        iconName=os.path.basename(iconFile)
        os.system('cp {} {}/{}/icon/.'.format(iconFile,basePath,directory))
        f.write('ICON = "icon/{}"\n'.format(iconName))
        f.write('BACKGROUND ="{}"\n'.format(background))

def removeCategoryFromToolDock(basePath,category,directory):
    setupFile='{}/setup.py'.format(basePath)
    os.system(''' sed -i '/setup(name="{}"/,/)/d' {}'''.format(category,setupFile))
    os.system('cd {} && rm -rf {} && rm -rf {}.egg'.format(basePath,directory,directory))




    
    




