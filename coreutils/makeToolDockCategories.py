#!/usr/bin/env python3
import sys, os, re

# pip insists on changing underscores to dashes
# however import will not work with python package names with dashes - this includes the paths
# so until we get rid of pip for dynamic installation we have to support two forms of the name - a dash form for package names and an underscore form for widget names and their paths which is derived from the package name

defaultCategoryIcon = "/icons/misc.png"


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


def entryString(category, directory):
    ret = 'setup(name="{}",\n'.format(category)
    ret += '      packages=["{}"],\n'.format(directory)
    ret += '      package_data={{"{}": ["icons/*.svg"]}},\n'.format(directory)
    ret += '      entry_points={{"orange.widgets": "{} = {}"}},)\n'.format(
        category, directory
    )
    return ret


def makeNewDirectory(basePath, directory, iconFile, background="light-purple"):
    if os.path.exists("{}/{}".format(basePath, directory)):
        return
    os.system("mkdir -p {}/{}/icon".format(basePath, directory))
    os.system("touch {}/{}/__init__.py".format(basePath, directory))
    with open("{}/{}/__init__.py".format(basePath, directory), "w") as f:
        f.write("import sysconfig\n")
        if not iconFile or not os.path.exists(iconFile):
            iconFile = defaultCategoryIcon
        iconName = os.path.basename(iconFile)
        os.system("cp {} {}/{}/icon/.".format(iconFile, basePath, directory))
        f.write('ICON = "icon/{}"\n'.format(iconName))
        f.write('BACKGROUND ="{}"\n'.format(background))

#!/usr/bin/env python3
import sys, os, re

# pip insists on changing underscores to dashes
# however import will not work with python package names with dashes - this includes the paths
# so until we get rid of pip for dynamic installation we have to support two forms of the name - a dash form for package names and an underscore form for widget names and their paths which is derived from the package name

defaultCategoryIcon = "/icons/misc.png"


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


def entryString(category, directory):
    ret = 'setup(name="{}",\n'.format(category)
    ret += '      packages=["{}"],\n'.format(directory)
    ret += '      package_data={{"{}": ["icons/*.svg"]}},\n'.format(directory)
    ret += '      entry_points={{"orange.widgets": "{} = {}"}},)\n'.format(
        category, directory
    )
    return ret


def makeNewDirectory(basePath, directory, iconFile, background="light-purple"):
    if os.path.exists("{}/{}".format(basePath, directory)):
        return
    os.system("mkdir -p {}/{}/icon".format(basePath, directory))
    os.system("touch {}/{}/__init__.py".format(basePath, directory))
    with open("{}/{}/__init__.py".format(basePath, directory), "w") as f:
        f.write("import sysconfig\n")
        if not iconFile or not os.path.exists(iconFile):
            iconFile = defaultCategoryIcon
        iconName = os.path.basename(iconFile)
        os.system("cp {} {}/{}/icon/.".format(iconFile, basePath, directory))
        f.write('ICON = "icon/{}"\n'.format(iconName))
        f.write('BACKGROUND ="{}"\n'.format(background))
        
def removeCategoriesFromSetupFile(categories,setupFile):
    name_strings=[]
    tempFile=setupFile+"_tmp"
    with open(setupFile, 'r') as f:
        content = f.read()
    parts=content.split("setup(")
    for category in categories:
        name_strings.append('name="{}"'.format(category))
    with open(tempFile,"w") as f:
        for part in parts[1:len(parts)-1]:
            write_part=True
            for name_string in name_strings:
                if (name_string in part):
                    write_part=False
                    break
            if write_part:
                f.write("setup({}\n".format(part))

    os.system("cp {} {}".format(tempFile,setupFile))
    os.system("rm {}".format(tempFile))
    
def removeCategoryFromToolDock(basePath,category,directory):
    setupFile='{}/setup.py'.format(basePath)
    categories=[category]
    removeCategoriesFromSetupFile(categories,setupFile)
    os.system('cd {} && rm -rf {} && rm -rf {}.egg*'.format(basePath,directory,directory))
