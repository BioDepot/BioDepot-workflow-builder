#!/usr/bin/env python3
import sys, os, re
def keepWhiteListinToolDock(whiteList,basePath,setupFile):
    tempFile=setupFile+"_tmp"
    with open(setupFile, 'r') as f:
        content = f.read()
    parts=content.split("setup(")
    with open(tempFile,"w") as f:
        for part in parts:
            category=""
            package=""
            categoryMatch=re.search('name\=\"(.*?)\"', part,re.MULTILINE)
            if categoryMatch:
                category=categoryMatch.group(1)
            packageMatch=re.search('packages\=\[\"(.*?)\"', part,re.MULTILINE)
            if packageMatch:
                package=packageMatch.group(1)
            if (category and category in whiteList) or (package and package in whiteList):
                sys.stderr.write("whitelisting {}\n".format(category))
                f.write("setup({}".format(part))
            elif package:
                os.system('cd {} && rm -rf {} && rm -rf {}.egg*'.format(basePath,package,package))
    os.system("cp {} {}".format(tempFile,setupFile))
    os.system("rm {}".format(tempFile))
    return

    
def whiteListToolDock(basePath,categoryList):
    whiteList=(open(categoryList).read().splitlines())
    setupFile='{}/setup.py'.format(basePath)
    keepWhiteListinToolDock(whiteList,basePath,setupFile)

if __name__ == "__main__":
    if len(sys.argv) == 2:
        basePath="/biodepot"
        whiteListToolDock(basePath,sys.argv[1])
