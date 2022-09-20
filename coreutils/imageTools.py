import os
import shutil
from glob import glob
from tifffile import imwrite, imread, TiffFile
def splitFilename(origFile):
    (basename,extension)=os.path.splitext(origFile)
    if (extension == ".gz" or extension == ".bz2"):
        (basename,extension2)=os.path.splitext(basename)
        if extension2:
            extension=extension2+extension
    return (basename,extension)
def concatenate(pattern,output):
    filenames=glob(pattern)
    if filenames:
        with open(output,'wb') as wfd:
            for f in filenames:
                with open(f,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)

def splitTiff(filename,numberOfSplitFiles):
    basename,extension=splitFilename(filename)
    numberOfSplitFiles=int(numberOfSplitFiles)
    tif=TiffFile(filename)
    axes = tif.series[0].axes
    imagej_metadata = tif.imagej_metadata
    numberOfSlices=len(tif.pages)
    imagej_metadata['axes']=axes
    sliceSize=int(numberOfSlices/numberOfSplitFiles + 0.5)
    if sliceSize < 1:
        sliceSize=1
    startSlice=0
    fileCount=0
    outputFiles=[]
    while startSlice < (numberOfSlices -1):
        endSlice=startSlice+sliceSize
        if endSlice > (numberOfSlices):
            endSlice=numberOfSlices
        countString="{}".format(fileCount)
        countString=countString.zfill(6)
        outputFile=basename+"_"+countString+extension
        outputFiles.append(outputFile)
        print("working on {}-{} {}".format(startSlice,endSlice,outputFile))
        img=imread(filename,key=range(startSlice,endSlice))
        imwrite(outputFile,img,photometric='minisblack',imagej=True,metadata=imagej_metadata)
        startSlice=endSlice+1
        fileCount+=1
    return outputFiles
