import os
import sys
import json
import jsonpickle
import pickle

def unPickleData( filename, jsonFlag=True):
    if jsonFlag:
        sys.stderr.write("opening file {}\n".format(filename))
        try:
             with open(filename, "r") as f:
                data = jsonpickle.decode(f.read())
                f.close()
        except Exception as e:
            with open(filename, "rb") as f:
                data = pickle.load(f)
                f.close()
    else:
        with open(filename, "rb") as f:
            data = pickle.load(f)
        f.close()
    return data


def findWidgetParams(workflowDir):
    widgetsDir=os.path.join(workflowDir, 'widgets')
    for workflowDrawer in os.listdir(widgetsDir):
        workflowDrawerDir=os.path.join(widgetsDir,workflowDrawer)
        for widget in os.listdir(workflowDrawerDir):
            fullWidgetPath=os.path.join(workflowDrawerDir,widget)
            if(os.path.isfile(os.path.join(fullWidgetPath,widget+'.json'))):
                print(widget)
                jsonFile=os.path.join(fullWidgetPath,widget+'.json')
                jsonData=unPickleData(jsonFile)
                print(jsonData['parameters'])
                return

workflowDir='/home/lhhung/ssh_repos/workflows/DNA/GATK_germline_variant/'
findWidgetParams(workflowDir)
