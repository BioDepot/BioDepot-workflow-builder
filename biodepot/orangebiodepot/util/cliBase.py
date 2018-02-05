#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  cliBase.py
#
#  Copyright 2018 lhhung <lhhunghimself@gmail.com>
#  

def generateCmd (executables = [], flags = {}, args = []):
    #flags are key value pairs - values begin with = if they are to be joined with the key i.e. --name=John and not --name John 
    cmdStr=''
    for executable in executables:
        cmdStr += executable + ' '
    #short flags no args first
    #short flags args
    #long flags no args
    #long flags args
    #flags no dashs
    #flags no args
    for flagName, flagValue in flags.items():
        if (flagName[:1] == '-') and (flagName[:2] != '--') and (not flagValue):
            cmdStr +=  flagName + ' '
    for flagName,flagValue in flags.items():
        if (flagName[:1] == '-') and (flagName[:2] != '--') and flagValue:
            if(flagValue[:1] == '='):
                cmdStr +=  flagName + flagValue + ' '
            else:
                cmdStr +=  flagName + ' ' + flagValue + ' '
    for flagName, flagValue in flags.items():
        if (flagName[:2] == '--')  and (not flagValue):
            cmdStr +=  flagName + ' '
    for flagName,flagValue in flags.items():
        if (flagName[:2] == '--')  and flagValue:
            if(flagValue[:1] == '='):
                cmdStr +=  flagName + flagValue + ' '
            else:
                cmdStr +=  flagName + ' ' + flagValue + ' '
    for flagName, flagValue in flags.items():
        if (flagName[:1] != '-')  and (not flagValue):
            cmdStr +=  flagName + ' '
    for flagName,flagValue in flags.items():
        if (flagName[:1] == '-')  and flagValue:
            if(flagValue[:1] == '='):
                cmdStr +=  flagName + flagValue + ' '
            else:
                cmdStr +=  flagName + ' ' + flagValue + ' '          
    for arg in args:
        cmdStr += arg + ' '
    #remove any extra space
    if cmdStr:
        cmdStr[:-1]
    return cmdStr
    
    
def main(args):
    executables=('kallisto','quant')
    flags={'-b':'10', '-index':'=transcripts.idx','-kh':'','--pseudoBam':''}
    args={"r1.fq","r2.fq"}
    cmdStr=generateCmd(executables,flags,args)
    print (str(cmdStr))
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
