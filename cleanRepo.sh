#!/bin/bash 
chown -R  $1:$2 . 
find .  -name __pycache__ -type d -exec rm -r {} +
find .  -name *.pyc -type f -exec rm {} +
#find coreutils -type f -name *.py -exec black {} \;
