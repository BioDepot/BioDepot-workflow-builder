#!/bin/bash 
chown -R  $1:$2 . 
find ./biodepot  -name __pycache__ -type d -exec rm -r {} +
find ./biodepot  -name *.pyc -type f -exec rm {} +
