#!/bin/bash 
find ./biodepot  -name __pycache__ -type d -exec rm -r {} +
find ./biodepot  -name *.pyc -type f -exec rm {} +
