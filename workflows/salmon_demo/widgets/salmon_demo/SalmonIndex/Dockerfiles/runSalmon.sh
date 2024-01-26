#!/bin/bash
directory="$(dirname $inputfile)"
baseDirectory="$(basename $directory)"
echo "salmon $@ -a $inputfile -o $prefix/$baseDirectory/"
eval "salmon $@ -a $inputfile -o $prefix/$baseDirectory/"