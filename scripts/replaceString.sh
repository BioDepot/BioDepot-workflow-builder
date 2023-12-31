#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 [-d] directory pattern1 pattern2"
    echo "  -d: Dry run. Only display matches without replacing."
}

# Default values
dry_run=0

# Parse options
while getopts ":d" opt; do
  case $opt in
    d)
      dry_run=1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

# Assign the directory and patterns to variables
dir=$1
pat1=$2
pat2=$3

# Check if enough arguments are provided
if [ "$#" -ne 3 ]; then
    usage
    exit 1
fi

# Function to process files
process_files() {
    find "$dir" -type f -exec grep -lZ "$pat1" {} + |
    while IFS= read -r -d '' file; do
        if [ "$dry_run" -eq 1 ]; then
            echo "Would process $file"
            grep --color=always "$pat1" "$file"
        else
            sed -i "s/$pat1/$pat2/g" "$file"
            echo "Processed $file"
        fi
    done
}

# Process files based on dry run flag
process_files
