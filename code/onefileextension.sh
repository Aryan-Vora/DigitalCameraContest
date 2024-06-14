#!/bin/bash

# Check if directory is provided
if [ -z "$1" ]
then
    echo "Please provide a directory."
    exit 1
fi

# Iterate over all files in the provided directory
for file in "$1"/*
do
    # Get the base name of the file without path
    base=$(basename "$file")
    
    # Remove all extensions except the last one
    newname="${base%%.*}.${base##*.}"
    
    # Rename the file
    mv -- "$file" "${file%/*}/$newname"
done
