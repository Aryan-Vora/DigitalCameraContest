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
    echo "Processing file: $file"

    # Get the base name of the file without path
    base=$(basename "$file")
    echo "Base name of the file: $base"
    
    # Remove all extensions except the last one
    newname="${base%%.*}.${base##*.}"
    echo "New name for the file: $newname"
    
    # Rename the file only if new name is different
    if [ "$newname" != "$base" ]
    then
        echo "Renaming $file to $newname"
        mv -- "$file" "${file%/*}/$newname"
    else
        echo "Skipping $file as it only has one extension"
    fi
done