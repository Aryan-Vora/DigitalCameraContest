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
    # Get the file extension
    ext="${file##*.}"

    # Create the directory for this file type if it doesn't exist
    mkdir -p "$1/$ext"

    # Move the file to the corresponding directory
    mv "$file" "$1/$ext/"
done