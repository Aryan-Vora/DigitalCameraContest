#!/bin/bash

# Check if file extension is provided
if [ -z "$1" ]
then
    echo "Please provide a file extension."
    exit 1
fi

# Create the directory if it doesn't exist
mkdir -p "$1"

# Move all files with the provided extension to the corresponding directory
for file in *."$1"
do
    mv "$file" "$1"/
done

