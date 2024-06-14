#!/bin/bash

# Check if directory is provided
if [ -z "$1" ]
then
    echo "Please provide a directory."
    exit 1
fi



# Call ./sony_clear.exe for each .SRF file in the provided directory
for file in "$1"/SRF/*.SRF
do
    ./sony_clear.exe "$file"
    echo "Processing $file"
done

#move all clear files to the clear directory
# Create the clear directory if it doesn't exist
mkdir -p "$1"/clear
for file in "$1"/SRF/*.clear
do
    base=$(basename "$file")
    newname="${base%%.*}.clear"
    mv -- "$file" "$1/clear/$newname"
done
