#!/bin/bash

# Check if directory is provided
if [ -z "$1" ]
then
  echo "Please provide a directory."
  exit 1
fi
make
# Create the bmp directory if it doesn't exist
mkdir -p $1/bmp
# Call ./image for each file in the clear directory
for file in $1/clear/*.clear
do
  ./image "$file"
done

# Move all .bmp files to the bmp directory
for file in $1/clear/*.bmp
do
  mv -- "$file" $1/bmp/
done

#remove all of the .clear parts of the bmp files in the bmp directory
./onefileextension.sh $1/bmp
