#!/bin/bash

# Create the bmp directory if it doesn't exist
mkdir -p ../1/bmp

# Call ./image for each file in the clear directory and then move it to the bmp directory
for file in ../1/clear/*
do
  ./image "$file"
  mv "$file" ../1/bmp/
done
