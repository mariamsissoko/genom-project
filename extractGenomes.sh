#!/bin/bash

export repertory=$1


files=$(ls $repertory)
echo files 
cd $repertory

for file in $files;
do
    gunzip $file
done

cd..

