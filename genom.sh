#!/bin/bash

export listeGenomes=$1
echo $1

mkdir genomes 

while read line  
do   
   IFS=' ' read a b <<< $line
   cd genomes 
   wget $b 
   cd ..
done < $1
