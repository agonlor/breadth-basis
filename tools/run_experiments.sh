#!/bin/bash

######################################################################
# Script Name: run_experiments
# Description: Compute many homology bases
# Author     : Aldo Gonzalez-Lorenzo                                                
# Email      : aldo.gonzalez-lorenzo@univ-amu.fr                                          
######################################################################

for file in  "Eight" "Neptune" "Dancing" "Pegasus" "Fertility" "Amphora"
do
  for n in 200 #40 80 160
  do
    echo $file $n
    # tools/make_object.sh data/mesh/${file}* $n

    build/shorthomologybasis data/${file}_${n}.pgm
    sleep 2
    echo; echo
    
    for seed in {10..19}
    do
     build/shorthomologybasis data/${file}_${n}.pgm --dey -s $seed
     sleep 2
     echo; echo
    done
  done
done