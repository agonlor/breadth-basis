#!/bin/bash

######################################################################
# Script Name: script_tb_basis
# Description: Take a mesh and convert it into a PGM file
# Args       : The mesh file name and the size (integer)
# Author     : Aldo Gonzalez-Lorenzo                                                
# Email      : aldo.gonzalez-lorenzo@univ-amu.fr                                          
######################################################################

# -- Arguments
if [ $# -ne 2 ]
then
  echo "Usage: $0 <mesh file> <precision>"
  exit 1
fi
echo "make_object.sh: file = "$1" precision = "$2

filename=$(basename -- "$1")
filename="${filename%.*}_"$2.pgm              # output filename


# -- Voxelize the mesh
tools/binvox -d $2 -v -fit $1                 # mesh -> binvox
tools/binvox2pgm data/mesh/*.binvox $filename # binvox -> pgm
python3 tools/pgm2obj.py $filename            # pgm -> obj
rm data/mesh/*.binvox
build/breadth-basis $filename --voxdist