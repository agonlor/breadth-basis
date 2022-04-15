# Tools

This folder contains several tools for computing the thickness and breadth pairs.

- **make_object** `make_object.sh` is a bash script that takes a mesh and a size and converts it into a digital object in PGM (3D) format. It uses the three following programs.
- **binvox** [Binvox](https://www.patrickmin.com/binvox/) is a software developed by Patrick Min to voxelize surface meshes.
- **binvox2pgm** `binvox2pgm.py` is a Python script that converts a binvox file into PGM 3D format.
- **pgm2obj** `pgm2obj.py` is a Python script that generates a OBJ mesh from a PGM file.
- **run_experiments** Bash script for computing cycle with different objects, algorithms and parameters
- **lscycles** `lscycles.py` is a Python script that outputs statistics of a list of solutions. Do `python3 lscycles.py ../data/*.json`