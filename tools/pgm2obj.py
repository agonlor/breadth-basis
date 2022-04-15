#!/usr/bin/python3
# Script for making a OBJ mesh from a PGM 3D file

import sys
import os


def read(fn):
    """
    Read a PGM 3D file and return a 3D array of values
    """
    with open(fn) as f:
        lines = f.readlines()
    for l in list(lines):
        if l[0] == '#':
            lines.remove(l)
    assert lines[0].strip() == 'P2' 
    size = [int(c) for c in lines[1].split()]
    assert lines[2].strip() == '255'
    data = []
    for line in lines[3:]:
        data.extend([int(c) for c in line.split()])
    voxel = []
    n = 0
    for z in range(size[2]):
        voxel_z = []
        for y in range(size[1]):
            voxel_y = data[n:n+size[0]]
            voxel_z.append(voxel_y)
            n += size[0]
        voxel.append(voxel_z)
    assert n == size[0]*size[1]*size[2]
    return (*size, voxel)

def make_mesh(fn):
    """
    Make a OBJ mesh file from a PGM file
    """
    size_x, size_y, size_z, voxel = read(sys.argv[1])
    vertices = list()
    faces = list()
    for x in range(size_x):
        for y in range(size_y):
            for z in range(size_z):
                if voxel[z][y][x] > 0:
                    if x-1 < 0 or voxel[z][y][x-1] == 0: # back face
                        n = len(vertices) + 1
                        vertices.append([x  , y  , z  ])
                        vertices.append([x  , y  , z+1])
                        vertices.append([x  , y+1, z+1])
                        vertices.append([x  , y+1, z  ])
                        faces.append([n, n+1, n+2, n+3])
                    if x+1 >= size_x or voxel[z][y][x+1] == 0: # front face
                        n = len(vertices) + 1
                        vertices.append([x+1, y  , z  ])
                        vertices.append([x+1, y+1, z  ])
                        vertices.append([x+1, y+1, z+1])
                        vertices.append([x+1, y  , z+1])
                        faces.append([n, n+1, n+2, n+3])
                    if y-1 < 0 or voxel[z][y-1][x] == 0: # left face
                        n = len(vertices) + 1
                        vertices.append([x  , y  , z  ])
                        vertices.append([x+1, y  , z  ])
                        vertices.append([x+1, y  , z+1])
                        vertices.append([x  , y  , z+1])
                        faces.append([n, n+1, n+2, n+3])
                    if y+1 >= size_y or voxel[z][y+1][x]== 0: # right face
                        n = len(vertices) + 1
                        vertices.append([x  , y+1, z  ])
                        vertices.append([x  , y+1, z+1])
                        vertices.append([x+1, y+1, z+1])
                        vertices.append([x+1, y+1, z  ])
                        faces.append([n, n+1, n+2, n+3])
                    if z-1 < 0 or voxel[z-1][y][x] == 0: # bottom face
                        n = len(vertices) + 1
                        vertices.append([x  , y  , z  ])
                        vertices.append([x  , y+1, z  ])
                        vertices.append([x+1, y+1, z  ])
                        vertices.append([x+1, y  , z  ])
                        faces.append([n, n+1, n+2, n+3])
                    if z+1 >= size_z or voxel[z+1][y][x] == 0: # top face
                        n = len(vertices) + 1
                        vertices.append([x  , y  , z+1])
                        vertices.append([x+1, y  , z+1])
                        vertices.append([x+1, y+1, z+1])
                        vertices.append([x  , y+1, z+1])
                        faces.append([n, n+1, n+2, n+3])
    fn = os.path.basename(fn)[:-3] + 'obj'
    with open(fn, 'w') as fi:
        fi.write(f'# {len(vertices)} vertices and {len(faces)} faces\n')
        for v in vertices:
            fi.write(f'v {v[0]} {v[1]} {v[2]}\n')
        for f in faces:
            fi.write(f'f {f[0]} {f[1]} {f[2]} {f[3]}\n')
    print(f'OBJ file written in {fn}')


if len(sys.argv) < 1:
    sys.exit("Usage: python3 pgm2obj.py object.pgm")
make_mesh(sys.argv[1])