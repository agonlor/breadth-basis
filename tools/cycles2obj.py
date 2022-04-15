#!/usr/bin/python3
# Script for making a OBJ mesh from a cycles file in JSON format

import sys
import json

def make_mesh(fn, r):
    """
    Make a OBJ mesh file from a PGM file
    """
    with open(sys.argv[1]) as f:
        cycles = json.load(f)
    fn = fn[:-5] + '_cycles.obj'
    with open(fn, 'w') as fi:
        fi.write(f'# OBJ file made with cycles2obj.py with r = {r}\n')
        n = 0
        i = 1
        for cycle in cycles["cycles"]:
            fi.write(f'o cycle_{i}\n')
            vertices, faces = [], []
            for cube in cycle:
                c = [cube['x'], cube['y'], cube['z']]
                u, v = [], []
                for k in range(3):
                    if c[k] % 2 == 0:
                        u.append(c[k]/2 - r)
                        v.append(c[k]/2 + r)
                    else:
                        u.append(c[k]/2 - 0.5 - r)
                        v.append(c[k]/2 + 0.5 + r)
                vertices.append([u[0], u[1], u[2]])
                vertices.append([v[0], u[1], u[2]])
                vertices.append([u[0], v[1], u[2]])
                vertices.append([v[0], v[1], u[2]])
                vertices.append([u[0], u[1], v[2]])
                vertices.append([v[0], u[1], v[2]])
                vertices.append([u[0], v[1], v[2]])
                vertices.append([v[0], v[1], v[2]])
                faces.append([n+2, n+4, n+8, n+6])
                faces.append([n+1, n+5, n+7, n+3])
                faces.append([n+3, n+7, n+8, n+4])
                faces.append([n+1, n+2, n+6, n+5])
                faces.append([n+5, n+6, n+8, n+7])
                faces.append([n+1, n+3, n+4, n+2])
                n += 8
            for v in vertices:
                fi.write(f'v {v[0]} {v[1]} {v[2]}\n')
            for f in faces:
                fi.write(f'f {f[0]} {f[1]} {f[2]} {f[3]}\n')
            i += 1
    print(f'OBJ file written in {fn} (with r={r})')


if len(sys.argv) < 1:
    sys.exit("Usage: python3 cycles2obj.py object_cycles.json")

r = 0.25
if len(sys.argv) > 2:
    r = float(sys.argv[2])
make_mesh(sys.argv[1], r)