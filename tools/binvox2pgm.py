#!/usr/bin/python3

#  Code adapted from binvox-rw-py (https://github.com/dimatura/binvox-rw-py)
#  Copyright (C) 2012 Daniel Maturana

import sys
import numpy as np

def read_header(fp):
    """ Read binvox header """
    line = fp.readline().strip()
    if not line.startswith(b'#binvox'):
        raise IOError('Not a binvox file')
    dims = list(map(int, fp.readline().strip().split(b' ')[1:]))
    dims.reverse()
    translate = list(map(float, fp.readline().strip().split(b' ')[1:]))
    scale = list(map(float, fp.readline().strip().split(b' ')[1:]))[0]
    line = fp.readline()
    return dims, translate, scale

def read(fp):
    dims, translate, scale = read_header(fp)
    raw_data = np.frombuffer(fp.read(), dtype=np.uint8)
    values, counts = raw_data[::2], raw_data[1::2]
    data = np.repeat(values, counts).astype(np.bool)
    return data, dims, translate, scale
    

if __name__ == '__main__':
    if len(sys.argv) < 1:
        sys.exit("Usage: python3 binvox2pgm.py input.pgm output.pgm")

    with open(sys.argv[1], 'rb') as f:
        data, dims, translate, scale = read(f)

    output_fn = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1][:-6] + "pgm"        
    with open(output_fn, 'w') as f:
        f.write("P2\n")
        f.write(' '.join(f"{dim}" for dim in dims) + '\n')
        f.write('255\n')
        count = 0
        for v in data:
            c = 1 if v else 0
            f.write(f"{c} ")
            count += 1
            if count == dims[0]:
                f.write("\n")
                count = 0
    print(f"PGM 3D file written to {output_fn}")