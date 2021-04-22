#!/usr/bin/env python3

import sys
import numpy as np


def pack_matrix(matrix_filename, output_filename):
    infile_fp = open(matrix_filename, 'r')
    outfile_fp = open(output_filename, 'wb')

    infile_fp.readline()
    infile_fp.readline()

    lines = []

    while line := infile_fp.readline():
        if not line.strip():
            break
        cell_id, cluster_id, expr = line.split()

        lines.append((
            int(cell_id),
            int(cluster_id),
            float(expr),
        ))

    # Sort lines by transcript_id
    lines = sorted(lines, key=lambda x: x[0])

    d = np.dtype([
        ('trancript_id', '<u2'),  # np.int16 -> now np.uint16
        ('cell_id', '<u2'),  # np.int16 -> now np.uint16
        ('expression', '<f4')  # np.float32
    ])

    for line in lines:
        row = np.array([line], dtype=d)
        outfile_fp.write(row.tobytes())

    infile_fp.close()
    outfile_fp.close()


if __name__ == '__main__':
    matrix_filename = sys.argv[1]
    output_filename = sys.argv[2]
    pack_matrix(matrix_filename, output_filename)
