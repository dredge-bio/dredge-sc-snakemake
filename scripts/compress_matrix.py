#!/usr/bin/env python3

import sys
import numpy as np


# Binary file format:
# uint16: the version number of the file
# uint16: the number of float lookups
# ...lookup table of float32 values
# ...table of:
# | uint16: transcriptID | uint16: cellID | uint16: expression_measure_offset |


VERSION_NUMBER = 2


def pack_matrix(matrix_filename, output_filename):
    infile_fp = open(matrix_filename, 'r')
    outfile_fp = open(output_filename, 'wb')

    infile_fp.readline()
    infile_fp.readline()

    lines = []

    expr_measure_offsets = {}

    while line := infile_fp.readline():
        if not line.strip():
            break

        cell_id, cluster_id, expr = line.split()

        expr_measure = float(expr)
        expr_measure = float('{:0.3f}'.format(expr_measure))

        if expr_measure not in expr_measure_offsets:
            expr_measure_offsets[expr_measure] = len(expr_measure_offsets)

        lines.append((
            int(cell_id),
            int(cluster_id),
            expr_measure_offsets[expr_measure],
        ))

    print('expressions: ' + str(len(lines)))
    print('unique expression measurements: ' + str(len(expr_measure_offsets)))
    print('ratio: ' + str(len(expr_measure_offsets) / len(lines)))

    # Sort lines first by expression measure, and then by transcript_id
    # (the former is to improve compression, since there are a lot of
    # repeated expression measures in a row)
    lines = sorted(lines, key=lambda x: x[2])
    lines = sorted(lines, key=lambda x: x[0])

    transcript_expr_entry = np.dtype([
        ('transcript_id', '<u2'),  # np.int16 -> now np.uint16
        ('cell_id', '<u2'),  # np.int16 -> now np.uint16
        ('expression_measure_offset', '<u2')  # np.float32
    ])

    # Version of file format
    outfile_fp.write(
        np.uint16(VERSION_NUMBER).tobytes())

    # Length of expression measure lookup table
    outfile_fp.write(
        np.uint16(len(expr_measure_offsets)).tobytes())

    outfile_fp.write(
        np.array(list(expr_measure_offsets.keys()), dtype='<f4').tobytes())

    for line in lines:
        row = np.array([line], dtype=transcript_expr_entry)
        outfile_fp.write(row.tobytes())

    infile_fp.close()
    outfile_fp.close()


if __name__ == '__main__':
    matrix_filename = sys.argv[1]
    output_filename = sys.argv[2]
    pack_matrix(matrix_filename, output_filename)
