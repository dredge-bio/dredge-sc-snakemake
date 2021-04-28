#!/usr/bin/env python3

import csv
import sys

infile_name = sys.argv[1]

with open(infile_name, 'r') as fp:
    # Skip the header
    fp.readline()

    records = [
        ','.join(reversed((row[1], *row[1].split('|'))))
        for row in csv.reader(fp)
    ]

print('\n'.join(records))
