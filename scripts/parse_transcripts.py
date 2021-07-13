#!/usr/bin/env python3

import argparse
import csv
import re

parser = argparse.ArgumentParser(
    description='Create a file with different names for transcripts')

parser.add_argument('transcript_file')

parser.add_argument('-p', '--transcript-prefix', required=False)

if __name__ == '__main__':
    args = parser.parse_args()

    if args.transcript_prefix:
        regex = rf'{args.transcript_prefix}'

        def split(x):
            return re.split(regex, x)
    else:
        def split(x):
            return []

    with open(args.transcript_file, 'r') as fp:
        # Skip the header
        fp.readline()

        records = [
            ','.join(x for x in reversed((row[1], *split(row[1]))) if x)
            for row in csv.reader(fp)
        ]

    print('\n'.join(records))
