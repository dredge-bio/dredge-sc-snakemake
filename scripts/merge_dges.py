#!/usr/bin/env python3

import json
import sys
from collections import defaultdict

infile = sys.argv[1]
with open(infile, 'r') as fp:
    data = json.load(fp)

measurements_by_treatment = defaultdict(dict)

for cluster, measurements in data.items():
    for measurement in measurements:
        treatment = measurement.pop('_row').replace('-', '|', 1)

        measurements_by_treatment[treatment][cluster] = measurement

print(json.dumps(measurements_by_treatment))
