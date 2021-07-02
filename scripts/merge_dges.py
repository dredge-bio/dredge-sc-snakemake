#!/usr/bin/env python3

import json
import sys

infile = sys.argv[1]
with open(infile, 'r') as fp:
    data = json.load(fp)

parsed_measurements = []

for cluster, measurements in data.items():
    for measurement in measurements:
        transcript = measurement.pop('_row').replace('-', '|', 1)

        parsed_measurements.append({
            "clusterID": cluster,
            "transcriptID": transcript,
            "logFC": measurement["avg_log2FC"],
            "pctExpressedCluster": measurement["pct.1"],
            "pctExpressedOther": measurement["pct.2"],
            "pValue": measurement["p_val_adj"],
        })

print(json.dumps(parsed_measurements))
