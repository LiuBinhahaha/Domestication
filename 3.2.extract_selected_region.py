import numpy as np

import argparse
import subprocess

parser = argparse.ArgumentParser(description="Extract selective region in Fst result")

parser.add_argument('--input', '-i', type=str, required=True, help="Input xpclr result file name (e.g., broomcorn_win20k_2k.xpclr")
parser.add_argument('--output', '-o', type=str, required=True, help="Ouput file for selected region (e.g. broomcorn_win20k_2k.xpclr_selected_region.txt)")

args = parser.parse_args()

file = open(args.input, 'r')
out = open(args.output, 'w')

def calu_threshold(list):
    threshold = np.quantile(list, 0.99)
    return threshold

xpclr = []
for line in file:
    if line.startswith("id"):
        out.write(line)
    else:
        _line = line.rstrip().split('\t')
        if len(_line) == 13:
            xpclr.append(float(_line[11]))
file.seek(0,0)

threshold_99 = calu_threshold(xpclr)
print(threshold_99)

for line in file:
    if not line.startswith("id"):
        _line = line.rstrip().split('\t')
        if len(_line) == 13:
            if float(_line[11]) > 0 and float(_line[11]) > threshold_99:
                out.write('\t'.join(_line) + '\n')

file.close()
out.close()
