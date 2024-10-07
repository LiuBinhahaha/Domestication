import numpy as np

import argparse
import subprocess

parser = argparse.ArgumentParser(description="Extract selective region in Fst result")

parser.add_argument('--input', '-i', type=str, required=True, help="Input fst result file name (e.g., 1094_wild_cul.windowed.weir.fst")
parser.add_argument('--output', '-o', type=str, required=True, help="Ouput file for CMplot (e.g. 1094_wild_cul.windowed.weir.fst_selected_region.txt)")

args = parser.parse_args()

file = open(args.input, 'r')
out = open(args.output, 'w')

def calu_threshold(list):
    threshold = np.quantile(list, 0.99)
    return threshold

fst = []
for line in file:
    if line.startswith("CHROM"):
        out.write(line)
    else:
        _line = line.rstrip().split('\t')
        fst.append(float(_line[5]))
file.seek(0,0)

threshold_99 = calu_threshold(fst)
print(threshold_99)

for line in file:
    if not line.startswith("CHROM"):
        _line = line.rstrip().split('\t')
        if float(_line[5]) > 0 and float(_line[5]) > threshold_99:
            out.write('\t'.join(_line) + '\n')

file.close()
out.close()
