import numpy as np
import argparse
import subprocess

parser = argparse.ArgumentParser(description="Extract selective region in Pi result")

parser.add_argument('--input', '-i', type=str, required=True, help="Input Pi result file name (e.g., wild_cul_rod.txt")
parser.add_argument('--output', '-o', type=str, required=True, help="Selected region in Pi result (e.g. wild_cul_rod.selected_region.txt)")

args = parser.parse_args()

file = open(args.input, 'r')
out = open(args.output, 'w')

def calu_threshold(list):
    threshold = np.quantile(list, 0.99)
    return threshold

ROD = []
for line in file:
    if line.startswith("CHROM"):
        out.write(line)
    else:
        _line = line.rstrip().split('\t')
        ROD.append(float(_line[3]))
file.seek(0,0)

threshold_99 = calu_threshold(ROD)
print(threshold_99)

for line in file:
    if not line.startswith("CHROM"):
        _line = line.rstrip().split('\t')
        if float(_line[3]) > threshold_99:
            out.write('\t'.join(_line) + '\n')

file.close()
out.close()
