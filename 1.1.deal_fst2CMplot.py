import argparse
import subprocess

parser = argparse.ArgumentParser(description="Deal Fst result for CMplot based on input files.")

parser.add_argument('--input', '-i', type=str, required=True, help="Input fst result file name (e.g., 1094_wild_cul.windowed.weir.fst")
parser.add_argument('--output', '-o', type=str, required=True, help="Ouput file for CMplot (e.g. 1094_wild_cul.windowed.weir.for_plot.fst)")

args = parser.parse_args()

file = open(args.input, 'r')
out = open(args.output, 'w')

for line in file:
    if  line.startswith("CHROM"):
        _line = line.rstrip().split('\t')
        out.write('ID\tChr\tPos\tFST\n')
    else:
        _line = line.rstrip().split('\t')
        if float(_line[5]) >0:
            out.write(f'{_line[0]}_{_line[1]}' + '\t' + _line[0] + '\t' + _line[1] + '\t' + _line[5] + '\n')
file.close()
out.close()


