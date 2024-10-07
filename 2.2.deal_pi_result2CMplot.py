import argparse

parser = argparse.ArgumentParser(description="Deal Pi result for CMplot based on input files.")

parser.add_argument('--input', '-i', type=str, required=True, help="Input Pi result file name (e.g., wild_cul_rod.txt")
parser.add_argument('--output', '-o', type=str, required=True, help="Ouput file for CMplot (e.g. wild_cul.for_plot.ROD)")

args = parser.parse_args()

file = open(args.input, 'r')
out = open(args.output, 'w')

for line in file:
    if line.startswith("CHROM"):
        _line = line.rstrip().split('\t')
        out.write('ID\tChr\tPos\tROD\n')
    else:
        _line = line.rstrip().split('\t')
        out.write(f'{_line[0]}_{_line[1]}' + '\t' + _line[0] + '\t' + _line[1] + '\t' + _line[3] + '\n')
file.close()
out.close()

        
