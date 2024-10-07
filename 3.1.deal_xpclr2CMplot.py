import argparse

parser = argparse.ArgumentParser(description="Deal xpclr result for CMplot based on input files.")

parser.add_argument('--input', '-i', type=str, required=True, help="Input xpclr result file name (e.g., broomcorn_win200k_2k.xplclr")
parser.add_argument('--output', '-o', type=str, required=True, help="Ouput file for CMplot (e.g. broomcorn_win200k_2k.for_plot.xpclr)")

args = parser.parse_args()

file = open(args.input, 'r')
out = open(args.output, 'w')

for line in file:
    if line.startswith("id"):
        _line = line.rstrip().split('\t')
        out.write('ID\tChr\tPos_start\txpclr\n')
    else:
        _line = line.rstrip().split('\t')
        if len(_line) == 13:  # Because the legth of some line were less than 13
            out.write(f'{_line[1]}_{_line[4]}' + '\t' + _line[1] + '\t' + _line[4] + '\t' + _line[11] + '\n')
file.close()
out.close()
