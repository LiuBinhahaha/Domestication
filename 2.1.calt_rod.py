import re
import os
import glob
import numpy as np

import argparse

# 创建解析器
parser = argparse.ArgumentParser(description="calculate Pi(wild)/Pi(cultivated) based on input files.")

# 添加参数
parser.add_argument('--wild', type=str, required=True, help="Input wild file name (e.g., wild_20k/100k_2k.pi.windowed.pi")
parser.add_argument('--cul', type=str, required=True, help="Input cultivated file name (e.g. cultivated_10k/100k_2k.pi.window.pi)")
parser.add_argument('--out', type=str, required=True, help='Output file nmae')

# 解析参数
args = parser.parse_args()

wild = open(args.wild, 'r')
cul = open(args.cul, 'r')
out = open(args.out, 'w')

dic = {}  # wild

for line in wild:
    if line.startswith("CHROM"):
        _line = line.rstrip().split('\t')
        out.write('\t'.join(_line[0:3]) + '\t' + 'ROD' + '\n')
    else:
        _line = line.rstrip().split('\t')
        chr_pos = f'{_line[0]}_{_line[1]}_{_line[2]}'
        # Chr  BIN_start  Bin_end  ROD

        if chr_pos not in dic:
            dic[chr_pos] = []
        dic[chr_pos].append(_line[4])

for line in cul:
    if not l.startswith("CHROM"):
        _l = l.rstrip().split('\t')
        chr_pos = f'{_l[0]}_{_l[1]}_{_l[2]}'
        if chr_pos in dic:
            rod = float(dic[chr_pos][0]) / float(_l[4])
            out.write('\t'.join(chr_pos.split('_')) + '\t' + str(rod) + '\n')

cul.close()
out.close()
wild.close()

