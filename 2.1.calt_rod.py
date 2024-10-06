import re
import os
import glob
import numpy as np

wild = open("wild_20k_2k.pi.windowed.pi", 'r')
out = open("wild_cul_rod.txt", 'w')

dic = {}

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

fs = ['cul_20k_2k.pi.windowed.pi']

for files in fs:
    print(files)

    file = open(f'{files}', 'r')
    for l in file:
        if not l.startswith("CHROM"):
            _l = l.rstrip().split('\t')
            chr_pos = f'{_l[0]}_{_l[1]}_{_l[2]}'
            if chr_pos in dic:
                rod = float(dic[chr_pos][0])/float(_l[4])
                out.write('\t'.join(chr_pos.split('_')) + '\t' + str(rod) + '\n')
    file.close()
out.close()
wild.close()

