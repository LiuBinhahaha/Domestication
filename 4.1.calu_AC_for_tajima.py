#!/bin/python

# vcf文件添加等位基因计数（Allele count）
# AC，AF 和 AN：AC(Allele Count) 表示该Allele的数目；AF(Allele Frequency) 表示Allele的频率； AN(Allele Number) 表示Allele的总数目。
# 对于1个diploid sample（双倍体）而言：则基因型 0/1 表示sample为杂合子，Allele数为1 (双倍体的sample在该位点只有1个等位基因发生了突变)，Allele的频率为0.5 (双倍体的sample在该位点只有50%的等位基因发生了突变)，总的Allele为2；
# 基因型 1/1 则表示sample为纯合的，Allele数为2，Allele的频率为1，总的Allele为2。

import numpy as np

vcf = open("/vol/liubin/data/NG_1844_selection_analysis/Seita.8G238800/tajimaD/cultivar.vcf", 'r')
out = open("/vol/liubin/data/NG_1844_selection_analysis/Seita.8G238800/tajimaD/cultivar_add_AC.vcf", 'w')

for line in vcf:
    line = line.replace('|', '/')

    if line.startswith("#"):
        out.write(line)
    else:
        _line = line.rstrip().split('\t')
        allele = []
        for i in _line[9:]:
            _i = i.rstrip().split(':')[0]
            allele.append(_i)
        het_AC = allele.count('0/1') + allele.count('1/0')
        hom_AC = 2 * allele.count('1/1')
        AC = str(het_AC + hom_AC)
        out.write('\t'.join(_line[0:7]) + '\t' + f'AC={AC};' + '\t' + '\t'.join(_line[8:]) + '\n')
vcf.close()
out.close()
