import numpy as np

import argparse
import subprocess

parser = argparse.ArgumentParser(description="Merge overlap for selective region in Pi result")

parser.add_argument('--input', '-i', type=str, required=True, help="Input selected region Pi result file name (e.g. wild_cul_rod.selected_region.txt")
parser.add_argument('--output', '-o', type=str, required=True, help="Ouput file for merge overlap Pi result (e.g. merged_wild_cul_rod.selected_region.txt")

args = parser.parse_args()

file = open(args.input, 'r')
out = open(args.output, 'w')

dic = {}
for line in file:
    if line.startswith("CHROM"):
        out.write(line)
    else:
        _line = line.rstrip().split('\t')
        if _line[0] not in dic:
            dic[_line[0]] = []
        dic[_line[0]].append(_line[1:])

for chrom in dic:
    rod = dic[chrom]
    merged_intervals = []
    current_interval = rod[0]
    print(current_interval)  # 当前

    rod_value = [float(current_interval[2])]

    for p in range(1, len(rod)):
        prev_start, prev_end = int(current_interval[0]), int(current_interval[1])
        curr_start, curr_end = int(rod[p][0]), int(rod[p][1])
        
        if curr_start <= prev_end:  # 区间有重叠，进行合并
            current_interval[1] = str(max(prev_end, curr_end))  # 更新区间结束的位置
            #current_interval[4] = str((float(current_interval[4]) + float(fst[p][4])) / 2)  # 更新Fst值，取平均
            rod_value.append(float(rod[p][2]))

        else:
            # 不重叠时计算当前区间中的fst均值
            current_interval[2] = str(sum(rod_value) / len(rod_value))

            # 将不重叠的区间添加到merged_intervals
            merged_intervals.append(current_interval)

            # 将当前的 current_interval 更新为新的区间，
            # 也就是下一行的区间 fst[p]，准备在接下来的迭代中检查这个新的区间是否与后面的区间有重叠。
            current_interval = rod[p]
            print(current_interval)  # 更新后

            # 重置fst列表
            rod_value = [float(current_interval[2])]

    # 处理最后一个区间
    if rod_value:
        current_interval[2] = str(np.sum(rod_value) / len(rod_value))
        merged_intervals.append(current_interval)

    for interval in merged_intervals:
        out.write(chrom + '\t' + '\t'.join(interval) + '\n')

file.close()
out.close()



