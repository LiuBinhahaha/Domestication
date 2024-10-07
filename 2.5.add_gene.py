import gzip
import re
import argparse
gff = open("/home/liubin/ref/broomcorn/mizi.gff3", 'r')

dic = {}
for line in gff:
    _line = line.rstrip().split('\t')

    if re.match('gene', _line[2]):
        gene = re.search("ID=(.*?);", _line[8]).group(1)
        if _line[0] not in dic:
            dic[_line[0]] = []
        dic[_line[0]].append(f'{_line[3]}_{_line[4]}_{gene}')
gff.close()

parser = argparse.ArgumentParser(description="Add genes for the selected regions in Pi result")

parser.add_argument('--input', '-i', type=str, required=True, help="Input Pi selected region file name (e.g., merged_wild_cul_rod.selected_region.txt")
parser.add_argument('--output', '-o', type=str, required=True, help="Ouput file for selected regions add genes (e.g. merged_wild_cul_rod.selected_region_add_genes.txt)")

args = parser.parse_args()

file = open(args.input, 'r')
out = open(args.output, 'w')

for line in file:
    if line.startswith("CHROM"):
        _line = line.rstrip().split('\t')
        out.write('\t'.join(_line[0:]) + ' \t' + "gene_annotation" + '\n')
    else:
        genes = []
        _line = line.rstrip().split('\t')
        if _line[0] in dic:
            interval_start = int(_line[1])
            interval_end = int(_line[2])
            
            for v in dic[_line[0]]:
                gene_start = int(v.split('_')[0])
                gene_end = int(v.split('_')[1])

                if gene_start >= interval_start and gene_end <= interval_end:
                    genes.append(v)
        out.write('\t'.join(_line) + '\t' + '\t'.join(genes) + '\n')
 
file.close()
out.close()


