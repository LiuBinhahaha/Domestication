群体间遗传分化指数，是种群分化和遗传距离的一种衡量方法，分化指数越大，差异越大。适用于亚群体间多样性的比较。用于衡量种群分化程度，取值从0到1，为0则认为两个种群间是随机交配的，基因型完全相似；为1则表示是完全隔离的，完全不相似。它往往从基因的多样性来估计，比如SNP。是一种以哈温平衡为前提的种群遗传学统计方法。
在实际的研究中Fst值为0--0.05时说明群体间遗传分化很小，可以不做考虑；
  0.05--0.15时，表明群体间存在中等程度的遗传分化；
  0.15--0.25时群体间存在较大的遗传分化；
  0.25以上的时候群体间就存在很大的遗传分化了。
vcftools --vcf Seita.2G444300_up_down_500k_snp.vcf \
         --weir-fst-pop c1_sample.txt \
         --weir-fst-pop c2_sample.txt \
         --weir-fst-pop c3_sample.txt \
         --weir-fst-pop w1_sample.txt \
         --weir-fst-pop w2_sample.txt \
         --weir-fst-pop w3_sample.txt \
         --weir-fst-pop w4_sample.txt \
         --out k7_snps \
         --fst-window-size 20000 --fst-window-step 2000  # 窗口20k，步长2k
# --weir-fst-population 群体文件名称(txt格式，一个样本ID一行)
# --fst-window-size 20000 指定窗口大小
# --fst-window-step 2000 指定步长
# 如果计算每个位点的话就不加后面的窗口和步长
vcftools --vcf xxx.vcf --weir-fst-pop xxx.sample1  --weir-fst-pop xxx.sample2 --out xxx.fst --fst-window-size 20000 --fst-window-step 2000
