# domestication
计算群体驯化的方式分为：
1. 群体分歧度检验(Fst, vcftools)
2. 核苷酸多样性检验(Pi, vcftools)
3. 跨群体复合似然比检验(xpclr, vcftools)
4. 群体内选择检验(Tajima's D)
5. 根据SNP验证自然选择(XPEHH)
以上所有方式都是在基因组区域寻找在物种驯化过程中受到选择的基因组区域，区域中又包含了哪些基因，与文献报道的基因相结合。
全基因组驯化结果展示利用CMplot，部分区域的受选择情况展示用ggplot2。
