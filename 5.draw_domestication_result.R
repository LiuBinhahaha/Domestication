setwd("/vol/liubin/data/dom_analysis/sorghum/pi/win20k_2k/")
library(ggplot2)
library(CMplot)
df <- read.table("wild_cul.for_plot.ROD", sep = '\t', header = T)
head(df)

# ggplot手动画
ggplot(df, aes(x=BIN_START, y=Fst,color= Chr)) + 
  geom_point( aes(color=as.factor(BIN_START)), alpha=0.8, size=0.3)+
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous( label = df$BIN_START, breaks = df$BIN_START) +
  facet_wrap(~BIN_START, ncol = 18) +  # 分面
  ylab("Fst") +
  xlab("Chr") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(file="Fst_1844.pdf", dpi = 600, width = 20, height = 5)


# 阈值线确定(95/99%分位线)
threshold =  quantile(na.omit(df$ROD),probs=seq(0,1,0.01))['99%']
print(threshold)

CMplot(df, type="h",plot.type="m", band=0.9, LOG10=FALSE, ylab="Pi(wild/cultivated)",
       threshold=2.963754, threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.6,
       col= c("#E57373","#80CBC4","#D28FDE","#81D5F9"), chr.den.col=NULL, file="pdf", dpi=300, file.output=TRUE,verbose=TRUE, highlight.text.xadj)  #ylim = c(0.1,1)
# 画图前改变ylab标签、阈值


#################################################################################################################################################
# 画一个基因的部分区域,Fst和Pi, TajimaD, XP-CLR都适用
mycolor1 <- c("#3AB5B3","#954F97","#B3D5F1","#F0686C","#5567B1","#80C66D")  # 蓝、绿、紫、橙黄
mycolor2 <- c("#F2A1A7","#7DC69B","#9BD7F3","#FBDDDD")

setwd("/vol/liubin/data/NG_1844_selection_analysis/Seita.8G238800")
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(showtext)  # 导出字体
library(patchwork)
library(cowplot)
library(ggprism)  # 改变坐标轴

fst <- read.table("Seita.8G238800_up_down_500k.fst.txt", sep='\t', header = T)

head(fst)

threshold_fst =  quantile(na.omit(fst$FST),probs=seq(0,1,0.01))['95%'] 
print(threshold_fst)

# 调色
colourCount = length(unique(fst$Type))
getPalette = colorRampPalette(brewer.pal(12, "Paired")) # 9 Set1

# 箭头标注基因位置
arrow_x = 39.164001
arrow_y_start = threshold_fst
arrow_y_end = arrow_y_start * 1.3


p1 <- 
  ggplot(fst, aes(x=BIN_START/1000000, y=FST,color = Type)) +
  geom_line(aes(fill=Type), linewidth = 0.5) +
  geom_hline(yintercept = threshold_fst, colour="black", linetype="dashed", linewidth = 0.43) + # geom_vline(xintercept = 49.119592)添加垂线
  
  # 添加箭头
  annotate("segment", x = arrow_x, xend = arrow_x, y =arrow_y_end, yend = arrow_y_start, 
           arrow = arrow(length = unit(0.3, "cm")), color = "black")+
  
  # 色盘添加颜色
  scale_color_manual(values = mycolor1)+
  
  theme_bw() +
  theme(panel.grid = element_blank())+ # theme(legend.position = 'none')去除图例
  theme(axis.text.x = element_text(angle = 30,vjust = 0.55,size=9, color="black", hjust = 0.55), 
        panel.border = element_rect(color = "black", size=0.43),
        text=element_text(size=9, color="black"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.43))+ # 调整x轴
  theme(axis.ticks.x=element_line(color="black",size=0.43,lineend = 1))+  # 调整x轴刻度粗细
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.43))+ # 调整y轴
  theme(axis.ticks.y=element_line(color="black",size=0.2,lineend = 0.43)) + # 调整y轴刻度线
  theme(panel.border = element_blank(), #去除坐标轴线默认填充的灰色
        axis.text.x = element_text(angle = 0,  # family= "Arial"
                                   vjust = 0.55,
                                   size = 9, 
                                   colour = "black"),
        axis.text.y = element_text(vjust=0.55, 
                                   size=9, 
                                   color="black", 
                                   hjust=0.55),
        axis.line = element_line(colour = "black",size=0.2, lineend = "square"), #将x=0轴和y=0轴加粗显示(size=0.3), x、y交界处是完美交接(lineend = "square"语句
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))+
  theme(legend.position = "right") +
  theme(legend.key.size = unit(30, "pt")) +  # 图例大小
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title=element_blank()) +  #  隐藏图例标题
  #scale_y_continuous(limits = c(0, 0.8), guide = "prism_offset") +
  labs(x="Position(Mb)", y="FST") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

print(p1)
ggsave(p1, file = "c1-c2-c3-wild.fst.pdf", height = 5, width = 10)


# 核苷酸多样性
pi <- read.table("Seita.8G238800_up_down_500k.pi.txt", header = T, sep = '\t')

threshold_wild = quantile(na.omit(pi_wild$PI),probs=seq(0,1,0.01))['95%']
print(threshold_wild)
# 调色
colourCount = length(unique(pi$Type))
getPalette = colorRampPalette(brewer.pal(12, "Paired")) # 9 Set1

# 箭头标注基因位置
arrow_x = 39.164001
arrow_y_start = threshold_pi
arrow_y_end = arrow_y_start * 1.3


p2 <- 
  ggplot(pi, aes(x=BIN_START/1000000, y=Pi,color = Type)) +
  geom_line(aes(fill=Type), linewidth = 0.5) +
  geom_vline(xintercept = 39.164001) + #添加垂线
  
  # 添加箭头
  #annotate("segment", x = arrow_x, xend = arrow_x, y =arrow_y_end, yend = arrow_y_start, 
  #         arrow = arrow(length = unit(0.3, "cm")), color = "black")+
  
  # 色盘添加颜色
  scale_fill_manual(values = c("#898989", "#ea9da3", "#79bf97", "#99d1eb"))+
  
  theme_bw() +
  theme(panel.grid = element_blank())+ # theme(legend.position = 'none')去除图例
  theme(axis.text.x = element_text(angle = 30,vjust = 0.55,size=9, color="black", hjust = 0.55), 
        panel.border = element_rect(color = "black", size=0.43),
        text=element_text(size=9, color="black"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.43))+ # 调整x轴
  theme(axis.ticks.x=element_line(color="black",size=0.43,lineend = 1))+  # 调整x轴刻度粗细
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.43))+ # 调整y轴
  theme(axis.ticks.y=element_line(color="black",size=0.43,lineend = 10)) + # 调整y轴刻度线
  theme(panel.border = element_blank(), #去除坐标轴线默认填充的灰色
        axis.text.x = element_text(angle = 0,
                                   vjust = 0.55,
                                   size = 9, 
                                   colour = "black"),
        axis.text.y = element_text(vjust=0.55, 
                                   size=9, 
                                   color="black", 
                                   hjust=0.55),
        axis.line = element_line(colour = "black",size=0.2, lineend = "square"), #将x=0轴和y=0轴加粗显示(size=0.3), x、y交界处是完美交接(lineend = "square"语句
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))+
  theme(legend.position = "right") +
  theme(legend.key.size = unit(30, "pt")) +  # 图例大小
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title=element_blank()) +  #  隐藏图例标题
  #scale_y_continuous(limits = c(0, 0.08), guide = "prism_offset")
  labs(x="Position(Mb)", y="Pi") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

print(p2)

ggsave(p2, file="c1-c2-c3-wild.pi.pdf",  height=2.5, width = 4.5)

# 拼图
combine <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 2, align = "v")
print(combine)
ggsave("Seita.8G238800_up_down_500k_fst_pi.pdf", plot = combine, width = 12.7, height = 2.5)



# 两个hap
fst_hap <- read.table("Seita.2G444300_hap1vshap2.windowed.weir.fst", sep='\t', header = T)

threshold_fst_hap =  quantile(na.omit(fst_hap$MEAN_FST),probs=seq(0,1,0.01))['95%'] 
print(threshold_fst_hap)

# 箭头标注基因位置
arrow_x = 49.119592
arrow_y_start = threshold_fst_hap
arrow_y_end = arrow_y_start * 1.3

p3 <- 
  ggplot(fst_hap, aes(x=BIN_START/1000000, y=MEAN_FST)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = threshold_fst_hap, colour="black", linetype="dashed", linewidth = 0.43) + # geom_vline(xintercept = 49.119592)添加垂线
  
  # 添加箭头
  annotate("segment", x = arrow_x, xend = arrow_x, y =arrow_y_end, yend = arrow_y_start, 
           arrow = arrow(length = unit(0.3, "cm")), color = "black")+
  
  theme_bw() +
  theme(panel.grid = element_blank())+ # theme(legend.position = 'none')去除图例
  theme(axis.text.x = element_text(angle = 30,vjust = 0.55,size=9, color="black", hjust = 0.55), 
        panel.border = element_rect(color = "black", size=0.43),
        text=element_text(size=9, color="black"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.43))+ # 调整x轴
  theme(axis.ticks.x=element_line(color="black",size=0.43,lineend = 1))+  # 调整x轴刻度粗细
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.43))+ # 调整y轴
  theme(axis.ticks.y=element_line(color="black",size=0.2,lineend = 0.43)) + # 调整y轴刻度线
  theme(panel.border = element_blank(), #去除坐标轴线默认填充的灰色
        axis.text.x = element_text(angle = 0,  # family= "Arial"
                                   vjust = 0.55,
                                   size = 9, 
                                   colour = "black"),
        axis.text.y = element_text(vjust=0.55, 
                                   size=9, 
                                   color="black", 
                                   hjust=0.55),
        axis.line = element_line(colour = "black",size=0.2, lineend = "square"), #将x=0轴和y=0轴加粗显示(size=0.3), x、y交界处是完美交接(lineend = "square"语句
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))+
  theme(legend.position = "right") +
  theme(legend.key.size = unit(30, "pt")) +  # 图例大小
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title=element_blank()) +  #  隐藏图例标题
  #scale_y_continuous(limits = c(0, 0.8), guide = "prism_offset") +
  labs(x="Position(Mb)", y="FST(Hap1vsHap2)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
print(p3)


###############################################################

# 画多个亚群Tajima'D值的折线图在一张图上
mycolor3 <- c("#A0ADD0", "#E47178", "#F5DC75", "grey")
setwd("/vol/liubin/data/NG_1844_selection_analysis/select_sweep_analysis/tajimaD/")
tajima <- read.csv("c1_c2_c3_wild_tajima_20k_2k.csv", header = T, sep = ",")

threshold_t =  quantile(na.omit(tajima$TajimaD),probs=seq(0,1,0.01))['95%'] 
print(threshold_t)

# 箭头标注基因位置
arrow_x = 49.114000
arrow_y_start = threshold_t
arrow_y_end = arrow_y_start * 1.3
head(tajima)

p5 <- 
  ggplot(tajima, aes(x=BIN_START/1000000, y=TajimaD, color = Type)) +
  geom_line(aes(color = Type), linewidth = 0.5) +
  geom_vline(xintercept = 49.114000, linetype = "dashed", linewidth = 0.43) +
  geom_hline(yintercept = 0, linewidth = 0.43) +
  
  # 添加箭头
  #annotate("segment", x = arrow_x, xend = arrow_x, y =arrow_y_end, yend = arrow_y_start, 
  #         arrow = arrow(length = unit(0.3, "cm")), color = "black")+
  
  # 色盘添加颜色
  scale_color_manual(values = mycolor3)+
  
  theme_bw() +
  theme(panel.grid = element_blank())+ # theme(legend.position = 'none')去除图例
  theme(axis.text.x = element_text(angle = 30,vjust = 0.55,size=9, color="black", hjust = 0.55), 
        panel.border = element_rect(color = "black", size=0.43),
        text=element_text(size=9, color="black"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.43))+ # 调整x轴
  theme(axis.ticks.x=element_line(color="black",size=0.43,lineend = 1))+  # 调整x轴刻度粗细
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.43))+ # 调整y轴
  theme(axis.ticks.y=element_line(color="black",size=0.43,lineend = 10)) + # 调整y轴刻度线
  theme(panel.border = element_blank(), #去除坐标轴线默认填充的灰色
        axis.text.x = element_text(angle = 0,
                                   vjust = 0.55,
                                   size = 9, 
                                   colour = "black"),
        axis.text.y = element_text(vjust=0.55, 
                                   size=9, 
                                   color="black", 
                                   hjust=0.55),
        axis.line = element_line(colour = "black",size=0.2, lineend = "square"), #将x=0轴和y=0轴加粗显示(size=0.3), x、y交界处是完美交接(lineend = "square"语句
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))+
  theme(legend.position = "right") +
  theme(legend.key.size = unit(30, "pt")) +  # 图例大小
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title=element_blank()) +  #  隐藏图例标题
  #scale_y_continuous(limits = c(0, 0.08), guide = "prism_offset")
  labs(x="Position(Mb)", y="Tajima'D") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
print(p5)


mycolor3 <- c("#A0ADD0", "grey")
setwd("/vol/liubin/data/NG_1844_selection_analysis/select_sweep_analysis/tajimaD/")
tajima <- read.csv("cultivar_wild_tajima_20k_2k.csv", header = T, sep = ",")

threshold_t =  quantile(na.omit(tajima$TajimaD),probs=seq(0,1,0.01))['95%'] 
print(threshold_t)

# 箭头标注基因位置
arrow_x = 49.114000
arrow_y_start = threshold_t
arrow_y_end = arrow_y_start * 1.3
head(tajima)

p6 <- 
  ggplot(tajima, aes(x=BIN_START/1000000, y=TajimaD, color = Type)) +
  geom_line(aes(color = Type), linewidth = 0.5) +
  geom_vline(xintercept = 49.114000, linetype = "dashed", linewidth = 0.43) +
  geom_hline(yintercept = 0, linewidth = 0.43) +
  
  # 添加箭头
  #annotate("segment", x = arrow_x, xend = arrow_x, y =arrow_y_end, yend = arrow_y_start, 
  #         arrow = arrow(length = unit(0.3, "cm")), color = "black")+
  
  # 色盘添加颜色
  scale_color_manual(values = mycolor3)+
  
  theme_bw() +
  theme(panel.grid = element_blank())+ # theme(legend.position = 'none')去除图例
  theme(axis.text.x = element_text(angle = 30,vjust = 0.55,size=9, color="black", hjust = 0.55), 
        panel.border = element_rect(color = "black", size=0.43),
        text=element_text(size=9, color="black"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.43))+ # 调整x轴
  theme(axis.ticks.x=element_line(color="black",size=0.43,lineend = 1))+  # 调整x轴刻度粗细
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.43))+ # 调整y轴
  theme(axis.ticks.y=element_line(color="black",size=0.43,lineend = 10)) + # 调整y轴刻度线
  theme(panel.border = element_blank(), #去除坐标轴线默认填充的灰色
        axis.text.x = element_text(angle = 0,
                                   vjust = 0.55,
                                   size = 9, 
                                   colour = "black"),
        axis.text.y = element_text(vjust=0.55, 
                                   size=9, 
                                   color="black", 
                                   hjust=0.55),
        axis.line = element_line(colour = "black",size=0.2, lineend = "square"), #将x=0轴和y=0轴加粗显示(size=0.3), x、y交界处是完美交接(lineend = "square"语句
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))+
  theme(legend.position = "right") +
  theme(legend.key.size = unit(30, "pt")) +  # 图例大小
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title=element_blank()) +  #  隐藏图例标题
  #scale_y_continuous(limits = c(0, 0.08), guide = "prism_offset")
  labs(x="Position(Mb)", y="Tajima'D") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

print(p6)

combine <- plot_grid(p5, p6, labels = c("A", "B"), ncol = 2, align = "v")
print(combine)
ggsave("Seita.2G444300_up_down_500k_tajimaD_c1_c2_c3_cultivar_wild_20k_2k.pdf", plot = combine, width = 12.2, height = 2.5)






























##################################################################################################################################################


# 画Tajima'D值图单个亚群一张图
setwd("/vol/liubin/data/NG_1844_selection_analysis/select_sweep_analysis/tajimaD")
tajima1 <- read.table("c1_19k_TajimaD.Tajima.D_add_95_interval.txt", header = T, sep = "\t")


threshold_t1 =  quantile(na.omit(tajima1$TajimaD),probs=seq(0,1,0.01))['95%'] 
print(threshold_t1)

# 箭头标注基因位置
arrow_x = 49.115000
arrow_y_start = threshold_t1
arrow_y_end = arrow_y_start * 0.5

tajima1$ci_lower <- as.numeric(tajima1$ci_lower)
tajima1$ci_upper <- as.numeric(tajima1$ci_upper)
p5 <- 
  ggplot(tajima1, aes(x=BIN_START/1000000, y=TajimaD)) +
  geom_line(color = "#A0ADD0", linewidth = 0.5) +
  
  # 手动添加95%置信区间
  geom_ribbon(data = tajima1, aes(x=BIN_START/1000000, ymin=ci_lower, ymax=ci_upper), fill="#A0ADD0", alpha=0.5) +

  # geom_smooth自动计算添加95%置信区间的平滑曲线
  #geom_smooth(method = "loess", se = TRUE, color = "blue", fill = "#A0ADD0", alpha = 0.3) +
  geom_vline(xintercept = 49.119592, linetype = "dashed")+
  
  # 添加箭头
  #annotate("segment", x = arrow_x, xend = arrow_x, y =arrow_y_end, yend = arrow_y_start, 
           #arrow = arrow(length = unit(0.3, "cm")), color = "black")+
  
  
  
  theme_bw() +# theme(legend.position = 'none')去除图例
  theme(axis.text.x = element_text(angle = 30,vjust = 0.55,size=9, color="black", hjust = 0.55), 
        panel.border = element_rect(color = "black", size=0.43),
        text=element_text(size=9, color="black"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.43))+ # 调整x轴
  theme(axis.ticks.x=element_line(color="black",size=0.43,lineend = 1))+  # 调整x轴刻度粗细
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.43))+ # 调整y轴
  theme(axis.ticks.y=element_line(color="black",size=0.43,lineend = 10)) + # 调整y轴刻度线
  theme(panel.border = element_blank(), #去除坐标轴线默认填充的灰色
        axis.text.x = element_text(angle = 0,
                                   vjust = 0.55,
                                   size = 9, 
                                   colour = "black"),
        axis.text.y = element_text(vjust=0.55, 
                                   size=9, 
                                   color="black", 
                                   hjust=0.55),
        axis.line = element_line(colour = "black",size=0.2, lineend = "square"), #将x=0轴和y=0轴加粗显示(size=0.3), x、y交界处是完美交接(lineend = "square"语句
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))+
  theme(legend.position = "right") +
  theme(legend.key.size = unit(30, "pt")) +  # 图例大小
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title=element_blank()) +  #  隐藏图例标题
  #scale_y_continuous(limits = c(0, 0.08), guide = "prism_offset")
  labs(x="Position(Mb)", y="Tajima'D") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

print(p5)
# ggsave("c1_tajimaD.pdf", plot = p5, width = 7.5, height = 4)



tajima2 <- read.table("c2_19k_TajimaD.Tajima.D_add_95_interval.txt", header = T, sep = "\t")

threshold_t2 =  quantile(na.omit(tajima2$TajimaD),probs=seq(0,1,0.01))['95%'] 
print(threshold_t2)

# 箭头标注基因位置
arrow_x = 49.115000
arrow_y_start = threshold_t2
arrow_y_end = arrow_y_start * 0.5


p6 <- 
  ggplot(tajima2, aes(x=BIN_START/1000000, y=TajimaD)) +
  geom_line(color = "#E47178", linewidth = 0.5) +
  
  #geom_smooth自动计算添加95%的置信区间
  #geom_smooth(method = "loess", se = TRUE, color = "blue", fill = "#A0ADD0", alpha = 0.3) +
  geom_ribbon(data = tajima2, aes(ymin=ci_lower, ymax=ci_upper), fill="#E47178", alpha=0.5) +
  geom_vline(xintercept = 49.1195952, linetype = "dashed") +
  
  # 添加箭头
  #annotate("segment", x = arrow_x, xend = arrow_x, y =arrow_y_end, yend = arrow_y_start, 
  #         arrow = arrow(length = unit(0.3, "cm")), color = "black")+
  
  
  theme_bw() +
  #theme(panel.grid = element_blank())+去除背景网格线 # theme(legend.position = 'none')去除图例
  theme(axis.text.x = element_text(angle = 30,vjust = 0.55,size=9, color="black", hjust = 0.55), 
        panel.border = element_rect(color = "black", size=0.43),
        text=element_text(size=9, color="black"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.43))+ # 调整x轴
  theme(axis.ticks.x=element_line(color="black",size=0.43,lineend = 1))+  # 调整x轴刻度粗细
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.43))+ # 调整y轴
  theme(axis.ticks.y=element_line(color="black",size=0.43,lineend = 10)) + # 调整y轴刻度线
  theme(panel.border = element_blank(), #去除坐标轴线默认填充的灰色
        axis.text.x = element_text(angle = 0,
                                   vjust = 0.55,
                                   size = 9, 
                                   colour = "black"),
        axis.text.y = element_text(vjust=0.55, 
                                   size=9, 
                                   color="black", 
                                   hjust=0.55),
        axis.line = element_line(colour = "black",size=0.2, lineend = "square"), #将x=0轴和y=0轴加粗显示(size=0.3), x、y交界处是完美交接(lineend = "square"语句
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))+
  theme(legend.position = "right") +
  theme(legend.key.size = unit(30, "pt")) +  # 图例大小
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title=element_blank()) +  #  隐藏图例标题
  #scale_y_continuous(limits = c(0, 0.08), guide = "prism_offset")
  labs(x="Position(Mb)", y="Tajima'D") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

print(p6)
# ggsave("c2_tajimaD.pdf", plot = p6, width = 7.5, height = 4)



tajima3 <- read.table("c3_19k_TajimaD.Tajima.D_add_95_interval.txt", header = T, sep = "\t")

threshold_t3 =  quantile(na.omit(tajima3$TajimaD),probs=seq(0,1,0.01))['95%'] 
print(threshold_t3)

# 箭头标注基因位置
arrow_x = 49.115000
arrow_y_start = threshold_t3
arrow_y_end = arrow_y_start * 1.3

p7 <- 
  ggplot(tajima3, aes(x=BIN_START/1000000, y=TajimaD)) +
  geom_line(color = "#F5DC75", linewidth = 0.5) +
  
  #geom_smooth自动计算添加95%的置信区间
  #geom_smooth(method = "loess", se = TRUE, color = "blue", fill = "#A0ADD0", alpha = 0.3) +
  geom_ribbon(data = tajima3, aes(ymin=ci_lower, ymax=ci_upper), fill="#F5DC75", alpha=0.5) +
  geom_vline(xintercept = 49.1195952, linetype = "dashed") +
  
  # 添加箭头
  #annotate("segment", x = arrow_x, xend = arrow_x, y =arrow_y_end, yend = arrow_y_start, 
  #         arrow = arrow(length = unit(0.3, "cm")), color = "black")+
  
  
  theme_bw() +
  #theme(panel.grid = element_blank())+去除背景网格线 # theme(legend.position = 'none')去除图例
  theme(axis.text.x = element_text(angle = 30,vjust = 0.55,size=9, color="black", hjust = 0.55), 
        panel.border = element_rect(color = "black", size=0.43),
        text=element_text(size=9, color="black"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.43))+ # 调整x轴
  theme(axis.ticks.x=element_line(color="black",size=0.43,lineend = 1))+  # 调整x轴刻度粗细
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.43))+ # 调整y轴
  theme(axis.ticks.y=element_line(color="black",size=0.43,lineend = 10)) + # 调整y轴刻度线
  theme(panel.border = element_blank(), #去除坐标轴线默认填充的灰色
        axis.text.x = element_text(angle = 0,
                                   vjust = 0.55,
                                   size = 9, 
                                   colour = "black"),
        axis.text.y = element_text(vjust=0.55, 
                                   size=9, 
                                   color="black", 
                                   hjust=0.55),
        axis.line = element_line(colour = "black",size=0.2, lineend = "square"), #将x=0轴和y=0轴加粗显示(size=0.3), x、y交界处是完美交接(lineend = "square"语句
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))+
  theme(legend.position = "right") +
  theme(legend.key.size = unit(30, "pt")) +  # 图例大小
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title=element_blank()) +  #  隐藏图例标题
  #scale_y_continuous(limits = c(0, 0.08), guide = "prism_offset")
  labs(x="Position(Mb)", y="Tajima'D") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

print(p7)
# ggsave("c3_tajimaD.pdf", plot = p7, width = 7.5, height = 4)


tajima4 <- read.table("wild_100k_TajimaD.Tajima.D", header = T, sep = "\t")


threshold_t4 =  quantile(na.omit(tajima4$TajimaD),probs=seq(0,1,0.01))['95%'] 
print(threshold_t4)

# 箭头标注基因位置
arrow_x = 39.164922
arrow_y_start = threshold_t3
arrow_y_end = arrow_y_start * 1.3


p8 <- 
  ggplot(tajima4, aes(x=BIN_START/1000000, y=TajimaD)) +
  geom_line(color = "#59B78F", linewidth = 0.5) +
  
  # 添加箭头
  annotate("segment", x = arrow_x, xend = arrow_x, y =arrow_y_end, yend = arrow_y_start, 
           arrow = arrow(length = unit(0.3, "cm")), color = "black")+
  
  
  theme_bw() +
  theme(panel.grid = element_blank())+ # theme(legend.position = 'none')去除图例
  theme(axis.text.x = element_text(angle = 30,vjust = 0.55,size=9, color="black", hjust = 0.55), 
        panel.border = element_rect(color = "black", size=0.43),
        text=element_text(size=9, color="black"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.43))+ # 调整x轴
  theme(axis.ticks.x=element_line(color="black",size=0.43,lineend = 1))+  # 调整x轴刻度粗细
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.43))+ # 调整y轴
  theme(axis.ticks.y=element_line(color="black",size=0.43,lineend = 10)) + # 调整y轴刻度线
  theme(panel.border = element_blank(), #去除坐标轴线默认填充的灰色
        axis.text.x = element_text(angle = 0,
                                   vjust = 0.55,
                                   size = 9, 
                                   colour = "black"),
        axis.text.y = element_text(vjust=0.55, 
                                   size=9, 
                                   color="black", 
                                   hjust=0.55),
        axis.line = element_line(colour = "black",size=0.2, lineend = "square"), #将x=0轴和y=0轴加粗显示(size=0.3), x、y交界处是完美交接(lineend = "square"语句
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))+
  theme(legend.position = "right") +
  theme(legend.key.size = unit(30, "pt")) +  # 图例大小
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title=element_blank()) +  #  隐藏图例标题
  #scale_y_continuous(limits = c(0, 0.08), guide = "prism_offset")
  labs(x="Position(Mb)", y="Tajima'D") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

print(p8)


combine <- plot_grid(p1, p5, p2, p6, p3, p7, labels = c("A", "B", "C", "D", "E", "F"), nrow = 3, align = "v")
combine <- plot_grid(p5,p5, p7,p7, p6,p6, labels = "AUTO", ncol = 2, align = "h", hjust = 0, vjust = 0)
print(combine)
ggsave("Seita.2G444300_up_down_500k_window_tajimaD_win19k.pdf", plot = combine, width =9.2, height = 6)


###### 画XP-CLR
setwd("/vol/liubin/data/NG_1844_selection_analysis/select_sweep_analysis/xpclr")
df <- read.table("Seita.2G444300_up_down_50k_20k.xpclr", header = T, sep = "\t")


threshold_t1 =  quantile(na.omit(df$xpclr),probs=seq(0,1,0.01))['95%'] 
print(threshold_t1)

# 箭头标注基因位置
arrow_x = 49.118750
arrow_y_start = threshold_t1
arrow_y_end = arrow_y_start * 0.5

p9 <- 
  ggplot(df, aes(x=start/1000000, y=xpclr)) +
  geom_bar(fill = "#326291", stat = "identity") +
  geom_vline(xintercept = 49.118750, linetype = "dashed", linewidth = 0.43)+
  geom_hline(yintercept = threshold_t1, color= "red", linewidth = 0.43) +
  
  # 添加箭头
  #annotate("segment", x = arrow_x, xend = arrow_x, y =arrow_y_end, yend = arrow_y_start, 
  #arrow = arrow(length = unit(0.3, "cm")), color = "black")+
  
  
  
  theme_bw() + # theme(legend.position = 'none')去除图例
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.55,size=9, color="black", hjust = 0.55), 
        panel.border = element_rect(color = "black", size=0.43),
        text=element_text(size=9, color="black"))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.43))+ # 调整x轴
  theme(axis.ticks.x=element_line(color="black",size=0.43,lineend = 1))+  # 调整x轴刻度粗细
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.43))+ # 调整y轴
  theme(axis.ticks.y=element_line(color="black",size=0.43,lineend = 10)) + # 调整y轴刻度线
  theme(panel.border = element_blank(), #去除坐标轴线默认填充的灰色
        axis.text.x = element_text(angle = 0,
                                   vjust = 0.55,
                                   size = 9, 
                                   colour = "black"),
        axis.text.y = element_text(vjust=0.55, 
                                   size=9, 
                                   color="black", 
                                   hjust=0.55),
        axis.line = element_line(colour = "black",size=0.2, lineend = "square"), #将x=0轴和y=0轴加粗显示(size=0.3), x、y交界处是完美交接(lineend = "square"语句
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))+
  theme(legend.position = "right") +
  theme(legend.key.size = unit(30, "pt")) +  # 图例大小
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title=element_blank()) +  #  隐藏图例标题
  #scale_y_continuous(limits = c(0, 0.08), guide = "prism_offset")
  labs(x="Position(Mb)", y="XP-CLR") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(expand = c(0, 0)) #scale_y_continuous(expand = c(0, 0)) 是一个用于控制 y 轴扩展范围的参数设置。在 ggplot2 中，expand 参数的作用是为绘图区域的轴添加一些额外的空间，避免图形元素紧贴轴线或绘图边界

print(p9)
ggsave("Seita.2G444300_c1_c2_50k_xpclr.pdf", plot = p9, width = 6.5, height = 3)

