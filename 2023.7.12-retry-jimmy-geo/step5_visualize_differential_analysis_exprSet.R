#接着上一步拿到的nrDEG差异分析矩阵，现在对这个矩阵进行可视化分析：


# 先把nrDEG打开看一眼，熟悉一下一个常规的差异分析表达矩阵都有哪些元素。
##### 贴一段对差异分析矩阵中列名的分析：
# logFC (log Fold Change)：这是差异表达分析中一个非常重要的参数。它表示在两个条件或分组（如疾病组和对照组）之间的基因表达量的比值（fold change）的对数。取对数的目的主要是为了处理数据的偏度和使数据更接近正态分布。如果logFC > 0，表示在实验组相比对照组基因上调；如果logFC < 0，表示基因在实验组相比对照组下调。
#p-value (P值)：这是检验基因表达差异是否显著的一个统计参数。一般来说，P值越小，表示这个基因在两个组间的表达差异越显著。
#adj.P.Val (校正后P值或FDR)：在多重比较中，需要对P值进行校正以防止假阳性结果的出现。校正后的P值又被称为FDR（False Discovery Rate，假发现率）。这个值越小，说明差异表达更有可能是真实的。
#aveexpr：这个参数通常代表“平均表达水平”(Average Expression Level)，用来描述一个基因在所有样本中的平均表达水平。对于微阵列数据，它通常是所有样本的原始强度值的对数转换后的平均值。对于RNA-seq数据，它可能是所有样本的原始计数的对数转换后的平均值。
#t：这个参数代表的是t统计值。在生物信息学的上下文中，这通常用于描述基因表达量在不同实验组之间的差异的统计显著性。t统计值基于样本均值、标准差以及样本大小来计算，如果t值较大，那么两个样本组之间的基因表达差异就更有可能是真实存在的，而非偶然产生的。
#B：这个参数可能表示贝叶斯统计值，但它的具体意义可能取决于你正在使用的具体生物信息学工具。在一些工具中，例如limma包，B值被用来描述基因在不同实验组之间表达差异的后验概率，也就是说，给定观测到的数据，一个基因在两个实验组之间有显著差异的概率。B值较大意味着我们对存在显著差异的信念更强。
#以上几个参数通常被用来筛选和鉴定差异表达基因。例如，研究者们可能会选择logFC大于2（或小于-2），并且FDR小于0.05的基因作为显著差异表达的基因。具体的阈值会依据研究背景和问题而变化。
library(ggplot2)

#################火山图
###  初次尝试：
plot(nrDEG$logFC,-log10(nrDEG$P.Value))

###  画图升级版，直接做成函数：
draw_volcano_plot<-function(DEG){
  logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
  DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
  )
  this_tile
  head(DEG)
  g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
  print(g)
}
draw_volcano_plot(nrDEG)

#################热图
### 浅画个热图,同样做成函数
## heatmap 
draw_top25gene_pheatmap <- function(new_exprSet,nrDEG){
  library(pheatmap)
  # 选取头部的25个差异最显著的基因，做一个小表达矩阵
  choose_gene=head(rownames(nrDEG),25)
  choose_matrix=new_exprSet[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  pheatmap(choose_matrix)
  #颜色深浅就是表达量的高低
}
draw_top25gene_pheatmap(new_exprSet,nrDEG)


