### 首先要搞懂富集分析要喂进去什么：
### 要喂进去火山图里得到的400多加400多差异显著的基因片段，对这些基因片段做富集分析，即注释。
### 要喂的格式不是基因名，是对应的entrezID
### 富集分析的目的，找到这个药物特殊作用哪几个通路。

# 先跟着视频走吧，这里有点难的。
library(cluster)
library(org.Hs.eg.db)
library(DOSE)
library(clusterProfiler)
#先随便取点基因
gene = head(rownames(nrDEG),1000)
# 怎么用R取到那几百个差异表达的基因是第一个点。

#拿到了之后，还要做转换时第二个点。
data(geneList, package="DOSE")# 这个genelist查看一下，就是一个基因名对应它的logFC，没了。
gene <- names(geneList)[abs(geneList) > 2]# 我懂了，rownames不止是字符串，后面的值还在。
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
# 已经处理好了。解释一下这段代码。
# 这段代码的意思是：data(geneList, package=“DOSE”)是载入DOSE包中的数据。
# gene <- names(geneList)[abs(geneList) > 2]是将基因名存储在gene中，
# 其中abs(geneList) > 2是筛选出foldchange绝对值大于2的基因。
# gene.df <- bitr(gene, fromType = “ENTREZID”, toType = c(“ENSEMBL”, “SYMBOL”), OrgDb = org.Hs.eg.db)是将基因名从ENTREZID转换为ENSEMBL和SYMBOL。head(gene.df)是查看转换后的前几行数据。
# 这个步骤好像已经把挑选也做好了。

# 首先要喂的应该是entrezID形式的基因名,用bitr函数来做。

# 开始进行KEGG通路分析
kk <- enrichKEGG(gene = gene,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)
# enrichKEGG(gene = gene, organism = 'hsa', pvalueCutoff = 0.05) 是一个 R 语言的函数，用于对基因集进行 KEGG 通路富集分析。其中，gene 是一个包含基因 ID 的向量，organism 是指定物种，pvalueCutoff 是指定显著性水平的阈值.
# hsa 是指人类（Homo sapiens）的缩写，是 KEGG 数据库中人类基因组的代号。
head(kk)[,1:6]
# 分析结果也能出来