###  先打开网址:http://www.bio-info-trainee.com/bioconductor_China/software/limma.html

###  前置准备：得到待处理的表达矩阵
suppressPackageStartupMessages(library(CLL))
data(sCLLex)
exprSet=exprs(sCLLex)   ##sCLLex是依赖于CLL这个package的一个对象
samples=sampleNames(sCLLex)
pdata=pData(sCLLex)
group_list=as.character(pdata[,2])
dim(exprSet)

###  查一下矩阵是不是一样，注意！这里视频里的跟这个不一样，progress和stable是对的
exprSet
dim(exprSet)
group_list

### 设计矩阵
library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = colnames(exprSet)
design

### 代码手工复现矩阵
tmp = data.frame(case = c(0,0,0,1,1,1),control=c(1,1,1,0,0,0))

###  做比较矩阵
contrast.matrix<-makeContrasts(paste0(unique(group_list,collapse = "-")),levels = design)
contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较

### 表达 矩阵做完后做后续的分析
##step1
fit <- lmFit(exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)


### 浅画个热图
## heatmap 
library(pheatmap)
# 选取头部的25个差异最显著的基因，做一个小表达矩阵
choose_gene=head(rownames(nrDEG),25)
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix)

###最后，差异分析很简单，只要搞懂input和output。
