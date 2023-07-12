# 0 项目总览
# 对应代码：https://github.com/jmzeng1314/GEO/blob/master/GSE42872_main/step1-download.R

f='GSE42872_eSet.Rdata'
#在R语言中，eSet是一个用于存储基因表达数据的类。它是Bioconductor为基因表达数据格式所定制的标准之一，也是所有涉及基因表达量相关数据在Bioconductor中进行操作的基础数据类型之一1.

library(GEOquery)

if(!file.exists(f)){
  gset <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
#### 这两个GPL为F之后，其实就下载了一个gz文件。
#你提供的代码是用R编程语言编写的。getGEO()函数用于从NCBI的GEO数据库下载和解析GEO（基因表达式综合数据库）数据1。
#GSE42872是你正在尝试下载的数据集的GEO系列访问号1。destdir参数指定任何下载的目标目录。它默认为与体系结构相关的tempdir。如果您想将文件保存以供以后使用，您可能需要指定不同的目录1。
#AnnotGPL参数指定是否下载实验中使用的平台的注释数据。默认值为TRUE，但您已将其设置为FALSE1。
#getGPL参数指定是否下载实验的平台数据。默认值为TRUE，但您已将其设置为FALSE1。

load('GSE42872_eSet.Rdata')  ## 载入数据

b = gset[[1]]
#gset[[1]]是一个列表中的第一个元素。在你的代码中，你已经使用了getGEO()函数下载了GEO数据集，并将其存储在变量gset中。这个变量是一个列表，其中包含有关数据集的信息。通过使用双方括号（[[）和数字1，您正在提取列表中的第一个元素12。

raw_exprSet = exprs(b)
#在R语言中，exprs()是一个函数，用于从一个ExpressionSet对象中提取表达矩阵。在你的代码中，raw_exprSet = exprs(b)的意思是将b中的表达矩阵提取出来并将其存储在名为raw_exprSet的变量中1。
#在R语言中，exprs()函数从一个ExpressionSet对象中提取表达矩阵。如果b是一个ExpressionSet对象，那么exprs(b)将返回一个数值矩阵，其中行名是基因，列名是细胞编号1。

group_list = c(rep('control',3),rep('case',3))
#在R语言中，c()函数用于将一组值组合成向量。在你的代码中，group_list = c(rep('control',3),rep('case',3))的意思是创建一个名为group_list的向量，其中包含3个值为’control’和3个值为’case’的元素。



















