# 1 下载表达矩阵并做最基本的预处理

# 对应代码：https://github.com/jmzeng1314/GEO/blob/master/GSE42872_main/step1-download.R

f='GSE42872_eSet.Rdata'
#在R语言中，eSet是一个用于存储基因表达数据的类。它是Bioconductor为基因表达数据格式所定制的标准之一，也是所有涉及基因表达量相关数据在Bioconductor中进行操作的基础数据类型之一1.

library(GEOquery)

### 这里就是从零开始的起点，找到目标GSE后的第一步：

## 这次演示的是代码下载




################## 方法一：
if(!file.exists(f)){
  gset <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,     ## 不下注释文件
                 getGPL = F)       ## 不下平台文件
  save(gset,file=f)   ## 保存到本地
}

##注意：这里得到的gset是一个对象！
#里面包含了很多东西，不止表达矩阵。各种杂七杂八的信息都在里面

#### 这两个GPL为F之后，其实就下载了一个gz文件。
#你提供的代码是用R编程语言编写的。getGEO()函数用于从NCBI的GEO数据库下载和解析GEO（基因表达式综合数据库）数据1。
#GSE42872是你正在尝试下载的数据集的GEO系列访问号。destdir参数指定任何下载的目标目录。它默认为与体系结构相关的tempdir。如果您想将文件保存以供以后使用，您可能需要指定不同的目录1。
#AnnotGPL参数指定是否下载实验中使用的平台的！注释！数据。默认值为TRUE，但您已将其设置为FALSE1。
#getGPL参数指定是否下载！实验的平台数据！。默认值为TRUE，但您已将其设置为FALSE1。

### 在GEOoverview中可以找到官方描述：
# 平台记录由阵列或序列器的摘要描述组成，对于基于阵列的平台， 定义数组模板的数据表。每个平台记录都分配了一个唯一且稳定的记录 GEO 入藏号（GPLxxx）。一个平台可以引用许多示例 已由多个提交者提交。
# 样本记录描述了处理单个样本的条件， 它所经历的操纵，以及每个操作的丰度测量 元素派生自它。每个示例记录都分配有唯一的和 稳定的GEO入藏号（GSMxxx）。示例实体必须引用 只有一个平台，可以包含在多个系列中。
#  系列记录将一组相关样本链接在一起，并提供整个研究的焦点和描述。 系列记录还可能包含描述提取数据的表， 总结结论或分析。每个系列记录都分配有一个 唯一且稳定的GEO入藏号（GSExxx）。
# 在生物信息里，GSE代表的是GEO Series，GPL代表的是GEO Platform，而GSM代表的是GEO Sample。

b = gset[[1]]
#注意！！！！b在这里仍然是一个对象！
#gset[[1]]是一个列表中的第一个元素。在你的代码中，你已经使用了getGEO()函数下载了GEO数据集，并将其存储在变量gset中。这个变量是一个列表，其中包含有关数据集的信息。通过使用双方括号（[[）和数字1，您正在提取列表中的第一个元素12。
#如果要从b中再提取出表达矩阵，那要如下操作：
#exprs的目标只能是object
raw_exprSet = exprs(b)
#在R语言中，exprs()是一个函数，用于从一个ExpressionSet对象中提取表达矩阵。在你的代码中，raw_exprSet = exprs(b)的意思是将b中的表达矩阵提取出来并将其存储在名为raw_exprSet的变量中1。
#在R语言中，exprs()函数从一个ExpressionSet对象中提取表达矩阵。如果b是一个ExpressionSet对象，那么exprs(b)将返回一个数值矩阵，其中行名是基因，列名是细胞编号1。




###################方法二：
### 如果使用网页下载txt.gz文件，可以用如下代码打开tz压缩包中的表达文件：
a = read.table('GSE42872_series_matrix.txt.gz',
               sep = '\t',quote = '',fill = T,comment.char = '!',header = T)
#这段代码是R语言中的代码，其中a是一个数据框，read.table()函数用于读取文本文件，
#'GSE42872_series_matrix.txt.gz’是文件名，sep = '\t’表示使用制表符作为分隔符，quote = ''表示不使用引号，fill = T表示填充缺失值，comment.char = '!'表示使用感叹号作为注释符。header = T表示第一行是列名。
### 注意，这里的a就不是对象了，而是直接读取的数据框。
## 下面要对a进行进一步的处理，这也是比方法一多出来的步骤：
rownames(a) = a[,1]
a = a[,-1]#去掉第一列
## 这两步之后，a就和上一方法得到的raw_exprSet的一样了。（除了列名）
## 而上面的两步，也就是getgeo这个函数在做的事情，它也是在做预处理。这里可以用？看下官方文档。


#视频额外写了一个函数，把预处理好的矩阵转换为csv格式可以被excel打开，方便交流。
downGSE <- function(studyID = 'GSE1009',destdir = '.'){
  
  library(GEOquery)
  eSet <- getGEO(studyID,destdir = destdir,getGPL = F)
  
  exprSet = exprs(eSet[[1]])
  pData = pData(eSet[[1]])
  
  write.csv(exprSet,paste0(studyID,'exprSet.csv'))
  write.csv(pData,paste0(studyID,'_metadata.csv'))
}
#这段代码是R语言中的代码，其中downGSE是一个函数，studyID = 'GSE1009’表示数据集的ID，destdir = '.'表示下载文件的目录为当前目录。
#这个函数使用GEOquery包中的getGEO()函数下载数据集，然后将表达矩阵和样本信息分别保存为exprSet.csv和_metadata.csv文件。
downGSE('GSE42872')




##############################
##############################
## 带着简单下载后预处理好的表达矩阵，进入ID转换。##


load('GSE42872_eSet.Rdata')  ## 载入数据


group_list = c(rep('control',3),rep('case',3)) 
#在R语言中，c()函数用于将一组值组合成向量。在你的代码中，group_list = c(rep('control',3),rep('case',3))的意思是创建一个名为group_list的向量，其中包含3个值为’control’和3个值为’case’的元素。

save(raw_exprSet,group_list,
     file = 'GSE_raw_exprSet.Rdata')
### 重点：这个文件就要给第二步使用了！！！！！！！！
#这是R语言中的一种保存数据的方法。其中，raw_exprSet是一个表达矩阵，group_list是一个分组信息的列表。这个函数将这两个对象一起保存到名为’GSE_raw_exprSet.Rdata’的文件中。

### 中间讲解

































