## ID转换

## 安装第三方包（ID数据库）
BiocManager::install('hugene10sttranscriptcluster.db')

##  下载数据aaa
library(GEOquery)
gset <- getGEO('GSE42872', destdir=".")


##  去官网下载gz后缀文件，放入与Rproject同一文件夹下，大概几百kb


## 读取数据
a = read.table('GSE42872_series_matrix.txt.gz',
               sep = '\t',quote = "",fill = T,comment.char = '!',
               header = T)

##  获取a1
b = gset[[1]]
a1 = exprs(b)


##  把a的行名称改为ID号：
rowname(a) = a[,1]
a = a[,-1]


## 开始导入包
library(hugene10sttranscriptcluster.db)


## 找对应关系 
ids=toTable(hugene10sttranscriptclusterSYMBOL)
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))

##  查看不同探针个数下对应的基因数量
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))


##  看下哪些基因不在probr_id里
table(rownames(exprSet) %in% ids$probe_id)

## 查看探针数
dim(exprSet)

##  做一个过滤，过滤掉不在里面的
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]

##  再次查看
dim(exprSet)

##  改探针顺序为了和表达矩阵顺序一样
ids=ids[match(rownames(exprSet),ids$probe_id),]

## 核查是否一样
head(ids)
exprSet[1:5,1:5]


##通过ID的symbol对表达矩阵进行分类
tmp = by(exprSet,
         ids$symbol,
          function(x) rownames(x)[which.max(rowMeans(x))] )


##  再把tmp做成列表
probes = as.character(tmp)


##  老三样，这次是去除重复id
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)


##  提取目的探针(IGKC)的表达矩阵
x = exprSet[ids[,2]=='IGKC',]
