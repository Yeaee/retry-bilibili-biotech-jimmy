##  了解表达矩阵

##  内参即内部参照，一般是指由管家基因编码表达的蛋白，它们在各组织和细胞中的表达相对恒定，在检测蛋白的表达水平变化时常用它来做参照物，类似marker

## 首先要安装好必要的包，要哪个安装哪个，注意，安装的时候包名称区分大小写，如下：
BiocManager::install('CLL')
BiocManager::install('hgu95av2.db')

## 接着跑这个代码，得到作为本次材料的表达矩阵（到67行）

if(T){
  # BiocManager::install('CLL')
  suppressPackageStartupMessages(library(CLL))
  data(sCLLex)
  sCLLex
  exprSet=exprs(sCLLex)   
  ##sCLLex是依赖于CLL这个package的一个对象
  samples=sampleNames(sCLLex)
  pdata=pData(sCLLex)
  group_list=as.character(pdata[,2])
  dim(exprSet)
  exprSet[1:5,1:5]
  
  # BiocManager::install('hgu95av2.db')
  library(hgu95av2.db)
  ls("package:hgu95av2.db") 
  
  ids=toTable(hgu95av2SYMBOL)
  save(ids,exprSet,pdata,file = 'input.Rdata')
  length(unique(ids$symbol))
  tail(sort(table(ids$symbol)))
  table(sort(table(ids$symbol)))
  
  plot(table(sort(table(ids$symbol))))
  
  table(rownames(exprSet) %in% ids$probe_id)
  dim(exprSet)
  exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
  dim(exprSet)
  
  ids=ids[match(rownames(exprSet),ids$probe_id),]
  head(ids)
  exprSet[1:5,1:5]
  if(F){
    tmp = by(exprSet,ids$symbol,
             function(x) rownames(x)[which.max(rowMeans(x))] )
    probes = as.character(tmp)
    dim(exprSet)
    exprSet=exprSet[rownames(exprSet) %in% probes ,]
    dim(exprSet)
    
    rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
    exprSet[1:5,1:5]
  }
  
  
  identical(ids$probe_id,rownames(exprSet))
  dat=exprSet
  ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
  dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
  rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
  dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
  dim(dat)
  
}


## 进行操作，把数字的行名转换成基因名称：
exprSet=dat

## 查看基因表达矩阵的分布图
library(reshape2)
exprSet_L=melt(exprSet)
colnames(exprSet_L)=c('probe','sample','value')
exprSet_L$group=rep(group_list,each=nrow(exprSet))
head(exprSet_L)
### ggplot2 
library(ggplot2)

###  开始画各种图
p=ggplot(exprSet_L,
         aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)
###  这下面一堆图跟上面的图都是一个意思，形式不一样。
p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_violin()
print(p)
p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
print(p)
p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
print(p)
p=ggplot(exprSet_L,aes(value,col=group))+geom_density() 
print(p)
p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
p=p+theme_set(theme_set(theme_bw(base_size=20)))
p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
print(p)

### 画pca图
plot(hclust(dist(t(exprSet))))

### 改col names
colnames(exprSet) = paste(group_list,1:22,sep='')

###再次plot
plot(hclust(dist(t(exprSet))))

### 画PCA图

# BiocManager::install('ggfortify')
library(ggfortify)
df=as.data.frame(t(exprSet))
df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
