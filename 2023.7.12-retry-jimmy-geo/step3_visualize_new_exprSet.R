## 带着上一step处理好的new_exprSet,这一步对之前处理得到的new_exprSet进行检验：用可视化的图像观察矩阵是否合理。
## 检测处理好的new_exprSet是否有处理的意义：
## 查看样本分布是否合理，用PCA图，hclust图

###########！！视频里一开始用的是另一个表达矩阵。
## 先把要用的表达矩阵load进来,如果有上一步带进来的new_exprSet可以跳过：
load("GSE_raw_exprSet.Rdata")
b = gset[[1]]
#这里跟第一步一样，b是raw_exprSet,最后得到的就是new_exprSet.

######准备工作完成正式开始：

## 第一步：自己制作group_list:
###############！！！含金量的步骤：要先自己制作表达矩阵的group_list.
## 首先去官网gse42872的界面，看一下sample的分布；为三个control，三个list。
## 或者可以在R里直接调用pData来查看分组信息,注意pData的括号内为对象：
view_group = pData(b)
# 第一列就是分组情况
group_list = c(rep('control',3),rep('case',3))
# 制作分组完成。


## 第二步；有了new_exprSet和group_list，就可以进行下面的分析：
###################################
###################################下面的内容同样可以进行函数封装，放到后面。
###################################
library(reshape2)
exprSet_L=melt(new_exprSet)
colnames(exprSet_L)=c('probe','sample','value')
exprSet_L$group=rep(group_list,each=nrow(new_exprSet))
head(exprSet_L)
### ggplot2 
library(ggplot2)
###  开始画各种图
p=ggplot(exprSet_L,
         aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)
###  这下面一堆图跟上面的图都是一个意思，形式不一样，挑最后这个最好看的做成函数。
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

## hclust 
colnames(exprSet)=paste(group_list,1:22,sep='')
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
hc=hclust(dist(t(exprSet)))
par(mar=c(5,5,5,10)) 
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)


###################################
###################################
###################################

#####下面进行上述部分代码的函数封装：
#为了避免重复导入库，先把要用的库导进去：
library(ggplot2,reshape2)

#首先是画箱线图：（检测图里的线是不是差不多在一条线上，是代表矩阵有分析价值）
view_matrix_box_plot <- function(new_exprSet,group_list){
  exprSet_L=melt(new_exprSet)
  colnames(exprSet_L)=c('probe','sample','value')
  exprSet_L$group=rep(group_list,each=nrow(new_exprSet))
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  print(p)
}
view_matrix_box_plot(new_exprSet,group_list)
#制作完成。



#接着画hclust图：
#注意这里多了一个参数，样本数量
view_hclust_plot <- function(new_exprSet,group_list,sample_size){
  colnames(new_exprSet)=paste(group_list,1:sample_size,sep='')
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                  cex = 0.7, col = "blue")
  hc=hclust(dist(t(new_exprSet)))
  par(mar=c(5,5,5,10)) 
  plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)
}
view_hclust_plot(new_exprSet,group_list,6)
#制作完成。













