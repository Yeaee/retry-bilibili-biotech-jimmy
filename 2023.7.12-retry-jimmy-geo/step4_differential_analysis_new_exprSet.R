## 带着step2得到的new_exprSet,和step3做好的group_list进入差异分析。
## 在这一步的差异分析里，重要的还是制作两个新矩阵：分组矩阵和比较矩阵。

####第一步制作design矩阵：
design <- model.matrix(~0+factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = colnames(exprSet)
design


####老惯例，做成函数：
make_design_matrix <- function(new_exprSet,group_list){
  design <- model.matrix(~0+factor(group_list))
  colnames(design) = levels(factor(group_list))
  rownames(design) = colnames(new_exprSet)
  return(design)
}
design = make_design_matrix(new_exprSet,group_list)
####制作完成,得到design矩阵

##########################################################

###第二步制作比较矩阵：
contrast.matrix<-makeContrasts(paste0(unique(group_list,collapse = "-")),levels = design)
contrast.matrix ##这个矩阵声明，我们要把两组进行差异分析比较

####做成函数，输入要group_list和上一步的design矩阵：
make_contrast_matrix <- function(group_list,design){
  library(limma)
  contrast.matrix<-makeContrasts(paste0(unique(group_list,collapse = "-")),levels = design)
  return(contrast.matrix)
}
contrast_matrix = make_contrast_matrix(group_list,design)
#制作完成，得到比较矩阵

###########################################################

#第三步就可以把三个矩阵喂进去，得到差异分析结果了：
#这里直接做成函数：
differential_analysis<- function(new_exprSet,design,contrast_matrix){
  library(limma)
  fit <- lmFit(new_exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast_matrix) ##这一步很重要，大家可以自行看看效果
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  return(nrDEG)
}
nrDEG = differential_analysis(new_exprSet,design,contrast_matrix)
#制作完成。


#########下一步就是把结果可视化了，带着差异分析矩阵nrDEG进入下一步。











