#function

###################
#function1：GSEnumber_to_visualized_raw_exprSet
#输入：GSE号，纯GSE+数字
#输出：return该GSE的原始表达集raw_exprSet。
#输出：本地save下txt_gz_file_eSet.Rdata和GSE_raw_exprSet文件。
GSEnumber_to_visualize_raw_exprSet <- function(GSE_number){
  library(GEOquery)
  GSE_txt_gz_flie <- getGEO(GSE_number,destdir = '.',
                            AnnotGPL = F,
                            getGPL = F)
  save(GSE_txt_gz_flie,file=paste0(GSE_number,"_eSet.Rdata"))
  group_list = c(rep('control',3),rep('case',3)) 
  
  the_first_object <- GSE_txt_gz_flie[[1]]
  new_name = array(paste0(GSE_number,'_raw_exprSet'))
  raw_exprSet <- exprs(the_first_object)
  save(raw_exprSet,group_list,
       file = paste0(GSE_number,'_raw_exprSet.Rdata'))
  return(raw_exprSet)}
#使用案例：
GSE42872_raw_exprSet<-GSEnumber_to_visualize_raw_exprSet('GSE42872')


###################
#function2:GSEnumber_to_visualized_excel
#输入：GSE号，纯GSE+数字
#输出：默认当前文件夹下该GSE的可视化excel表格
GSEnumber_to_visualized_excel <- function(studyID,destdir = '.'){
  
  library(GEOquery)
  eSet <- getGEO(studyID,destdir = destdir,getGPL = F)
  
  exprSet = exprs(eSet[[1]])
  pData = pData(eSet[[1]])
  
  write.csv(exprSet,paste0(studyID,'_exprSet.csv'))
  write.csv(pData,paste0(studyID,'_metadata.csv'))
}
#使用案例（在该R脚本的默认文件夹下可以看到csv文件）：
GSEnumber_to_visualized_excel('GSE42872')



###################
#function3:make_GSEs_ids
#人工：搜索GSE号在GEO首页找到platform一栏，找到其对应platform为GPL****或其他.
#人工: 在网站http://www.bio-info-trainee.com/1399.html寻找该GPL对应的bioc_package
#输入：该GSE_number芯片平台对应的bioc_package的.db与SYMBOL形式。
#输出：制作好的ids
make_GSEs_ids <- function(symbol){
  ####if (!require("BiocManager", quietly = TRUE))
  ####  install.packages("BiocManager")
  #library(BiocManager)
  #BiocManager::install(array(paste0(bioc_package,".db")))
  #这里报错,怎么搞都导不进去，算了，这步就手动挡导入库吧，放在最前面就行
  #library(db)
  ids = toTable(symbol)
  return(ids)
}
#使用案例：
##搜索GSE号在GEO首页找到platform一栏，找到其对应platform为GPL****或其他.
##在网站http://www.bio-info-trainee.com/1399.html
##寻找该GPL对应的bioc_package为hugene10sttranscriptcluster；
##手动输入.db格式和SYMBOL格式。
BiocManager::install("hugene10sttranscriptcluster.db")
library("hugene10sttranscriptcluster.db")
GSE_42872_ids <- make_GSEs_ids(hugene10sttranscriptclusterSYMBOL)



###################
#function4：raw_exprSet_to_final_exprSet
#输入：function1得到的raw_exprSet以及function3得到的ids。
#输出：final_exprSet
#函数简介：这个函数滤了两次,根据symbol分了个类,最后把数字换成了基因名。
raw_exprSet_to_final_exprSet <- function(exprSet,ids){
  table(rownames(exprSet)%in%ids$probe_id)
  exprSet = exprSet[rownames(exprSet)%in%ids$probe_id,]
  ids = ids[match(rownames(exprSet),ids$probe_id),]
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rownames(x))]) 
  probes = as.character(tmp)
  new_exprSet = exprSet[rownames(exprSet)%in% probes,]
  rownames(new_exprSet)=ids[match(rownames(new_exprSet),ids$probe_id),2]
  return(new_exprSet)
}
#使用范例：
GSE42872_final_exprSet <- raw_exprSet_to_final_exprSet(GSE42872_raw_exprSet,GSE_42872_ids)


###################
#function5:make_group_list_of_GSE_exprSet
#人工：去该GSE主页看下sample一栏中的样本情况。
#输入：控制组个数，对照组个数。
#输出：根据样本情况写出group_list。
make_group_list_of_GSE_exprSet <- function(number_of_control,number_of_case){
  group_list = c(rep('control',number_of_control),rep('case',number_of_case))
  return(group_list)
}
#使用案例：
GSE_42872_group_list <- make_group_list_of_GSE_exprSet(3,3)


###################
#function6:check_the_analysis_value_of_final_exprSets
#输入：GSE*****_final_exprSet,GSE_*****_group_list
#输出：箱线图
#解释：就是画箱线图（检测图里的线是不是差不多在一条线上，是则代表矩阵有分析价值）
check_the_analysis_value_of_final_exprSets <- function(new_exprSet,group_list){
  library(reshape2)
  library(ggplot2)
  exprSet_L=melt(new_exprSet)
  colnames(exprSet_L)=c('probe','sample','value')
  exprSet_L$group=rep(group_list,each=nrow(new_exprSet))
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank(),panel.background = element_rect(fill = 'lightgrey'))
  print(p)
}
#案例：
check_the_analysis_value_of_final_exprSets(GSE42872_final_exprSet,GSE_42872_group_list)



###################
#function7:check_if_difference_exsits_between_control_and_case
#输入：GSE*****_final_exprSet,GSE_*****_group_list,sample_size(为control加case数)
#输出：hclust_plot，一种树状图，显示数据的分层结构和聚类关系。
#解释：看图里的case和control分的开不开，分的开就是好，说明存在差异。
check_if_difference_exsits_between_control_and_case <- function(new_exprSet,group_list,sample_size){
  colnames(new_exprSet)=paste(group_list,1:sample_size,sep='')
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                  cex = 0.7, col = "blue")
  hc=hclust(dist(t(new_exprSet)))
  par(mar=c(5,5,5,10)) 
  plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)
}
#案例：
check_if_difference_exsits_between_control_and_case(GSE42872_final_exprSet,GSE_42872_group_list,6)



###################
#function8:make_sure_the_final_exprSet_can_be_used
#输入：GSE*****_final_exprSet,GSE_*****_group_list。
#输出：PCA_plot，用来绘制主成分分析（PCA）的结果，也就是一种散点图。
#解释：PCA图：也是看分地开不开，三个图已经够了，可以确定矩阵没什么问题。
make_sure_the_final_exprSet_can_be_used <- function(new_exprSet,group_list){
  library(ggfortify)
  df=as.data.frame(t(new_exprSet))
  df$group=group_list 
  autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
}
#案例：
make_sure_the_final_exprSet_can_be_used(GSE42872_final_exprSet,GSE_42872_group_list)
