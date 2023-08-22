#function

#####################
#function1：GSEnumber_to_visualized_raw_exprSet
#输入：GSE号，纯GSE+数字
#输出：return该GSE的原始表达集raw_exprSet。
#输出：本地save下txt_gz_file_eSet.Rdata和GSE_raw_exprSet文件。
GSEnumber_to_visualize_raw_exprSet <- function(GSE_number){
  library(GEOquery)
  GSE_txt_gz_flie <- getGEO(GSE_number,destdir = '.',
                            AnnotGPL = F,
                            getGPL = F)
  save(GSE_txt_gz_flie,file=paste(GSE_number,"_eSet.Rdata"))
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
#function3：raw_exprSet_to_final_exprSet
#输入：raw_exprSet
#输出：final_exprSet
raw_exprSet_to_final_exprSet <- function(raw_exprSet){
  
}




















