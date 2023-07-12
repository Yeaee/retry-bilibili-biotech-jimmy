## 这一步ID转换就是要把探针的数字编号转换为熟知的基因名称。
## 拿取上一步得到的原始矩阵，进行这一环节。



#################前置准备：

#首先要加载一些特定的三方数据库包，这些包里的数据库是基础素材。
#这里要去找原数据集使用芯片平台对应的包。
#搜索GSE42872在首页找到platform，找到其对应为GPL6244.
#edge搜索“GPL6244 r package”，找到bioconductor的官网，它会直接给出下载代码。
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pd.hugene.1.0.st.v1")
## 但是官网给出的好像不全，这里就可以去生信的博客右上角搜索框里搜索：http://www.bio-info-trainee.com/
## 找到对应关系如下：67  GPL6244       Homo sapiens   hugene10sttranscriptcluster
## 安装导入，记得加.db,博客给出的只是前缀。
BiocManager::install("hugene10sttranscriptcluster.db")



################## 开始寻找对应关系：
library(hugene10sttranscriptcluster.db)






