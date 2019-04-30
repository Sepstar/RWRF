####20180111byhes，本脚本目的是根据前期初筛结果，对重点关注的几种肿瘤的融合聚类效果复现。

#source("E:/博士/随机游走融合/RWR_fusion/RWR_fusion_neighbor.R")
#source("E:/博士/随机游走融合/RWR_fusion/RWR_fusion.R")
source("/home/wenyuqi/RWR/RWR_fusion_neighbor_new.R")
source("/home/wenyuqi/RWR/RWR_fusion_new.R")
.libPaths("/home/wenyuqi/Rpackage/")# 不需要是Rpackage-R-3.3.3
library(SNFtool)
library(survival)
library(clValid)
library(kernlab)
library(foreach)
library(doParallel)
cl<-makeCluster(8)  
registerDoParallel(cl) 
clusterEvalQ(cl, .libPaths("/home/wenyuqi/Rpackage/"))
# setwd("/home/huangxin/yujijun/originaldata2")
setwd("/home/wenyuqi/RWR/robust/")
set.seed(7)

####2ACC设置参数####
subpath="2ACC"
CLUSTER_NUM=3
gama_list=seq(from=0.1,to=0.9,by=0.1)
neighbor_num_list=seq(from=5,to=30,by=5)
alpha_list=seq(from=0.1,to=0.9,by=0.1)
beta_list=seq(from=0.1,to=0.9,by=0.1)
####2ACC####
mrna_data=as.matrix(read.table(file = paste("/home/100/Cancer_originaldata2/",subpath,"/HiSeqV2_PANCAN2.txt",sep = ""),header = T,sep = "\t",quote = "",check.names = F))
mirna_data=as.matrix(read.table(file = paste("/home/100/Cancer_originaldata2/",subpath,"/miRNA_HiSeq_gene1.txt",sep = ""),header = T,sep = "\t",quote = "",check.names = F))
meth_data=as.matrix(read.table(file = paste("/home/100/Cancer_originaldata2/",subpath,"/HumanMethylation4501.txt",sep = ""),header = T,sep = "\t",quote = "",check.names = F))
mrna_data=t(mrna_data)
mirna_data=t(mirna_data)
meth_data=t(meth_data)
dataL=list(mrna_data)
dataL=c(dataL,list(mirna_data))
dataL=c(dataL,list(meth_data))

survival_data=read.table(file = paste("/home/100/Cancer_originaldata2/",subpath,"/clinicalMatrix15.txt",sep = ""),header = TRUE, sep="\t")
survival_data[,1]=survival_data[,1]/30
survival_data=cbind(survival_data,survival_data[,2])
colnames(survival_data)[3]="Cluster"

##利用SNFtool包构建相似性矩阵
distL = lapply(dataL, function(x) dist2(x, x))
affinityL = lapply(distL, function(x) affinityMatrix(x, 20, 0.5))

####gama####
gama_result_matrix=matrix(data = 0,nrow = 9,ncol = 4)
colnames(gama_result_matrix)=c("RWR-p-value","RWR-Dunn","RWRN-p-value","RWRN-Dunn")
rownames(gama_result_matrix)=gama_list
RWR_similarity_fusion_list=NULL
RWRN_similarity_fusion_list=NULL
RWR_similarity_fusion_temp_list=NULL
RWRN_similarity_fusion_temp_list=NULL
for (gama_index in 1:length(gama_list))
{
  ##融合
  time1=Sys.time()
  RWR_similarity_fusion=RWR_fusion(sim_list = affinityL,gama=gama_list[gama_index])
  time2=Sys.time()
  print(time2-time1)
  RWR_similarity_fusion_temp=RWR_similarity_fusion
  diag(RWR_similarity_fusion)=0
  
  time1=Sys.time()
  RWRN_similarity_fusion=RWR_fusion_neighbor(sim_list = affinityL,gama=gama_list[gama_index])
  time2=Sys.time()
  print(time2-time1)
  RWRN_similarity_fusion_temp=RWRN_similarity_fusion
  diag(RWRN_similarity_fusion)=0
  RWR_similarity_fusion_list=c(RWR_similarity_fusion_list,list(RWR_similarity_fusion))
  RWRN_similarity_fusion_list=c(RWRN_similarity_fusion_list,list(RWRN_similarity_fusion))
  RWR_similarity_fusion_temp_list=c(RWR_similarity_fusion_temp_list,list(RWR_similarity_fusion_temp))
  RWRN_similarity_fusion_temp_list=c(RWRN_similarity_fusion_temp_list,list(RWRN_similarity_fusion_temp))
  
}

for (gama_index in 1:length(gama_list))
{
  ##聚类
  RWR_result=specc(RWR_similarity_fusion_list[[gama_index]],centers=CLUSTER_NUM)
  RWRN_result=specc(RWRN_similarity_fusion_list[[gama_index]],centers=CLUSTER_NUM)
  
  
  ##生存分析
  survival_data[,3]=RWR_result
  y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
  sdf=survdiff(y ~ as.character(survival_data[,3]))
  gama_result_matrix[gama_index,1]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  gama_result_matrix[gama_index,2]=dunn(-log(RWR_similarity_fusion_temp_list[[gama_index]]),RWR_result)
  
  survival_data[,3]=RWRN_result
  y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
  sdf=survdiff(y ~ as.character(survival_data[,3]))
  gama_result_matrix[gama_index,3]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  gama_result_matrix[gama_index,4]=dunn(-log(RWRN_similarity_fusion_temp_list[[gama_index]]),RWRN_result)
  
}

write.table(gama_result_matrix,file = paste("gama_result_matrix",subpath,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)

####neighbor_num####
neighbor_result_matrix=matrix(data = 0,nrow = 6,ncol = 2)
colnames(neighbor_result_matrix)=c("RWRN-p-value","RWRN-Dunn")
rownames(neighbor_result_matrix)=neighbor_num_list
RWRN_similarity_fusion_list=NULL
RWRN_similarity_fusion_temp_list=NULL
for (neighbor_num_index in 1:length(neighbor_num_list))
{
  ##融合
  time1=Sys.time()
  RWRN_similarity_fusion=RWR_fusion_neighbor(sim_list = affinityL,neighbor_num=neighbor_num_list[neighbor_num_index])
  time2=Sys.time()
  print(time2-time1)
  RWRN_similarity_fusion_temp=RWRN_similarity_fusion
  diag(RWRN_similarity_fusion)=0
  RWRN_similarity_fusion_list=c(RWRN_similarity_fusion_list,list(RWRN_similarity_fusion))
  RWRN_similarity_fusion_temp_list=c(RWRN_similarity_fusion_temp_list,list(RWRN_similarity_fusion_temp))
}
for (neighbor_num_index in 1:length(neighbor_num_list))
{
  ##聚类
  RWRN_result=specc(RWRN_similarity_fusion_list[[neighbor_num_index]],centers=CLUSTER_NUM)


  ##生存分析
  survival_data[,3]=RWRN_result
  y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
  sdf=survdiff(y ~ as.character(survival_data[,3]))
  neighbor_result_matrix[neighbor_num_index,1]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  neighbor_result_matrix[neighbor_num_index,2]=dunn(-log(RWRN_similarity_fusion_temp_list[[neighbor_num_index]]),RWRN_result)

}

write.table(neighbor_result_matrix,file = paste("neighbor_result_matrix",subpath,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)

####alpha####
alpha_result_matrix=matrix(data = 0,nrow = 9,ncol = 2)
colnames(alpha_result_matrix)=c("RWRN-p-value","RWRN-Dunn")
rownames(alpha_result_matrix)=alpha_list
RWRN_similarity_fusion_list=NULL
RWRN_similarity_fusion_temp_list=NULL
for (alpha_num_index in 1:length(alpha_list))
{
  ##融合
  time1=Sys.time()
  RWRN_similarity_fusion=RWR_fusion_neighbor(sim_list = affinityL,alpha=alpha_list[alpha_num_index])
  time2=Sys.time()
  print(time2-time1)
  RWRN_similarity_fusion_temp=RWRN_similarity_fusion
  diag(RWRN_similarity_fusion)=0
  RWRN_similarity_fusion_list=c(RWRN_similarity_fusion_list,list(RWRN_similarity_fusion))
  RWRN_similarity_fusion_temp_list=c(RWRN_similarity_fusion_temp_list,list(RWRN_similarity_fusion_temp))
}
for (alpha_num_index in 1:length(alpha_list))
{
  ##聚类
  RWRN_result=specc(RWRN_similarity_fusion_list[[alpha_num_index]],centers=CLUSTER_NUM)


  ##生存分析
  survival_data[,3]=RWRN_result
  y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
  sdf=survdiff(y ~ as.character(survival_data[,3]))
  alpha_result_matrix[alpha_num_index,1]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  alpha_result_matrix[alpha_num_index,2]=dunn(-log(RWRN_similarity_fusion_temp_list[[alpha_num_index]]),RWRN_result)

}

write.table(alpha_result_matrix,file = paste("alpha_result_matrix",subpath,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)

####beta####
beta_result_matrix=matrix(data = 0,nrow = 9,ncol = 2)
colnames(beta_result_matrix)=c("RWRN-p-value","RWRN-Dunn")
rownames(beta_result_matrix)=beta_list
RWRN_similarity_fusion_list=NULL
RWRN_similarity_fusion_temp_list=NULL
for (beta_num_index in 1:length(beta_list))
{
  ##融合
  time1=Sys.time()
  RWRN_similarity_fusion=RWR_fusion_neighbor(sim_list = affinityL,beta=beta_list[beta_num_index])
  time2=Sys.time()
  print(time2-time1)
  RWRN_similarity_fusion_temp=RWRN_similarity_fusion
  diag(RWRN_similarity_fusion)=0
  RWRN_similarity_fusion_list=c(RWRN_similarity_fusion_list,list(RWRN_similarity_fusion))
  RWRN_similarity_fusion_temp_list=c(RWRN_similarity_fusion_temp_list,list(RWRN_similarity_fusion_temp))
}
for (beta_num_index in 1:length(beta_list))
{
  ##聚类
  RWRN_result=specc(RWRN_similarity_fusion_list[[beta_num_index]],centers=CLUSTER_NUM)


  ##生存分析
  survival_data[,3]=RWRN_result
  y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
  sdf=survdiff(y ~ as.character(survival_data[,3]))
  beta_result_matrix[beta_num_index,1]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  beta_result_matrix[beta_num_index,2]=dunn(-log(RWRN_similarity_fusion_temp_list[[beta_num_index]]),RWRN_result)

}

write.table(beta_result_matrix,file = paste("beta_result_matrix",subpath,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)
stopCluster(cl)
