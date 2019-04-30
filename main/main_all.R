source("E:/研究生科研工作/2017异质网络随机游走做融合/4_跑一遍程序/1_main/RWR_fusion_neighbor.R",encoding = "UTF-8")
source("E:/研究生科研工作/2017异质网络随机游走做融合/4_跑一遍程序/1_main/RWR_fusion.R",encoding = "UTF-8")
library(SNFtool)
library(survival)
library(iCluster)
library(clValid)
library(kernlab)
library(speccalt)
library(cluster)
dataL=NULL

dataL=c(dataL,list(as.matrix(read.table(file = "E:/研究生科研工作/2017异质网络随机游走做融合(何松课题)/4_跑一遍程序/1_main/Lung(SNF test dataset)/LUNG_Gene_Expression.txt",header = T,sep = "\t",row.names = 1,quote = "",check.names = F))))
dataL=c(dataL,list(as.matrix(read.table(file = "E:/研究生科研工作/2017异质网络随机游走做融合(何松课题)/4_跑一遍程序/1_main/Lung(SNF test dataset)/LUNG_Methy_Expression.txt",header = T,sep = "\t",row.names = 1,quote = "",check.names = F))))
dataL=c(dataL,list(as.matrix(read.table(file = "E:/研究生科研工作/2017异质网络随机游走做融合(何松课题)/4_跑一遍程序/1_main/Lung(SNF test dataset)/LUNG_Mirna_Expression.txt",header = T,sep = "\t",row.names = 1,quote = "",check.names = F))))
survival_data=read.table('E:/研究生科研工作/2017异质网络随机游走做融合(何松课题)/4_跑一遍程序/1_main/Lung(SNF test dataset)/LUNG_Survival.txt',header = TRUE, sep="\t")

colnames(dataL[[2]])=colnames(dataL[[1]])
colnames(dataL[[3]])=colnames(dataL[[1]])# 病人的名字对齐了吗？
for (i in 1:3)
{
  dataL[[i]]=t(dataL[[i]])
}
survival_data[,2]=survival_data[,2]/30
survival_data=cbind(survival_data,survival_data[,3])
colnames(survival_data)[4]="Cluster"
CLUSTER_NUM=4

p_matrix=matrix(data = 0,nrow = 8,ncol = 3)
colnames(p_matrix)=c("SNFtool","Kernlab","Speccalt")
rownames(p_matrix)=c("RWR_Eu","RWR","RWR_neighbor_Eu","RWR_neighbor","SNF","Concatenation_Eu","Concatenation","iCluster")
time_matrix=p_matrix
dunn_matrix=p_matrix
sill_matrix=p_matrix

##基于简单的欧氏距离构建的相似性矩阵
affinity_Eu=NULL
for(i in 1:length(dataL))
{
  dist_temp=as.matrix(dist(dataL[[i]]))
  dist_temp=(dist_temp-min(dist_temp))/(max(dist_temp)-min(dist_temp))
  affinity_temp=1-dist_temp
  affinity_Eu=c(affinity_Eu,list(affinity_temp))
  rm(dist_temp)
  rm(affinity_temp)
}

##利用SNFtool包构建相似性矩阵
distL = lapply(dataL, function(x) dist2(x, x))
affinityL = lapply(distL, function(x) affinityMatrix(x, 20, 0.5))

##利用简单的欧氏距离构建用于Concatenation的相似性矩阵
dataL_Concatenation_Eu=NULL
for(i in 1:length(dataL))
{
  dataL_Concatenation_Eu=cbind(dataL_Concatenation_Eu,dataL[[i]])
}
dist_temp=as.matrix(dist(dataL_Concatenation_Eu))
dist_temp=(dist_temp-min(dist_temp))/(max(dist_temp)-min(dist_temp))
affinityL_Concatenation_Eu=1-dist_temp
rm(dist_temp)

##利用SNFtool包构建用于Concatenation的相似性矩阵
dataL_Concatenation=NULL
for(i in 1:length(dataL))
{
  dataL_Concatenation=cbind(dataL_Concatenation,dataL[[i]])
}
distL_Concatenation=dist2(dataL_Concatenation,dataL_Concatenation)
affinityL_Concatenation = affinityMatrix(distL_Concatenation, 20, 0.5)


similarity_fusion=NULL

time1=Sys.time()
similarity_fusion=c(similarity_fusion,list(RWR_fusion(sim_list = affinity_Eu)))
time2=Sys.time()
time_matrix[1,1]=time2-time1

time1=Sys.time()
similarity_fusion=c(similarity_fusion,list(RWR_fusion(sim_list = affinityL)))
time2=Sys.time()
time_matrix[2,1]=time2-time1

time1=Sys.time()
similarity_fusion=c(similarity_fusion,list(RWR_fusion_neighbor(sim_list = affinity_Eu)))
time2=Sys.time()
time_matrix[3,1]=time2-time1

time1=Sys.time()
similarity_fusion=c(similarity_fusion,list(RWR_fusion_neighbor(sim_list = affinityL)))
time2=Sys.time()
time_matrix[4,1]=time2-time1


time1=Sys.time()
similarity_fusion=c(similarity_fusion,list(SNF(affinityL, 20, 20)))
time2=Sys.time()
time_matrix[5,1]=time2-time1

time1=Sys.time()
similarity_fusion=c(similarity_fusion,list(affinityL_Concatenation_Eu))
time2=Sys.time()
time_matrix[6,1]=time2-time1

time1=Sys.time()
similarity_fusion=c(similarity_fusion,list(affinityL_Concatenation))
time2=Sys.time()
time_matrix[7,1]=time2-time1

similarity_fusion_temp=similarity_fusion
for (i in 1:7)
{
  diag(similarity_fusion[[i]])=0
}

for (i in 1:5)
{
  print(i)
  result_cluster1=spectralClustering(similarity_fusion[[i]],CLUSTER_NUM)
  result_cluster2=specc(similarity_fusion[[i]],centers=CLUSTER_NUM)
  result_cluster3=speccalt(similarity_fusion[[i]],k=CLUSTER_NUM)
  
  survival_data[,4]=result_cluster1
  y=Surv(time = as.numeric(survival_data[,2]), event = as.numeric(survival_data[,3]))
  sdf=survdiff(y ~ as.character(survival_data[,4]))
  p_matrix[i,1]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  dunn_matrix[i,1]=dunn(-log(similarity_fusion_temp[[i]]),result_cluster1)
  sill_matrix[i,1]=summary(silhouette(result_cluster1,-log(similarity_fusion_temp[[i]])))[[4]]
  
  survival_data[,4]=result_cluster2
  y=Surv(time = as.numeric(survival_data[,2]), event = as.numeric(survival_data[,3]))
  sdf=survdiff(y ~ as.character(survival_data[,4]))
  p_matrix[i,2]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  dunn_matrix[i,2]=dunn(-log(similarity_fusion_temp[[i]]),result_cluster2)
  sill_matrix[i,2]=summary(silhouette(result_cluster2,-log(similarity_fusion_temp[[i]])))[[4]]
  
  survival_data[,4]=result_cluster3
  y=Surv(time = as.numeric(survival_data[,2]), event = as.numeric(survival_data[,3]))
  sdf=survdiff(y ~ as.character(survival_data[,4]))
  p_matrix[i,3]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  dunn_matrix[i,3]=dunn(-log(similarity_fusion_temp[[i]]),result_cluster3)
  sill_matrix[i,3]=summary(silhouette(result_cluster3,-log(similarity_fusion_temp[[i]])))[[4]]
}

for (i in 6:7)
{
  print(i)
  result_cluster1=spectralClustering(similarity_fusion[[i]],CLUSTER_NUM)
  result_cluster2=specc(similarity_fusion[[i]],centers=CLUSTER_NUM)
  result_cluster3=speccalt(similarity_fusion[[i]],k=CLUSTER_NUM)
  
  survival_data[,4]=result_cluster1
  y=Surv(time = as.numeric(survival_data[,2]), event = as.numeric(survival_data[,3]))
  sdf=survdiff(y ~ as.character(survival_data[,4]))
  p_matrix[i,1]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  dunn_matrix[i,1]=dunn(distL_Concatenation,result_cluster1)
  sill_matrix[i,1]=summary(silhouette(result_cluster1,distL_Concatenation))[[4]]
  
  survival_data[,4]=result_cluster2
  y=Surv(time = as.numeric(survival_data[,2]), event = as.numeric(survival_data[,3]))
  sdf=survdiff(y ~ as.character(survival_data[,4]))
  p_matrix[i,2]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  dunn_matrix[i,2]=dunn(distL_Concatenation,result_cluster2)
  sill_matrix[i,2]=summary(silhouette(result_cluster2,distL_Concatenation))[[4]]
  
  survival_data[,4]=result_cluster3
  y=Surv(time = as.numeric(survival_data[,2]), event = as.numeric(survival_data[,3]))
  sdf=survdiff(y ~ as.character(survival_data[,4]))
  p_matrix[i,3]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
  dunn_matrix[i,3]=dunn(distL_Concatenation,result_cluster3)
  sill_matrix[i,3]=summary(silhouette(result_cluster3,distL_Concatenation))[[4]]
}







