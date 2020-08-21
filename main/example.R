source("/home/wenyuqi/RWR/RWR_fusion_neighbor_new.R")
source("/home/wenyuqi/RWR/RWR_fusion_new.R")
library(SNFtool)
library(survival)
library(clValid)
library(kernlab)
library(foreach)
library(doParallel)
cl<-makeCluster(8)  
registerDoParallel(cl) 

####Setting parameters####
subpath="2ACC"
CLUSTER_NUM_LIST=2:6
####2ACC####
mrna_data=as.matrix(read.table(file = paste(subpath,"/HiSeqV2_PANCAN2.txt",sep = ""),header = T,sep = "\t",quote = "",check.names = F))
print("read mrna")
mirna_data=as.matrix(read.table(file = paste(subpath,"/miRNA_HiSeq_gene1.txt",sep = ""),header = T,sep = "\t",quote = "",check.names = F))
print("read mirna")
meth_data=as.matrix(read.table(file = paste(subpath,"/HumanMethylation4501.txt",sep = ""),header = T,sep = "\t",quote = "",check.names = F))
print("read meth")
#cnv_data=as.matrix(read.table(file = paste(subpath,"/Gistic2_CopyNumber_Gistic2_all_data_by_genes1.txt",sep = ""),header = T,sep = "\t",quote = "",check.names = F))
mrna_data=t(mrna_data)
mirna_data=t(mirna_data)
meth_data=t(meth_data)
#cnv_data=t(cnv_data)
dataL=list(mrna_data)
dataL=c(dataL,list(mirna_data))
dataL=c(dataL,list(meth_data))
#dataL=c(dataL,list(cnv_data))

result_matrix=matrix(data = 0,nrow = 7,ncol = 10)
colnames(result_matrix)=c("p-value_2","Dunn_2","p-value_3","Dunn_3","p-value_4","Dunn_4","p-value_5","Dunn_5","p-value_6","Dunn_6")
rownames(result_matrix)=c("RWR_AFF","RWRN_AFF","SNF","Concatenation","mRNA","miRNA","Methylation")

survival_data=read.table(file = paste(subpath,"/clinicalMatrix15.txt",sep = ""),header = TRUE, sep="\t")
survival_data[,1]=survival_data[,1]/30
survival_data=cbind(survival_data,survival_data[,2])
colnames(survival_data)[3]="Cluster"
print("preprocess")

##Construct similarity matrix
distL = lapply(dataL, function(x) dist2(x, x))
affinityL = lapply(distL, function(x) affinityMatrix(x, 20, 0.5))
print("AFF")

##Constructing similarity matrix for Concatenation using simple Euclidean distance
dataL_Concatenation_Eu=NULL
for(i in 1:length(dataL))
{
  dataL_Concatenation_Eu=cbind(dataL_Concatenation_Eu,dataL[[i]])
}
distL_Concatenation_Eu=as.matrix(dist(dataL_Concatenation_Eu))
distL_Concatenation_Eu=(distL_Concatenation_Eu-min(distL_Concatenation_Eu))/(max(distL_Concatenation_Eu)-min(distL_Concatenation_Eu))
affinityL_Concatenation_Eu=1-distL_Concatenation_Eu

##Fusion
similarity_fusion_list=NULL

temp=RWR_fusion(sim_list = affinityL)
similarity_fusion_list=c(similarity_fusion_list,list(temp))
rm(temp)

temp=RWR_fusion_neighbor(sim_list = affinityL)
similarity_fusion_list=c(similarity_fusion_list,list(temp))
rm(temp)

time1=Sys.time()
temp=SNF(affinityL, 20, 20)
time2=Sys.time()
print(time2-time1)
similarity_fusion_list=c(similarity_fusion_list,list(temp))
rm(temp)

similarity_fusion_list=c(similarity_fusion_list,list(affinityL_Concatenation_Eu))

similarity_fusion_list=c(similarity_fusion_list,list(affinityL[[1]]))

similarity_fusion_list=c(similarity_fusion_list,list(affinityL[[2]]))

similarity_fusion_list=c(similarity_fusion_list,list(affinityL[[3]]))

similarity_fusion_list_temp=similarity_fusion_list
for (i in 1:length(similarity_fusion_list))
{
  diag(similarity_fusion_list[[i]])=0
}

##Clustering
set.seed(7)
print("Clustering")
cluster_list=NULL
for (CLUSTER_NUM in CLUSTER_NUM_LIST)
{
  print(CLUSTER_NUM)
  temp=NULL
  temp=c(temp,list(specc(similarity_fusion_list[[1]],centers=CLUSTER_NUM)@.Data))
  temp=c(temp,list(specc(similarity_fusion_list[[2]],centers=CLUSTER_NUM)@.Data))
  temp=c(temp,list(spectralClustering(similarity_fusion_list[[3]],CLUSTER_NUM)))
  temp=c(temp,list(spectralClustering(similarity_fusion_list[[4]],CLUSTER_NUM)))
  temp=c(temp,list(specc(similarity_fusion_list[[5]],centers=CLUSTER_NUM)@.Data))
  temp=c(temp,list(specc(similarity_fusion_list[[6]],centers=CLUSTER_NUM)@.Data))
  temp=c(temp,list(specc(similarity_fusion_list[[7]],centers=CLUSTER_NUM)@.Data))
  cluster_list=c(cluster_list,list(temp))
  rm(temp)
}

##Survival analysis and cluster validation
for (CLUSTER_NUM in CLUSTER_NUM_LIST)
{
  print(CLUSTER_NUM)
  for (i in 1:length(similarity_fusion_list))
  {
    survival_data[,3]=cluster_list[[CLUSTER_NUM-1]][[i]]
    y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
    sdf=survdiff(y ~ as.character(survival_data[,3]))
    result_matrix[i,(CLUSTER_NUM*2-3)]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
    if (i<=3)
    {
      result_matrix[i,(CLUSTER_NUM*2-2)]=dunn(-log(similarity_fusion_list_temp[[i]]),cluster_list[[CLUSTER_NUM-1]][[i]])
    }else if (i==4)
    {
      result_matrix[i,(CLUSTER_NUM*2-2)]=dunn(distL_Concatenation_Eu,cluster_list[[CLUSTER_NUM-1]][[i]])
    }else if (i==5)
    {
      result_matrix[i,(CLUSTER_NUM*2-2)]=dunn(distL[[1]],cluster_list[[CLUSTER_NUM-1]][[i]])
    }else if (i==6)
    {
      result_matrix[i,(CLUSTER_NUM*2-2)]=dunn(distL[[2]],cluster_list[[CLUSTER_NUM-1]][[i]])
    }else if (i==7)
    {
      result_matrix[i,(CLUSTER_NUM*2-2)]=dunn(distL[[3]],cluster_list[[CLUSTER_NUM-1]][[i]])
    }
    rm(y)
    rm(sdf)
    hhh=similarity_fusion_list[[i]]
    hhh1_RANK=NULL
    for (j in 1:length(unique(cluster_list[[CLUSTER_NUM-1]][[i]])))
    {
      hhh1_RANK=c(hhh1_RANK,which(cluster_list[[CLUSTER_NUM-1]][[i]]==j))
    }
    hhh=hhh[hhh1_RANK,hhh1_RANK]
    jpeg(file=paste("/home/100/wenyuqi/",subpath,"/","heatmap",i,"_","C",CLUSTER_NUM,"_",subpath,".jpeg",sep = ""))
    heatmap(hhh,col = colorRampPalette(c("navy", "grey", "firebrick3"))(50),Rowv =NA,Colv =NA)
    dev.off()
    rm(hhh)
    rm(hhh1_RANK) 
  }
}

write.table(result_matrix,file = paste("/home/",subpath,"/result_matrix","_",subpath,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)
