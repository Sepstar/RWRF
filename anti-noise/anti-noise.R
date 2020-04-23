source("home/RWR_fusion_neighbor_new.R",encoding = "UTF-8")
source("home/RWR_fusion_new.R",encoding = "UTF-8")
library(kernlab)
library(SNFtool)
library(pheatmap)
library(foreach)
library(doParallel)
cl<-makeCluster(8)  
registerDoParallel(cl) 
set.seed(666)
####Build raw virtual data####
v1=c(runif(100,min = 1,max = 5),runif(100,min = 5,max = 9))
v2=c(runif(100,min = 1,max = 5),runif(100,min = 5,max = 9))
data_original=cbind(v1,v2)
rownames(data_original)=c(paste("sample",seq(from=1,to=200),sep = ""))
truelabel = c(matrix(1,100,1),matrix(2,100,1))
plot(data_original,col=truelabel,main = 'data type')
original_dist=dist(data_original)
original_dist=(original_dist-min(original_dist))/(max(original_dist)-min(original_dist))
original_affinity=1-original_dist

sd_list=seq(from=0.2,to=3,by=0.2)
sd_NMI_matrix=matrix(data = 0,nrow = length(sd_list),ncol = 4)
colnames(sd_NMI_matrix)=c("RWR","RWRN","SNF","Concatenation")
rownames(sd_NMI_matrix)=sd_list
for (sd_index in 1:length(sd_list))
{
  print(sd_index)
  ####Add Gaussian noise####
  data_gaussian=data_original+rnorm(400, mean = 0, sd = sd_list[sd_index])
  plot(data_gaussian,col=truelabel,main = 'data type gaussian')
  
  ####Add Gama noise####
  gamma_noise=rgamma(400,shape = 1,scale = 1/2)
  gamma_noise[which(gamma_noise<0)]=0
  data_gamma=data_original+gamma_noise
  plot(data_gamma,col=truelabel,main = 'data type gamma')
  
  dataL=c(list(data_gaussian),list(data_gamma))
  
  ##Similarity matrix based on simple Euclidean distance
  affinity_Eu=NULL
  dist_Eu=NULL
  for(i in 1:length(dataL))
  {
    dist_temp=as.matrix(dist(dataL[[i]]))
    dist_temp=(dist_temp-min(dist_temp))/(max(dist_temp)-min(dist_temp))
    affinity_temp=1-dist_temp
    dist_Eu=c(dist_Eu,list(dist_temp))
    affinity_Eu=c(affinity_Eu,list(affinity_temp))
    rm(dist_temp)
    rm(affinity_temp)
  }
  
  ##Using SNFtool package to construct similarity matrix
  distL = lapply(dataL, function(x) dist2(x, x))
  affinityL = lapply(distL, function(x) affinityMatrix(x, 20, 0.5))
  
  ##Constructing similarity matrix for Concatenation using simple Euclidean distance
  dataL_Concatenation_Eu=NULL
  for(i in 1:length(dataL))
  {
    dataL_Concatenation_Eu=cbind(dataL_Concatenation_Eu,dataL[[i]])
  }
  distL_Concatenation_Eu=as.matrix(dist(dataL_Concatenation_Eu))
  distL_Concatenation_Eu=(distL_Concatenation_Eu-min(distL_Concatenation_Eu))/(max(distL_Concatenation_Eu)-min(distL_Concatenation_Eu))
  affinityL_Concatenation_Eu=1-distL_Concatenation_Eu
  
  ####Fusion of data with two types of noise####
  RWR_similarity_fusion=RWR_fusion(sim_list = affinityL)
  RWRN_similarity_fusion=RWR_fusion_neighbor(sim_list = affinityL)
  SNF_similarity_fusion=SNF(affinityL, 20, 20)
  
  diag(RWR_similarity_fusion)=0
  diag(RWRN_similarity_fusion)=0
  diag(SNF_similarity_fusion)=0
  
  RWR_result=specc(RWR_similarity_fusion,centers=2)
  RWRN_result=specc(RWRN_similarity_fusion,centers=2)
  SNF_result=spectralClustering(SNF_similarity_fusion,2)
  Concatenation_result=spectralClustering(affinityL_Concatenation_Eu,2)
  
  
  sd_NMI_matrix[sd_index,1]=calNMI(RWR_result@.Data,truelabel)
  sd_NMI_matrix[sd_index,2]=calNMI(RWRN_result@.Data,truelabel)
  sd_NMI_matrix[sd_index,3]=calNMI(SNF_result,truelabel)
  sd_NMI_matrix[sd_index,4]=calNMI(Concatenation_result,truelabel)
}


beta_list=seq(from=0.2,to=3,by=0.2)
beta_NMI_matrix=matrix(data = 0,nrow = length(beta_list),ncol = 4)
colnames(beta_NMI_matrix)=c("RWR","RWRN","SNF","Concatenation")
rownames(beta_NMI_matrix)=beta_list
for (beta_index in 1:length(beta_list))
{
  print(beta_index)
  ####Add Gaussian noise####
  data_gaussian=data_original+rnorm(400, mean = 0, sd = 2)
  plot(data_gaussian,col=truelabel,main = 'data type gaussian')
  
  ####Add Gama noise####
  gamma_noise=rgamma(400,shape = 1,scale = beta_list[beta_index])
  gamma_noise[which(gamma_noise<0)]=0
  data_gamma=data_original+gamma_noise
  plot(data_gamma,col=truelabel,main = 'data type gamma')
  
  ##Similarity matrix based on simple Euclidean distance
  dataL=c(list(data_gaussian),list(data_gamma))
  affinity_Eu=NULL
  dist_Eu=NULL
  for(i in 1:length(dataL))
  {
    dist_temp=as.matrix(dist(dataL[[i]]))
    dist_temp=(dist_temp-min(dist_temp))/(max(dist_temp)-min(dist_temp))
    affinity_temp=1-dist_temp
    dist_Eu=c(dist_Eu,list(dist_temp))
    affinity_Eu=c(affinity_Eu,list(affinity_temp))
    rm(dist_temp)
    rm(affinity_temp)
  }
  
  ##Using SNFtool package to construct similarity matrix
  distL = lapply(dataL, function(x) dist2(x, x))
  affinityL = lapply(distL, function(x) affinityMatrix(x, 20, 0.5))
  
  ##Constructing similarity matrix for Concatenation using simple Euclidean distance
  dataL_Concatenation_Eu=NULL
  for(i in 1:length(dataL))
  {
    dataL_Concatenation_Eu=cbind(dataL_Concatenation_Eu,dataL[[i]])
  }
  distL_Concatenation_Eu=as.matrix(dist(dataL_Concatenation_Eu))
  distL_Concatenation_Eu=(distL_Concatenation_Eu-min(distL_Concatenation_Eu))/(max(distL_Concatenation_Eu)-min(distL_Concatenation_Eu))
  affinityL_Concatenation_Eu=1-distL_Concatenation_Eu
  
  ####Fusion of data with two types of noise####
  RWR_similarity_fusion=RWR_fusion(sim_list = affinityL)
  RWRN_similarity_fusion=RWR_fusion_neighbor(sim_list = affinityL)
  SNF_similarity_fusion=SNF(affinityL, 20, 20)
  
  diag(RWR_similarity_fusion)=0
  diag(RWRN_similarity_fusion)=0
  diag(SNF_similarity_fusion)=0
  
  RWR_result=specc(RWR_similarity_fusion,centers=2)
  RWRN_result=specc(RWRN_similarity_fusion,centers=2)
  SNF_result=spectralClustering(SNF_similarity_fusion,2)
  Concatenation_result=spectralClustering(affinityL_Concatenation_Eu,2)
  
  
  beta_NMI_matrix[beta_index,1]=calNMI(RWR_result,truelabel)
  beta_NMI_matrix[beta_index,2]=calNMI(RWRN_result,truelabel)
  beta_NMI_matrix[beta_index,3]=calNMI(SNF_result,truelabel)
  beta_NMI_matrix[beta_index,4]=calNMI(Concatenation_result,truelabel)
}
save.image(file = "anti-noise.RData")

#load("home/anti-noise.RData")
colors_border_method=c( rgb(178/255,34/255,34/255,0.8), rgb(240/255,128/255,128/255,0.8), rgb(135/255,206/255,250/255,0.8), rgb(100/255,149/255,237/255,0.8), rgb(0.5,0.5,0.5,0.8) )
plot(sd_NMI_matrix[,1]~rownames(sd_NMI_matrix) , type="b" , bty="l" , xlab=expression(sigma), ylab="NMI" , col=colors_border_method[1] , lwd=3 , pch=16 , ylim=c(0,1) )
lines(sd_NMI_matrix[,2]~rownames(sd_NMI_matrix) , col=colors_border_method[2] , lwd=3 , pch=17 , type="b" )
# lines(sd_NMI_matrix[,3]~rownames(sd_NMI_matrix) , col=colors_border_method[3] , lwd=3 , pch=18 , type="b" )
# lines(sd_NMI_matrix[,4]~rownames(sd_NMI_matrix) , col=colors_border_method[4] , lwd=3 , pch=18 , type="b" )
legend("bottomleft",legend = c("RWRF", "RWRNF"),col = colors_border_method,
       pch = c(16,17,18,19),bty = "n",pt.cex = 1,cex = 1.2,text.col = "black",horiz = F ,inset = c(0.7, 0.1))
# legend("bottomleft",legend = c("RWRF", "RWRNF","SNF","Concatenation"),col = colors_border_method,
#        pch = c(16,17,18,19),bty = "n",pt.cex = 1,cex = 1.2,text.col = "black",horiz = F ,inset = c(0.7, 0.1))


plot(beta_NMI_matrix[,1]~rownames(beta_NMI_matrix) , type="b" , bty="l" , xlab=expression(beta), ylab="NMI" , col=colors_border_method[1] , lwd=3 , pch=16 , ylim=c(0,1) )
lines(beta_NMI_matrix[,2]~rownames(beta_NMI_matrix) , col=colors_border_method[2] , lwd=3 , pch=17 , type="b" )
# lines(beta_NMI_matrix[,3]~rownames(beta_NMI_matrix) , col=colors_border_method[3] , lwd=3 , pch=18 , type="b" )
# lines(beta_NMI_matrix[,4]~rownames(beta_NMI_matrix) , col=colors_border_method[4] , lwd=3 , pch=18 , type="b" )
legend("bottomleft",legend = c("RWRF", "RWRNF"),col = colors_border_method,
       pch = c(16,17,18,19),bty = "n",pt.cex = 1,cex = 1.2,text.col = "black",horiz = F ,inset = c(0.7, 0.1))
# legend("bottomleft",legend = c("RWRF", "RWRNF","SNF","Concatenation"),col = colors_border_method,
#        pch = c(16,17,18,19),bty = "n",pt.cex = 1,cex = 1.2,text.col = "black",horiz = F ,inset = c(0.7, 0.1))


data_original=cbind(v1,v2)
rownames(data_original)=c(paste("sample",seq(from=1,to=200),sep = ""))
truelabel = c(matrix(1,100,1),matrix(2,100,1))
plot(data_original,col=truelabel,main = 'data type')

data_gaussian=data_original+rnorm(400, mean = 0, sd = 1)
plot(data_gaussian,col=truelabel,main = 'data type gaussian')

####Add Gama noise####
gamma_noise=rgamma(400,shape = 1,scale = 1/2)
gamma_noise[which(gamma_noise<0)]=0
data_gamma=data_original+gamma_noise
plot(data_gamma,col=truelabel,main = 'data type gamma')

dataL=c(list(data_gaussian),list(data_gamma))

distL = lapply(dataL, function(x) dist2(x, x))
affinityL = lapply(distL, function(x) affinityMatrix(x, 20, 0.5))
RWR_similarity_fusion=RWR_fusion(sim_list = affinityL)
RWRN_similarity_fusion=RWR_fusion_neighbor(sim_list = affinityL)

diag(RWR_similarity_fusion)=0
diag(RWRN_similarity_fusion)=0

pheatmap(original_affinity,cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
pheatmap(affinity_Eu[[1]],cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
pheatmap(affinity_Eu[[2]],cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
pheatmap(RWR_similarity_fusion,cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
pheatmap(RWRN_similarity_fusion,cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
pheatmap(SNF_similarity_fusion,cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
pheatmap(affinityL_Concatenation_Eu,cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))

stopCluster(cl)

