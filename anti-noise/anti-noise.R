source("home/RWR_fusion_neighbor_new.R",encoding = "UTF-8")
source("home/RWR_fusion_new.R",encoding = "UTF-8")
library(kernlab)
library(SNFtool)
library(pheatmap)

####Build raw virtual data####
v1=c(runif(100,min = 1,max = 5),runif(100,min = 5,max = 9))
v2=c(runif(100,min = 1,max = 5),runif(100,min = 5,max = 9))

# plot
data_original=cbind(v1,v2)
rownames(data_original)=c(paste("sample",seq(from=1,to=200),sep = ""))
truelabel = c(matrix(1,100,1),matrix(2,100,1))
plot(data_original,col=truelabel,main = 'data type')

data_gaussian=data_original+rnorm(400, mean = 0, sd = 1.5)
plot(data_gaussian,col=truelabel,main = 'data type gaussian')

gamma_noise=rgamma(400,shape = 3,scale = 1)
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

pheatmap(original_affinity,cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
pheatmap(affinity_Eu[[1]],cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
pheatmap(affinity_Eu[[2]],cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
pheatmap(RWR_similarity_fusion,cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
pheatmap(RWRN_similarity_fusion,cluster_rows=F,cluster_cols = F,show_rownames=F,show_colnames=F,color = colorRampPalette(c("skyblue", "firebrick1", "firebrick4"))(50))
                   
a = affinity_Eu[[1]]
b = affinity_Eu[[2]]
diag(a) = 0
diag(b) = 0

data_gaussian_result = specc(a,centers=2)
data_gamma_result = specc(b,centers=2)
RWR_result=specc(RWR_similarity_fusion,centers=2)
RWRN_result=specc(RWRN_similarity_fusion,centers=2)

calNMI(RWR_result@.Data,truelabel)#0.7563153
calNMI(RWRN_result@.Data,truelabel)#0.7563153
calNMI(data_gaussian_result@.Data,truelabel)#0.6009723
calNMI(data_gamma_result@.Data,truelabel)#0.5550035

