source("home/RWR_fusion_neighbor_new.R",encoding = "UTF-8")
source("home/RWR_fusion_new.R",encoding = "UTF-8")
library(kernlab)
library(SNFtool)
library(pheatmap)
library(survival)

if (0)
{
  ####Build raw virtual data####
  v1=c(runif(100,min = 1,max = 5),runif(100,min = 5,max = 9))
  v2=c(runif(100,min = 1,max = 5),runif(100,min = 5,max = 9))
  data_original=cbind(v1,v2)
  rownames(data_original)=c(paste("sample",seq(from=1,to=200),sep = ""))
  truelabel = c(matrix(1,100,1),matrix(2,100,1))
  
  
  data_gaussian=data_original+rnorm(400, mean = 0, sd = 2)
  gamma_noise=rgamma(400,shape = 1,scale = 2)
  gamma_noise[which(gamma_noise<0)]=0
  data_gamma=data_original+gamma_noise
  
  dataL=c(list(data_gaussian),list(data_gamma))
  distL = lapply(dataL, function(x) dist2(x, x))
  affinityL = lapply(distL, function(x) affinityMatrix(x, 20, 0.5))
  
  
  
  NMI_matrix=matrix(data = 0,nrow = 30,ncol = 2)
  colnames(NMI_matrix)=c("RWR","RWRN")
  rownames(NMI_matrix)=1:30
  ####Test RWR convergence####
  for (i in 1:5)
  {
    print(i)
    RWR_similarity_fusion=RWR_fusion(sim_list = affinityL,iteration_max = i)[[1]]
    RWRN_similarity_fusion=RWR_fusion_neighbor(sim_list = affinityL,iteration_max = i)[[1]]
    diag(RWR_similarity_fusion)=0
    diag(RWRN_similarity_fusion)=0
    RWR_result=specc(RWR_similarity_fusion,centers=2)
    RWRN_result=specc(RWRN_similarity_fusion,centers=2)
    NMI_matrix[i,1]=calNMI(RWR_result,truelabel)
    NMI_matrix[i,2]=calNMI(RWRN_result,truelabel)
  }
}
####ACC####
load("home/result_2ACC.RData")
rm(RWR_fusion)
rm(RWR_fusion_neighbor)
library(foreach)
library(doParallel)
cl<-makeCluster(8)  
registerDoParallel(cl) 

library(clValid)
p_matrix=matrix(data = 0,nrow = 20,ncol = 2)
colnames(p_matrix)=c("RWR","RWRN")
rownames(p_matrix)=1:20
dunn_matrix=p_matrix
relative_error_matrix=p_matrix
RWR_similarity_fusion_list=NULL
RWRN_similarity_fusion_list=NULL
RWR_similarity_fusion_temp_list=NULL
RWRN_similarity_fusion_temp_list=NULL
####Test RWR convergence####
for (i in 1:20)
{
  print(i)
  RWR_similarity_fusion=RWR_fusion(sim_list = affinityL,iteration_max = i)
  RWRN_similarity_fusion=RWR_fusion_neighbor(sim_list = affinityL,iteration_max = i)
  relative_error_matrix[i,1]=RWR_similarity_fusion[[2]]
  relative_error_matrix[i,2]=RWRN_similarity_fusion[[2]]
  RWR_similarity_fusion=RWR_similarity_fusion[[1]]
  RWRN_similarity_fusion=RWRN_similarity_fusion[[1]]
  RWR_similarity_fusion_temp=RWR_similarity_fusion
  RWRN_similarity_fusion_temp=RWRN_similarity_fusion
  diag(RWR_similarity_fusion)=0
  diag(RWRN_similarity_fusion)=0
  RWR_similarity_fusion_list=c(RWR_similarity_fusion_list,list(RWR_similarity_fusion))
  RWRN_similarity_fusion_list=c(RWRN_similarity_fusion_list,list(RWRN_similarity_fusion))
  RWR_similarity_fusion_temp_list=c(RWR_similarity_fusion_temp_list,list(RWR_similarity_fusion_temp))
  RWRN_similarity_fusion_temp_list=c(RWRN_similarity_fusion_temp_list,list(RWRN_similarity_fusion_temp))
}

if (0)
{
  for (i in 1:20)
  {
    print(i)  
    RWR_result=specc(RWR_similarity_fusion_list[[i]],centers=2)
    RWRN_result=specc(RWRN_similarity_fusion_list[[i]],centers=2)
    
    survival_data[,3]=RWR_result
    y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
    sdf=survdiff(y ~ as.character(survival_data[,3]))
    p_matrix[i,1]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
    dunn_matrix[i,1]=dunn(-log(RWR_similarity_fusion_temp_list[[i]]),RWR_result)
    
    survival_data[,3]=RWRN_result
    y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
    sdf=survdiff(y ~ as.character(survival_data[,3]))
    p_matrix[i,2]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
    dunn_matrix[i,2]=dunn(-log(RWRN_similarity_fusion_temp_list[[i]]),RWRN_result)
  }
}
log10_relative_error_matrix=-log10(relative_error_matrix)
plot( log10_relative_error_matrix[,1]~rownames(log10_relative_error_matrix) , type="b" , bty="l" , xlab="Number of Iterations" , ylab="-lg(Relative Error)" , col=rgb(178/255,34/255,34/255,0.8) , lwd=3 , pch=16 , ylim=c(0,12) )
lines(log10_relative_error_matrix[,2]~rownames(log10_relative_error_matrix) , col=rgb(240/255,128/255,128/255,0.8) , lwd=3 , pch=16 , type="b" )
legend("bottomleft",legend = c("RWRF", "RWRNF"),col = c(rgb(178/255,34/255,34/255,0.8), rgb(240/255,128/255,128/255,0.8)),
       pch = c(16,16),bty = "n",pt.cex = 1,cex = 1.2,text.col = "black",horiz = F ,inset = c(0.8, 0.1))

####4BLCA####
load("E:/研究生科研工作/2017异质网络随机游走做融合/4_跑一遍程序/4_重新海选20180710/results/筛选结果/4BLCA/result_4BLCA.RData")

library(clValid)
p_matrix=matrix(data = 0,nrow = 20,ncol = 2)
colnames(p_matrix)=c("RWR","RWRN")
rownames(p_matrix)=1:20
dunn_matrix=p_matrix
relative_error_matrix=p_matrix
RWR_similarity_fusion_list=NULL
RWRN_similarity_fusion_list=NULL
RWR_similarity_fusion_temp_list=NULL
RWRN_similarity_fusion_temp_list=NULL
####测试RWR的收敛性####
for (i in 1:20)
{
  print(i)
  RWR_similarity_fusion=RWR_fusion(sim_list = affinityL,iteration_max = i)
  RWRN_similarity_fusion=RWR_fusion_neighbor(sim_list = affinityL,iteration_max = i)
  relative_error_matrix[i,1]=RWR_similarity_fusion[[2]]
  relative_error_matrix[i,2]=RWRN_similarity_fusion[[2]]
  RWR_similarity_fusion=RWR_similarity_fusion[[1]]
  RWRN_similarity_fusion=RWRN_similarity_fusion[[1]]
  RWR_similarity_fusion_temp=RWR_similarity_fusion
  RWRN_similarity_fusion_temp=RWRN_similarity_fusion
  diag(RWR_similarity_fusion)=0
  diag(RWRN_similarity_fusion)=0
  RWR_similarity_fusion_list=c(RWR_similarity_fusion_list,list(RWR_similarity_fusion))
  RWRN_similarity_fusion_list=c(RWRN_similarity_fusion_list,list(RWRN_similarity_fusion))
  RWR_similarity_fusion_temp_list=c(RWR_similarity_fusion_temp_list,list(RWR_similarity_fusion_temp))
  RWRN_similarity_fusion_temp_list=c(RWRN_similarity_fusion_temp_list,list(RWRN_similarity_fusion_temp))
}

if (0)
{
  for (i in 1:20)
  {
    print(i)  
    RWR_result=specc(RWR_similarity_fusion_list[[i]],centers=2)
    RWRN_result=specc(RWRN_similarity_fusion_list[[i]],centers=2)
    
    survival_data[,3]=RWR_result
    y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
    sdf=survdiff(y ~ as.character(survival_data[,3]))
    p_matrix[i,1]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
    dunn_matrix[i,1]=dunn(-log(RWR_similarity_fusion_temp_list[[i]]),RWR_result)
    
    survival_data[,3]=RWRN_result
    y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
    sdf=survdiff(y ~ as.character(survival_data[,3]))
    p_matrix[i,2]=1-pchisq(sdf$chisq, length(sdf$n) - 1)
    dunn_matrix[i,2]=dunn(-log(RWRN_similarity_fusion_temp_list[[i]]),RWRN_result)
  }
}
log10_relative_error_matrix=-log10(relative_error_matrix)
plot( log10_relative_error_matrix[,1]~rownames(log10_relative_error_matrix) , type="b" , bty="l" , xlab="Number of Iterations" , ylab="-lg(Relative Error)" , col=rgb(178/255,34/255,34/255,0.8) , lwd=3 , pch=16 , ylim=c(0,12) )
lines(log10_relative_error_matrix[,2]~rownames(log10_relative_error_matrix) , col=rgb(240/255,128/255,128/255,0.8) , lwd=3 , pch=16 , type="b" )
legend("bottomleft",legend = c("RWRF", "RWRNF"),col = c(rgb(178/255,34/255,34/255,0.8), rgb(240/255,128/255,128/255,0.8)),
       pch = c(16,16),bty = "n",pt.cex = 1,cex = 1.2,text.col = "black",horiz = F ,inset = c(0.8, 0.1))

