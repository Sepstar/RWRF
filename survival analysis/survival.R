library(FactoMineR)
library(rgl)
library(survival)

####2ACC####
load("home/result_2ACC.RData")

survival_data[,3]=cluster_list[[2]][[1]]
y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
kaplan_meier <- survfit(y ~ as.character(survival_data[,3]))
colour <- c('steelblue','yellow4', 'firebrick3', 'skyblue', "salmon", "darkseagreen3", "springgreen4","tan4","whitesmoke","orchid2"  )
plot(kaplan_meier, col=colour,xlab='Survival Time in Months',ylab='Survival Probabilities', lwd=2,xlim=c(0,160))
leg.txt<-c("Subtype 1","Subtype 2","Subtype 3")
legend("topright",leg.txt,col=colour,lwd=3,lty=1)
sdf=survdiff(y ~ as.character(survival_data[,3]))
1-pchisq(sdf$chisq, length(sdf$n) - 1)

survival_data[,3]=cluster_list[[2]][[2]]
y=Surv(time = as.numeric(survival_data[,1]), event = as.numeric(survival_data[,2]))
kaplan_meier <- survfit(y ~ as.character(survival_data[,3]))
colour <- c('steelblue','yellow4', 'firebrick3', 'skyblue', "salmon", "darkseagreen3", "springgreen4","tan4","whitesmoke","orchid2"  )
plot(kaplan_meier, col=colour,xlab='Survival Time in Months',ylab='Survival Probabilities', lwd=2,xlim=c(0,160))
leg.txt<-c("Subtype 1","Subtype 2","Subtype 3")
legend("topright",leg.txt,col=colour,lwd=3,lty=1)
sdf=survdiff(y ~ as.character(survival_data[,3]))
1-pchisq(sdf$chisq, length(sdf$n) - 1)
#rm(list = ls())

data_PCA=cbind(as.data.frame(similarity_fusion_list[[1]]),cluster_list[[2]][[1]])
res.PCA = PCA(data_PCA, scale.unit=TRUE, ncp=4, graph=F ,quali.sup=nrow(similarity_fusion_list[[1]])+1) 
plot(res.PCA$ind$coord[,1] , res.PCA$ind$coord[,2]  , xlab="PC1" , ylab="PC2" , pch=20 , cex=1 , 
     col=colour[as.numeric(res.PCA$call$quali.sup$quali.sup[,1])] )
open3d()
plot3d(res.PCA$ind$coord[,1] , res.PCA$ind$coord[,2],res.PCA$ind$coord[,3],size=5 , 
       xlab="PC1" , ylab="PC2",zlab = "PC3" ,col=colour[as.numeric(res.PCA$call$quali.sup$quali.sup[,1])])

data_PCA=cbind(as.data.frame(similarity_fusion_list[[2]]),cluster_list[[2]][[2]])
res.PCA = PCA(data_PCA, scale.unit=TRUE, ncp=4, graph=F ,quali.sup=nrow(similarity_fusion_list[[2]])+1) 
plot(res.PCA$ind$coord[,1] , res.PCA$ind$coord[,2]  , xlab="PC1" , ylab="PC2" , pch=20 , cex=1 , 
     col=colour[as.numeric(res.PCA$call$quali.sup$quali.sup[,1])] )
open3d()
plot3d(res.PCA$ind$coord[,1] , res.PCA$ind$coord[,2],res.PCA$ind$coord[,3],size=5 , 
       xlab="PC1" , ylab="PC2",zlab = "PC3" ,col=colour[as.numeric(res.PCA$call$quali.sup$quali.sup[,1])])

rgl.snapshot("plot1.png")
rgl.postscript("plot1.pdf", "pdf", drawText = TRUE)

rm(list = ls())

# relapse-free survival and survival after relapse
load("home/result_2ACC.RData")
# RWR
survival_data[,3]=cluster_list[[2]][[1]]
ACC_clinical=t(read.table(file="home/ACC_clinicalMatrix",header = T,sep = "\t",row.names = 1,check.names = F))
new_survival_data=NULL
for (i in 1:nrow(survival_data))
{
  num=which(rownames(survival_data)[i]==colnames(ACC_clinical))
  temp1=ACC_clinical[7,num]
  temp2=ACC_clinical[8,num]
  tempall = cbind(survival_data[i,],t(temp1),t(temp2))
  new_survival_data=rbind(new_survival_data,tempall)
  rm(temp1,temp2,tempall)
  rm(num)
}
new_survival_data=as.data.frame(new_survival_data)
colnames(new_survival_data)[3:5]=c("Subtype","_RFS","_RFS_IND")

# relapse-free survival
new_survival_data2 = data.frame(sapply(new_survival_data,trimws,which="both"), stringsAsFactors=FALSE)
new_survival_data3 = new_survival_data2[-which(new_survival_data2[,5] == "1"),]

y=Surv(time = as.numeric(new_survival_data3[,1]), event = as.numeric(new_survival_data3[,2]))
kaplan_meier <- survfit(y ~ as.character(new_survival_data3[,3]))
colour <- c('steelblue','yellow4', 'firebrick3', 'skyblue', "salmon", "darkseagreen3", "springgreen4","tan4","whitesmoke","orchid2"  )
plot(kaplan_meier, col=colour,xlab='Survival Time in Months',ylab='Survival Probabilities', lwd=2,xlim=c(0,160))
leg.txt<-c("Subtype 1","Subtype 2","Subtype 3")
legend("topright",leg.txt,col=colour,lwd=3,lty=1)
sdf=survdiff(y ~ as.character(new_survival_data3[,3]))
1-pchisq(sdf$chisq, length(sdf$n) - 1)

# survival after relapse
new_survival_data2 = data.frame(sapply(new_survival_data,trimws,which="both"), stringsAsFactors=FALSE)
new_survival_data3 = new_survival_data2[which(new_survival_data2[,5] == "1"),]
new_survival_data3[,4]=as.numeric(as.character(new_survival_data3[,4]))/30
new_survival_data3[,1]=as.numeric(as.character(new_survival_data3[,1]))
new_survival_data3=cbind(new_survival_data3,(new_survival_data3[,1]-new_survival_data3[,4]))
colnames(new_survival_data3)[6]="survival after relapse"

y=Surv(time = as.numeric(new_survival_data3[,6]), event = as.numeric(new_survival_data3[,2]))
kaplan_meier <- survfit(y ~ as.character(new_survival_data3[,3]))
colour <- c('steelblue','yellow4', 'firebrick3', 'skyblue', "salmon", "darkseagreen3", "springgreen4","tan4","whitesmoke","orchid2"  )
plot(kaplan_meier, col=colour,xlab='Survival Time in Months',ylab='Survival Probabilities', lwd=2,xlim=c(0,160))
leg.txt<-c("Subtype 1","Subtype 2","Subtype 3")
legend("topright",leg.txt,col=colour,lwd=3,lty=1)
sdf=survdiff(y ~ as.character(new_survival_data3[,3]))
1-pchisq(sdf$chisq, length(sdf$n) - 1)

