####gama####
gama_result=read.table(file = "gama_result_matrix2ACC.txt",header = T,sep = "\t")
gama_result[,1]=-log10(gama_result[,1])
gama_result[,3]=-log10(gama_result[,3])
plot( gama_result[,1]~rownames(gama_result) , type="b" , bty="l" , xlab=expression(gamma) , ylab=expression(-log[10](P)) , col=rgb(178/255,34/255,34/255,0.8) , lwd=3 , pch=16 , ylim=c(0,7) )
lines(gama_result[,3]~rownames(gama_result) , col=rgb(240/255,128/255,128/255,0.8) , lwd=3 , pch=16 , type="b" )
legend("bottomleft",legend = c("RWRF", "RWRNF"),col = c(rgb(178/255,34/255,34/255,0.8),rgb(240/255,128/255,128/255,0.8)),
       pch = c(16,16),bty = "n",pt.cex = 1,cex = 1.2,text.col = "black",horiz = F ,inset = c(0.8, 0.1))


plot( gama_result[,2]~rownames(gama_result) , type="b" , bty="l" , xlab=expression(gamma) , ylab="Dunn" , col=rgb(178/255,34/255,34/255,0.8) , lwd=3 , pch=16 , ylim=c(0,1) )
lines(gama_result[,4]~rownames(gama_result) , col=rgb(240/255,128/255,128/255,0.8) , lwd=3 , pch=16 , type="b" )
legend("bottomleft",legend = c("RWRF", "RWRNF"),col = c(rgb(178/255,34/255,34/255,0.8),rgb(240/255,128/255,128/255,0.8)),
       pch = c(16,16),bty = "n",pt.cex = 1,cex = 1.2,text.col = "black",horiz = F ,inset = c(0.8, 0.1))

####neighbor####
neighbor_result=read.table(file = "neighbor_result_matrix2ACC.txt",header = T,sep = "\t")
neighbor_result=neighbor_result[1:5,]
neighbor_result[,1]=-log10(neighbor_result[,1])
plot( neighbor_result[,1]~rownames(neighbor_result) , type="b" , bty="l" , xlab="m" , ylab=expression(-log[10](P)) , col=rgb(240/255,128/255,128/255,0.8) , lwd=3 , pch=16 , ylim=c(1,9) )
plot( neighbor_result[,2]~rownames(neighbor_result) , type="b" , bty="l" , xlab="m" , ylab="Dunn" , col=rgb(240/255,128/255,128/255,0.8) , lwd=3 , pch=16 , ylim=c(0,1) )

####alpha####
alpha_result=read.table(file = "alpha_result_matrix2ACC.txt",header = T,sep = "\t")
alpha_result[,1]=-log10(alpha_result[,1])
plot( alpha_result[,1]~rownames(alpha_result) , type="b" , bty="l" , xlab=expression(alpha) , ylab=expression(-log[10](P)) , col=rgb(240/255,128/255,128/255,0.8) , lwd=3 , pch=16 , ylim=c(0,8) )
plot( alpha_result[,2]~rownames(alpha_result) , type="b" , bty="l" , xlab=expression(alpha) , ylab="Dunn" , col=rgb(240/255,128/255,128/255,0.8) , lwd=3 , pch=16 , ylim=c(0,1) )

####beta####
beta_result=read.table(file = "beta_result_matrix2ACC.txt",header = T,sep = "\t")
beta_result[,1]=-log10(beta_result[,1])
plot( beta_result[,1]~rownames(beta_result) , type="b" , bty="l" , xlab=expression(beta) , ylab=expression(-log[10](P)) , col=rgb(240/255,128/255,128/255,0.8) , lwd=3 , pch=16 , ylim=c(0,8) )
plot( beta_result[,2]~rownames(beta_result) , type="b" , bty="l" , xlab=expression(beta) , ylab="Dunn" , col=rgb(240/255,128/255,128/255,0.8) , lwd=3 , pch=16 , ylim=c(0,1) )



