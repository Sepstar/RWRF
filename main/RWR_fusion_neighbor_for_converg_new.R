##20171104byhes
##本脚本是相似性网络融合的重启随机游走的函数

RWR_fusion_neighbor <- function(sim_list,iteration_max=1000,gama=0.7,neighbor_num=10,alpha=0.9,beta=0.9) 
{
  
  ########
  len_net=length(sim_list)#相似性网络的个数
  len_node=nrow(sim_list[[1]])#相似性网络中的节点数
  
  if (0)#相似性矩阵的归一化
  {
    for (i in 1:len_net)
    {
      for (m in 1:len_node)
      {
        for (n in 1:len_node)
        {
          if (m==n)
          {
            sim_list[[i]][m,n]=0.5
          }else
          {
            sim_list[[i]][m,n]=0.5*sim_list[[i]][m,n]/(sum(sim_list[[i]][m,])-sim_list[[i]][m,m])
          }
        }
      }
    }
  }
  
  
  if (0)
  {
    for (i in 1:len_net)
    {
      sim_list[[i]]=(sim_list[[i]])^2
    }
  }
  
  ####parameters####
  #alpha=0.5
  lamda=1/len_net
  #iteration_max=1000
  #gama=0.7
  #neighbor_num=10
  #beta=0.5
  ####output matrix####
  result_matrix=matrix(data = NA,nrow = len_net*len_node, ncol = len_net*len_node)
  rownames(result_matrix)=rep(rownames(sim_list[[1]]),len_net)
  colnames(result_matrix)=rownames(result_matrix)
  
  ####由最近邻构成的邻接矩阵####
  adj_list=NULL
  for (i in 1:len_net)
  {
    temp=sim_list[[i]]
    temp[which(temp!=0)]=0
    for (j in 1:len_node)
    {
      num=setdiff(which(rank(-sim_list[[i]][j,])<=neighbor_num+1),j)
      temp[j,num]=(1-beta)/neighbor_num
      temp[j,j]=beta#邻接矩阵中dj(net1)和dj(net2)之间的关联打分为beta，dj(net1)和dk(net2)（k!=j）之和为1-beta
      rm(num)
    }
    #temp=(temp+t(temp))/2
    adj_list[[i]]=temp
    rm(temp)
  }
  
  ####随机游走####
  end_clock_list=NULL
  start_time_all=Sys.time()
  for (i in 1:nrow(result_matrix))
  {
    #print(i)
    #start_time=Sys.time()
    index=ceiling(i/len_node)#判断此时循环到第几个网络
    #print(index)#用于调试，检查是否分类正确
    ####1设置初值####
    p0=rep(0,len_net*len_node)
    p0[i]=alpha
    temp=i%%len_node
    for (j in 1:len_net)
    {
      if (j!=index)
      {
        if (temp==0)
        {
          p0[(1+(j-1)*len_node):(j*len_node)]=(1-alpha)*(adj_list[[index]][len_node,])/(len_net-1)
        }else
        {
          p0[(1+(j-1)*len_node):(j*len_node)]=(1-alpha)*(adj_list[[index]][temp,])/(len_net-1)
        }
      }
    }
    #print(p0)#用于调试，检查初值是否设置正确
    
    ####2填写转移概率矩阵M（注意：迭代时要把M转置）####
    W_list=NULL
    for (m in 1:len_net)
    {
      for (n in 1:len_net)
      {
        if (m==n)#同一个相似性网络转移
        {
          temp=NULL
          for (j in 1:len_node)
          {
            temp=rbind(temp,lamda*sim_list[[m]][j,]/sum(sim_list[[m]][j,]))
          }
          W_list=c(W_list,list(temp))
          rm(temp)
        }else#不同的相似性网络之间转移m->n
        {
          #W_list=c(W_list,list(diag(x=lamda,nrow = len_node,ncol = len_node)))
          temp=NULL
          for (j in 1:len_node)
          {
            temp=rbind(temp,lamda*adj_list[[m]][j,]/sum(adj_list[[m]][j,]))
          }
          W_list=c(W_list,list(temp))
          rm(temp)
        }
      }
    }
    W=NULL
    k=1
    for (m in 1:len_net)
    {
      temp=NULL
      for (n in 1:len_net)
      {
        temp=cbind(temp,W_list[[k]])
        k=k+1
      }
      W=rbind(W,temp)
    }
    
    
    ####3迭代####
    temp=p0
    for (iteration_num in 1:iteration_max)
    {
      pt0=temp
      pt1=(1-gama)*t(W)%*%pt0+gama*p0
      temp=pt1
      end_clock=sum(abs(pt1-pt0))#迭代终止条件，L1范数（即该向量元素的绝对值之和）小于1e-10
      #print(end_clock/sum(pt0))
      if (end_clock<=1e-10)
      {
        #print(iteration_num)
        #print(end_clock)
        break
      }
    }
    end_clock_list=c(end_clock_list,end_clock/sum(pt0))
    rm(temp)
    rm(end_clock)
    pt1=pt1/sum(pt1)
    result_matrix[i,]=pt1
    #end_time=Sys.time()
    #print(end_time-start_time)
  }
  end_time_all=Sys.time()
  print(end_time_all-start_time_all)
  RWR_similarity=result_matrix[1:len_node,1:len_node]
  RWR_similarity[which(RWR_similarity!=0)]=0
  for (m in 1:len_net)
  {
    for (n in 1:len_net)
    {
      RWR_similarity=result_matrix[(1+(m-1)*len_node):(m*len_node),(1+(m-1)*len_node):(m*len_node)]+RWR_similarity
    }
  }
  RWR_similarity=RWR_similarity/(len_net)
  RWR_similarity=(RWR_similarity+t(RWR_similarity))/2
  
  if (0)#对结果矩阵进行行列归一化
  {
    for (i in 1:len_node)
    {
      for (j in 1:len_node)
      {
        RWR_similarity[i,j]=RWR_similarity[i,j]/(sum(RWR_similarity[i,])-RWR_similarity[i,i])
      }
    }
    diag(RWR_similarity)=1
  }
  
  if (1)#对结果进行最大值归一化
  {
    RWR_similarity=RWR_similarity/max(RWR_similarity)
    #diag(RWR_similarity)=1
  }
  
  return(c(list(RWR_similarity),list(mean(end_clock_list))))
}



