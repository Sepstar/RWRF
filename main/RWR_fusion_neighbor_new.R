RWR_fusion_neighbor <- function(sim_list,iteration_max=1000,gama=0.7,neighbor_num=10,alpha=0.9,beta=0.9) 
{
  
  ########
  len_net=length(sim_list)#Number of similarity networks
  len_node=nrow(sim_list[[1]])#Number of nodes in similarity network
  
  if (0)#Normalization of similarity matrix
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
  
  ####Adjacency matrix####
  adj_list=NULL
  for (i in 1:len_net)
  {
    temp=sim_list[[i]]
    temp[which(temp!=0)]=0
    for (j in 1:len_node)
    {
      num=setdiff(which(rank(-sim_list[[i]][j,])<=neighbor_num+1),j)
      temp[j,num]=(1-beta)/neighbor_num
      temp[j,j]=beta
      rm(num)
    }
    #temp=(temp+t(temp))/2
    adj_list[[i]]=temp
    rm(temp)
  }
  
  ####random walk####
  start_time_all=Sys.time()
  result_matrix = foreach (i = 1:nrow(result_matrix),.combine = rbind) %dopar%
  {
    #print(i)
    #start_time=Sys.time()
    index=ceiling(i/len_node)#Determine the number of networks to loop to at this time
    ####1Set initial value####
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
    
    ####2Fill in the transition probability matrix M####
    W_list=NULL
    for (m in 1:len_net)
    {
      for (n in 1:len_net)
      {
        if (m==n)# in the same similarity network
        {
          temp=NULL
          for (j in 1:len_node)
          {
            temp=rbind(temp,lamda*sim_list[[m]][j,]/sum(sim_list[[m]][j,]))
          }
          W_list=c(W_list,list(temp))
          rm(temp)
        }else#different similarity network
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
    
    ####3Iterate####
    temp=p0
    for (iteration_num in 1:iteration_max)
    {
      pt0=temp
      pt1=(1-gama)*t(W)%*%pt0+gama*p0
      temp=pt1
      end_clock=sum(abs(pt1-pt0))
      if (end_clock<=1e-10)
      {
        #print(iteration_num)
        #print(end_clock)
        break
      }
    }
    rm(temp)
    rm(end_clock)
    pt1=pt1/sum(pt1)
    result_matrix[i,]=pt1
    t(pt1)
    #end_time=Sys.time()
    #print(end_time-start_time)
  }
  end_time_all=Sys.time()
  print(end_time_all-start_time_all)
  rownames(result_matrix)=rep(rownames(sim_list[[1]]),len_net)
  colnames(result_matrix)=rownames(result_matrix)
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
  
  RWR_similarity=RWR_similarity/max(RWR_similarity)
  #diag(RWR_similarity)=1
  
  return(RWR_similarity)
}



