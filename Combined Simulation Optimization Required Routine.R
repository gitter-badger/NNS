Sim.opt.combined<- function(A,pop = ncol(A)*100,
                            threshold = minimum_weight_req,
                            upper.threshold = maximum_weight_req,
                            min.weight=NULL, 
                            max.weight = NULL,
                            degree_UPM=degree_UPM,degree_LPM=degree_LPM,
                            target_UPM=target_UPM,target_LPM=target_LPM,
                            objective.function, tol=tol){
  
  n<- ncol(A)
  weight = numeric(0L)
  weight.pct = numeric(0L)
  values = numeric(0L)
  Dominated_set = numeric(0L)
  
  minimum_weight = numeric(0L)
  maximum_weight = numeric(0L)
  
  #minimum_weight = as.integer(min.weight*100)
  #maximum_weight = as.integer(max.weight*100)
  minimum_weight = rep(0,n)
  maximum_weight = rep(100,n)

  mean_weight = rep(100,n)
  
  final_weight = numeric(0L)
  
  max.return.weights = data.frame(matrix(ncol = n+1, nrow = (100*n)))
  max.return.weights.iterations = data.frame(matrix(ncol = n+1, nrow = length(2*n)))
  
  
  l= 1L
  while (l > 0){
    
   
  
    for (k in 1:(1*n)){
    n<- ncol(A)/l
      for (j in 1:(min(c(pop,max(c(100*n,100)))))){
        
        n<- ncol(A)
        for (i in 1:n){ weight[i] = sample(minimum_weight[i]:maximum_weight[i],1,replace = TRUE)} 
        
        
        for (i in 1:length(weight)){ weight.pct[i] = (weight[i]/sum(weight))  
        
        
        }
        
        for (i in 1:n){
     if(weight.pct[i]<threshold) {weight.pct[i]=0} }
        
        
        
        #AA<- as.matrix(A)  
        
        port = rowSums(A %*% weight.pct)
        UPMLPM = UPM(degree_UPM,target_UPM,port)/LPM(degree_LPM,target_LPM,port)
        Var = 1/var(port)
        MV = mean(port)/var(port)
        GM.SV = exp(mean(log(port)))/LPM(2,target_LPM,port)
        
        if(objective.function==1) {max.return.weights[j,] = c(weight,UPMLPM)}
        if(objective.function==2) {max.return.weights[j,] = c(weight,MV)}
        if(objective.function==3){max.return.weights[j,] = c(weight,Var)}
        if(objective.function==4){max.return.weights[j,] = c(weight,GM.SV)}
        
        
        
      }
      print(c(j,k,l))
      highest = (which.max(max.return.weights[,(n+1)]))
      
      
      max.return.weights.iterations[k,] = c(max.return.weights[highest,])
      
    }#K
    
    
    
    
    for (i in 1:n){
      minimum_weight[i] = min(max.return.weights.iterations[,i])
      maximum_weight[i] = max(max.return.weights.iterations[,i])
      mean_weight[i] = mean(c(minimum_weight[i],maximum_weight[i]))
      if(mean_weight[i]<threshold*100){mean_weight[i]=0}                                                 
      
      if (maximum_weight[i] <= 1) {Dominated_set[i] = i} 
      
      #else {Dominated_set[i] = 0}
      #Min weight Threshold...
      if (minimum_weight[i]<= (threshold*100)) {
        minimum_weight[i] = 0 
        maximum_weight[i]= 0} else {minimum_weight[i] =minimum_weight[i]
        }
    }
     
    mean.weights.vector = mean_weight/sum(mean_weight)
    for(i in 1:n){
    if(mean_weight[i]/sum(mean_weight)<threshold){mean_weight[i]=0}  
    #if(mean.weights.vector[i]>upper.threshold){mean.weights.vector[i]=upper.threshold}
    }
    
    
    
    
    
    print(mean_weight/sum(mean_weight))
    
    PORTFOLIO_RETURNS=  rowSums(A %*% mean_weight/sum(mean_weight)) 

    remaining.vars = (matrix(c(colnames(A), minimum_weight,maximum_weight,mean_weight/sum(mean_weight)),nrow=ncol(A),ncol = 4))
    
    
    
    bar.plot.remaining.vars = mean_weight/sum(mean_weight)
    
    bar.plot.remaining.vars = bar.plot.remaining.vars[bar.plot.remaining.vars>0]
    
    var.names = subset(remaining.vars, remaining.vars[,4]>0)
    
    
    par(mfrow=c(1,2)) 
    hist(PORTFOLIO_RETURNS, main = "PORTFOLIO RETURNS")
    barplot(bar.plot.remaining.vars, 
            names.arg = c(var.names[,1]),
            las=2, main="Portfolio Weights")
    
    print(c("UPM/LPM" = UPM(degree_UPM,target_UPM,PORTFOLIO_RETURNS)/LPM(degree_LPM,target_LPM,PORTFOLIO_RETURNS)))
    print(mean(PORTFOLIO_RETURNS)/var(PORTFOLIO_RETURNS))
    print(var(PORTFOLIO_RETURNS))
    print(exp(mean(log(PORTFOLIO_RETURNS)))/LPM(2,target_LPM,PORTFOLIO_RETURNS))
    
    l = l + 1L 
    print(mean(maximum_weight)-mean(minimum_weight))
    if (mean(maximum_weight)-mean(minimum_weight)<tol)  
      
      return(SSD.MV.portfolio<<-PORTFOLIO_RETURNS)
      
      #return(print(matrix(c(colnames(A), minimum_weight,maximum_weight,mean_weight/sum(mean_weight)),nrow=ncol(A),ncol = 4), quote = FALSE))
    
    
     
    
  } #L 
  
  
  
}


