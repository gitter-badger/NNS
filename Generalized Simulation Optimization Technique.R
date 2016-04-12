Sim_opt<- function(A,pop = ncol(A)*100,threshold = 0){
  
  n<- ncol(A)
  weight = numeric(0L)
  weight.pct = numeric(0L)
  values = numeric(0L)
  Dominated_set = numeric(0L)
  
  minimum_weight = numeric(0L)
  maximum_weight = numeric(0L)
  
  minimum_weight = rep(0,n)
  maximum_weight = rep(100,n)
  
  mean_weight = rep(100,n)
  
  final_weight = numeric(0L)
  
  max.return.weights = data.frame(matrix(ncol = n+1, nrow = (100*n)))
  max.return.weights.iterations = data.frame(matrix(ncol = n+1, nrow = length(2*n)))
 l= 1L
 while (l > 0){
#for (l in 1:(2*n)){
  n<- ncol(A)
  for (k in 1:(1*n)){
    n<- ncol(A)/l
    ### Perhaps for (j in (min(c(pop,max(c(100*n,100))))) )
    for (j in 1:(min(c(pop,max(c(100*n,100)))))){
      
          n<- ncol(A)
    
          for (i in 1:n){ weight[i] = sample(minimum_weight[i]:maximum_weight[i],1,replace = TRUE)} 
     
      
          for (i in 1:length(weight)){ weight.pct[i] = (weight[i]/sum(weight))  }
    
         
          
          AA<- as.matrix(A)  
          
          port = rowSums(AA %*% diag(weight.pct))
          UPMLPM = UPM(1,1,port)/LPM(1,1,port)
    
          max.return.weights[j,] = c(weight,UPMLPM)
          
          }
     print(c(j,k,l))
          highest = (which.max(max.return.weights[,(n+1)]))
  
  
          max.return.weights.iterations[k,] = c(max.return.weights[highest,])
  
     }#K
    
    
  #print(max.return.weights.iterations) 
  
  
  for (i in 1:n){
      minimum_weight[i] = min(max.return.weights.iterations[,i])
      maximum_weight[i] = max(max.return.weights.iterations[,i])
      mean_weight[i] = mean(c(minimum_weight[i],maximum_weight[i]))
    
     if (minimum_weight[i] <= 1) {Dominated_set[i] = i} else {Dominated_set[i] = 0}
     #Min weight Threshold...
        if (minimum_weight[i]<= (threshold*100)) {
            minimum_weight[i] = 0 
            maximum_weight[i]= 0} else {minimum_weight[i] =minimum_weight[i]
            }
       }
  
        
      
          print(matrix(c(colnames(A), minimum_weight,maximum_weight,mean_weight/sum(mean_weight)),nrow=ncol(A),ncol = 4), quote = FALSE)
    
          A <- as.matrix(A) 
          
          PORTFOLIO_RETURNS=  rowSums(A %*% diag(mean_weight/sum(mean_weight)) )
           
          par(mfrow=c(1,2)) 
          hist(PORTFOLIO_RETURNS, main = "PORTFOLIO RETURNS")
          barplot(mean_weight/sum(mean_weight), 
          names.arg = colnames(A[,!colnames(A) %in% (Dominated_set)]  ),
                  las=2, main="Portfolio Weights")
          print(c("UPM/LPM" = UPM(1,1,PORTFOLIO_RETURNS)/LPM(1,1,PORTFOLIO_RETURNS)))
          
         l = l + 1L 
         
         if (mean(maximum_weight)-mean(minimum_weight)<1) {break}
   
   
   
         

  } #L 
      
}

 
  