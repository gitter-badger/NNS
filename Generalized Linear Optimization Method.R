#Step 1 Linear Optimization Method Generalized

#Objective.Function =
  
  

Lin_Opt <- function(A,objective.function = 1,
                    target_LPM=1,target_UPM=1,
                    degree_LPM=1,degree_UPM=1,
                    resolution=100,
                    minimum_weight_req=0,maximum_weight_req=1,
                    Objective = "max",method = "naive", 
                    shorts= "no", max_position = 1,
                    tol=1,pop=NULL){
  
    UPMs = 0L
    LPMs = 0L
    
    Current_set = numeric(0L)
    Dominated_set = numeric(0L)
    n = ncol(A)
 
    A_new= data.frame(matrix(ncol = n, nrow = length(A[,1])))
    A_new2= data.frame(matrix(ncol = n, nrow = resolution))
   
    weight = numeric(0L)
    w=numeric(0L)
    ww=numeric(0L)
    final_w = numeric(0L)
    z=numeric(0L)
    output = character(0L)
   
    if (shorts == "yes") {range <- -100} else {range <- 0}
  
for(k in 1:n){
  
  n = ncol(A)
  resolution = max(c(1000*(1/n),100))

  A_new= data.frame(matrix(ncol = n, nrow = length(A[,1])))
  A_new2= data.frame(matrix(ncol = n, nrow = resolution))
  #print(A[,1])
  #A_new= data.frame(matrix(ncol = n, nrow = 600))
  #A_new2= data.frame(matrix(ncol = n, nrow = resolution))
  
  weight = numeric(0L)
  w=numeric(0L)
  ww=numeric(0L)
  final_w = numeric(0L)
  
  output = character(0L)

  
     for (j in seq(range,100,((abs(range)+100)/resolution))){  
      for (i in 1:n){

          Key = A[,i]
          
          Others = A[,-i]
          
          w1<- j/100
          #Create the Portfolio Returns for each variable
                if (method == "naive") {
                
                w2<- (1-abs(w1))/(n-1)
              
                A_new[,i]=(w1*Key + rowSums(w2*Others))
                
                
                UPMs[i]=(UPM(degree_UPM,target_UPM,A_new[i]))
                LPMs[i]=(LPM(degree_LPM,target_LPM,A_new[i]))
                
              
                
                }
                
                
                
         
          
                if (method == "weighted") {
### Change weighting to match OBJECTIVE FUNCTION  
                 
                Z = apply(Others,2,function(Others) UPM(1,1,Others)/LPM(1,1,Others))/sum(apply(Others,2,function(Others) UPM(1,1,Others)/LPM(1,1,Others)))
                
                
                w2<- (1-abs(w1))
      
                A_new[i]=(rowSums(data.frame(w1*Key,w2*Z*Others)))
                }
          
                if (method == "random") {
### Change weighting to match OBJECTIVE FUNCTION  
                Z = apply(Others,2,function(Others) UPM(1,1,Others)/LPM(1,1,Others))/sum(apply(Others,2,function(Others) UPM(1,1,Others)/LPM(1,1,Others)))
               
                w2<- (1-abs(w1))
            
                A_new[i]=(rowSums(data.frame(w1*Key,w2*Z*Others)))
                }
      
      }
      ###IN J
      ###OBJECTIVE FUNCTION
     if (objective.function==1){
          objective = "max"
          A_new2[j+1,] =  UPMs/LPMs}
            
            #apply(A_new,2,function(A_new) UPM(degree_UPM,target_UPM,A_new)/LPM(degree_LPM,target_LPM,A_new))}
    
          #M/V
     if (objective.function==2){
          objective = "min"
          A_new2[j+1,] =   apply(A_new, 2, mean) / apply(A_new, 2, var)}
            
           
       
          # MIN VAR 
     if (objective.function==3){
          objective = "min"
          A_new2[j+1,] = apply(A_new, 2, var)}
            
            
            
         
          # GM / LPM
     if (objective.function==4){
          objective = "max"
          A_new2[j+1,] = apply(A_new, 2, function(A_new) exp(mean(log(A_new)))) / apply(A_new, 2, function(A_new) LPM(2,target_LPM,A_new))}
      
    }
    
     
    ###Select min or max observation from new columns and value weight it
 
  if (Objective == "max") {
  for (i in 1:n){ 
        if (shorts == "no"){ 
          if(which.max(na.omit(A_new2[,i]))==1){
        weight[i] = 0
        ww[i]=(which.max(na.omit(A_new2[,i])))
        } else {
          weight[i] = (which.max(na.omit(A_new2[,i])))*max(na.omit(A_new2[,i]))
          ww[i]=(which.max(na.omit(A_new2[,i])))}
          }
      
     
    
      ### NEED TO RESCALE THE W[i] -100
        if (shorts == "yes"){
        weight[i] = ((which.max(na.omit(A_new2[,i])-1)-(resolution/2))*max(na.omit(A_new2[,i])))
      print(which.max(na.omit(A_new2[,i]))-(resolution/2))  
          }
        
        }
   
   }
   
  
  if (Objective == "min") {
    for (i in 1:n){ 
      if (shorts == "no"){
        if(which.min(na.omit(A_new2[,i]))==1){
          weight[i] = 0} else {
        weight[i] = which.min(na.omit(A_new2[,i]))*min(na.omit(A_new2[,i]))}
        
      }
      
      ### NEED TO RESCALE THE W[i] -100
      if (shorts == "yes"){
        weight[i] = ((which.min(na.omit(A_new2[,i]))-(resolution/2))*min(na.omit(A_new2[,i])))
        print(which.min(na.omit(A_new2[,i]))-(resolution/2))  
      }
      
    }
    
  }
   
    ###Assign weights, if very low weight, discard variable via Dominated Set
    for (i in 1:n){
        if (shorts == "no"){
        w[i] = weight[i]/sum(weight)
    
            if (w[i]>minimum_weight_req) {w[i]=w[i]} 
            else {w[i]=0
                  Dominated_set[i] = i} }
      
      
      if (shorts == "yes"){
       
        w[i] = weight[i]/sum(abs(weight))
        
        if (abs(w[i])>minimum_weight_req && abs(w[i])< max_position ) {w[i]=w[i]} 
        else {w[i]=0
        Dominated_set[i] = i} 
        
        } #shorts = yes
      
      
        } #weighting assignment
  
      
      print(Dominated_set)
    
    ### If all variables kept, return output      
    if (length(na.omit(Dominated_set))==0){  
    
        ##For fully invested 0 cash 
        for (i in 1:n) {
        
          final_w[i] = w[i]/sum(abs(w))
     
          output[i] = print(paste(final_w[i],colnames(A[i]),sep="*"),quote = FALSE)
    
          }
    final_w.min = pmax(0,final_w-.5)
    final_w.max = pmin(1,final_w+.5)
   
   
    
    A <- as.matrix(A) 
    
    A.Lin.output = A
    
    PORTFOLIO_RETURNS=  rowSums(A %*% final_w )
   
   # print(PORTFOLIO_RETURNS)
    par(mfrow=c(1,2)) 
    hist(PORTFOLIO_RETURNS, main = "PORTFOLIO RETURNS")
    barplot(final_w, names.arg = colnames(A),las=2, main="Portfolio Weights")
    
    return(Sim.opt.combined(A,min.weight = final_w.min, 
                            max.weight =final_w.max ,
                           degree_UPM=degree_UPM, degree_LPM = degree_LPM, 
                            target_LPM=target_LPM, target_UPM=target_UPM,
                            objective.function=objective.function,
                            tol=tol,threshold=minimum_weight_req,
                            upper.threshold=maximum_weight_req,
                            pop=pop))
  
 #   return(c("UPM/LPM" = UPM(degree_UPM,target_UPM,PORTFOLIO_RETURNS)/LPM(degree_LPM,target_LPM,PORTFOLIO_RETURNS),
  # ("UPM/LPM_1" = UPM(1,1,PORTFOLIO_RETURNS)/LPM(1,1,PORTFOLIO_RETURNS)),exp(mean(log(PORTFOLIO_RETURNS)))/LPM(1,1,PORTFOLIO_RETURNS)))
  
   }
     
  
    ### If another variable Dominated, make a new matrix of remaining variables and repeat 
    else {
  
        A<- A[-na.omit(c(Dominated_set))]
        ### Reset the Dominated Set
        Dominated_set = numeric(0L)
        
  
      }
  
       
  } #K Loop
   
 
      
}   #END
