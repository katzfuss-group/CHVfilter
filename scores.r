RMSE = function(x, y){
  N = length(x)
  res = sqrt(sum((x-y)**2)/N)
  return(res)
}



calculateRMSPE = function(predsMRAstandard, predsMRAnoevolerror, 
                          #predsMRAnoevolerrornoapprox, 
                          x){

  Tmax = length(x)
  length.one = length(predsMRAstandard)
  length.two = length(predsMRAnoevolerror)
  RMSPEs = c()
  
  for(t in 1:Tmax){
    if(t <= length.one){
      
      MRA = RMSE(x[[t]], predsMRAstandard[[t]]$state)
       
    }else{
      MRA = NA
    }
    
    if(t <= length.two){
      
      LR = RMSE(x[[t]], predsMRAnoevolerror[[t]]$state)
    }else{
      LR = NA
    }
     
    #LR2 = RMSE(x[[t]], predsMRAnoevolerrornoapprox[[t]]$state)
    RMSPEs.t = c(t, MRA, LR#, LR2
                 )
    RMSPEs = c(RMSPEs, RMSPEs.t)
  }
  RMSPEs = matrix( RMSPEs, ncol=3, byrow=TRUE)
  colnames(RMSPEs) = c("time", "Standard", "NoErrorApp"#, "NoErrorNoApp"
                       )
  return( RMSPEs )
}



calculateRRMSPE = function(predsMRAstandard, predsMRAnoevolerror, predsE, x){
  
  Tmax = length(x)
  RRMSPEs = c()
  
  for(t in 1:Tmax){
    MRA = RMSE(x[[t]], predsMRAstandard[[t]]$state)
    LR = RMSE(x[[t]], predsMRAnoevolerror[[t]]$state)
    E  = RMSE(x[[t]], predsE[[t]]$state)
    RRMSPEs.t = c(t, MRA/E, LR/E)
    RRMSPEs = c(RRMSPEs, RRMSPEs.t)
  }
  RRMSPEs = matrix( RRMSPEs, ncol=3, byrow=TRUE)
  colnames(RRMSPEs) = c("time", "Standard", "NoError")
  return( RRMSPEs )
}



LogS = function(preds, y){
  
  Tmax = length(y)
  LS = rep(0, Tmax)
  length.one = length(preds)
  for( t in 1:Tmax ){
    if(t <= length.one){
      LS[t] = tryCatch(-dposterior(y[[t]], preds[[t]]), error = function(e) NA) 
    }else{
      LS[t] = NA
    }
  }
  return(LS)
}


dposterior = function(y, pred){
  
  mu = pred$state
  W= pred$W
  # VL calculates V, Laplace calculates W                           
  if("V" %in% names(pred)){
    #det_term = sum(log(diag(pred$V)))
    det_term = 0
    for(i in 1:((dim(pred$V))[1]))
    {
      det_term = det_term + log((pred$V)[i, i])
    }
  }else{
    V=t(chol(Matrix::forceSymmetric(pred$W)))
    #V = t(chol((W + t(W))/2))                                      
    det_term = sum(log(diag(V)))
    #det_term =log(det(W))/2                                        
  }
  quad_term = -Matrix::t(y-mu) %*% W %*% (y-mu)/2
  pi_term = -length(y)/2*log(2*pi)
  # value is of class "dgeMatrix"                                   
  return((quad_term+det_term+pi_term)[1,1])
}


calculateLSs = function(predsMRAstandard, predsMRAnoevolerror, x){
  LS.MRA = LogS(predsMRAstandard, x)
  LS.LR = LogS(predsMRAnoevolerror, x)
  #LS.E = LogS(predsE, x)
  time = seq(1:length(LS.MRA))
  results = cbind(time, LS.MRA#-LS.E
                  , LS.LR#-LS.E
                  )
  colnames(results) = c("time", "Standard", "NoError")
  return(results)
}
