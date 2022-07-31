rm(list = ls()); #foreach::registerDoSEQ()
library(Rcpp)
source("aux_funs_general.R")
Rfast::sourceR(path = "GPvecchia//R//")
#source("C://PhD Research//Self//UMRF//Resource codes//vecchiaFilter//aux-functions.R")
source("lorenz-simulations-functions//aux-functions-Lorenz.r")

library(VEnKF)
library(Matrix)
#library(GPvecchia)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)




######### set parameters #########
#set.seed(1988)
#n = 32
m = 20
frac.obs = 0.1
Tmax = 40

## evolution values ##
Force = 10
K = 35
dtvec = 0.0005
M = 30
b = 0.2

max.iter = 5#getDoParWorkers()



## covariance function
sig2 = 0.0; range = .15; smooth = 0.5;
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)


## likelihood settings
me.var = 0.1;
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  if (!(args[1] %in% c("gauss", "poisson", "logistic", "gamma"))) {
    stop("One of the models has to be passed as argument")
  } else {
    data.model = args[1]
  }
} else {
  data.model = "gauss"
}
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)


## generate grid of pred.locs
M.multires = 7; r = c(6, rep(3, M.multires)); J = c(2, rep(2, M.multires-1))
#load("auxoutput_circular.Rdata")
#grid.oneside = seq(0,1,length = round(n))
locsettings = generatePoints_v2(M.multires, r, J)
attach(locsettings)
n = length(locs)

## set initial state
#cat("Loading the moments of the long-run Lorenz\n")

#Sig0ordered = Sig0
#muordered = mu
####x0 = b*getX0(n, Force, K, dt)
#x0 = t(chol(Sig0)) %*% matrix(rnorm(n), ncol=1) + mu
#Sigt = sig2*Sig0
Sigt = diag(rep(10^-100, n)) #covfun(locs)


## define Vecchia approximation
cat("Calculating the approximations\n")

vecchia_obj.circ = vecchia_specify_modified(matrix(locsettings$locsOrderedPoints), 
                                            ordering = 'none', 
                         conditioning = 'mra', 
                         mra.options = list(M = M.multires, J = J, r = r)) 
approximations = list(mra = vecchia_obj.circ)
PointerData = list()
PointerData$DataComp = extractPointerData(M = M.multires, J = J, r = r, action = "compress")
PointerData$DataDeomp = extractPointerData(M = M.multires, J = J, r = r, action = "decompress")

data.model.list = c("gauss", "gamma")
for(data.model.dup in data.model.list){
  
  data.model = data.model.dup
  lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)
  
  
  for(i in 1:length(dtvec)){
    
    dt = dtvec[i]
    cl = parallel::makePSOCKcluster(5)
    registerDoParallel(cl)
    result.lorenz.ours = foreach( iter=1:max.iter) %dopar% {
      #for (iter in 1:max.iter) {
      library(VEnKF);library(GPvecchia)
      
      seed.value = 1988 + 100 * (iter - 1)
      
      moments = getLRMuCovariance(n, Force, dt, K, seed = seed.value)
      Sig0 = ((b**2)*moments[["Sigma"]] + diag(1e-10, n))
      mu = (b*moments[["mu"]])
      
      cat("Simulating data\n")
      set.seed(seed.value)
      x0 = (t(chol(Sig0)) %*% matrix(rnorm(n), ncol=1) + mu)
      XY = simulate.xy(x0, 
                       evolFun, Q = NULL, frac.obs, lik.params, Tmax, sig2 = 0, seed = seed.value)
      
      Sig0 = Sig0[locsettings$MatchLocationOrdering, 
                  locsettings$MatchLocationOrdering]
      mu = mu[locsettings$MatchLocationOrdering]
      
      for(i in 1:Tmax)
      {
        XY$x[[i]] = (XY$x[[i]])[locsettings$MatchLocationOrdering, , drop = FALSE]
        XY$y[[i]] = (XY$y[[i]])[locsettings$MatchLocationOrdering]
      }
      
      cat(paste("iteration: ", iter, " started.", "\n"))
      
      
      start = proc.time()
      predsMRAnoevolerror = 
        filterLorenz_v5('mra', XY, mod = "noevolerror", PointerData = PointerData, DataSetup.init = locsettings,
                        mu.tt = mu, Sig0 = Sig0, r = r)
      
      d2 = as.numeric(proc.time() - start)
      data.list = list()
      data.list$XY = XY
      #data.list$predsMRAstandard = predsMRAstandard
      #data.list$d = d
      data.list$predsMRA = predsMRAnoevolerror
      data.list$d2 = d2
      #data.list$predsMRAnoevolerrornoapprox = predsMRAnoevolerrornoapprox
      #data.list$d3 = d3
      cat(paste("iteration: ", iter, " ended.", "\n"))
      return(data.list)
    }
    if(data.model == "gauss"){
      save(locsettings, result.lorenz.ours, file = paste0("filtering-lorenz-chv-ekvl/b", b, "_lorenzrunours_N_", n, "_dt_", dt, ".Rdata")) 
    }else{
      save(locsettings, result.lorenz.ours, file = paste0("filtering-lorenz-chv-ekvl/b", b, "_lorenzrunours_N_", n, "_dt_", dt, "_", data.model, ".Rdata")) 
    }
    
    parallel::stopCluster(cl)
    registerDoSEQ()
    rm(cl)
    cat(paste0("\n\t\ti = ", i))
  }
  
  cat(paste0("\nCompleted our code for data.model = ", data.model))
  
  for(i in 1:length(dtvec)){
    
    dt = dtvec[i]
    cl = parallel::makePSOCKcluster(5)
    registerDoParallel(cl)
    result.lorenz.marcin = foreach( iter=1:max.iter) %dopar% {
      
      library(VEnKF); library(GPvecchia); 
      seed.value = 1988 + 100 * (iter - 1)
      moments = getLRMuCovariance(n, Force, dt, K, seed = seed.value)
      Sig0 = ((b**2)*moments[["Sigma"]] + diag(1e-8, n))
      mu = (b*moments[["mu"]])
      
      cat("Simulating data\n")
      set.seed(seed.value)
      x0 = (t(chol(Sig0)) %*% matrix(rnorm(n), ncol=1) + mu)
      XY = simulate.xy(x0, 
                       evolFun, Q = NULL, frac.obs, lik.params, Tmax, sig2 = 0, seed = seed.value)
      
      Sig0 = Sig0[locsettings$MatchLocationOrdering, 
                  locsettings$MatchLocationOrdering]
      mu = mu[locsettings$MatchLocationOrdering]
      
      for(i in 1:Tmax)
      {
        XY$x[[i]] = (XY$x[[i]])[locsettings$MatchLocationOrdering, , drop = FALSE]
        XY$y[[i]] = (XY$y[[i]])[locsettings$MatchLocationOrdering]
      }
      
      #cat(paste("iteration: ", iter, ", MRA", "\n", sep = ""))
      cat(paste("iteration: ", iter, " started.", "\n"))
      approximations = list(mra = vecchia_obj.circ)
      start = proc.time()
      predsMRAstandard = 
        filterLorenz(approximations, 'mra', XY, mu, Sig0, Sigt,
                     data.model, lik.params, Force, K, dt, M, DataSetup.init = locsettings)
      
      d = as.numeric(proc.time() - start)
      
      
      data.list = list()
      data.list$XY = XY
      data.list$predsMRA = predsMRAstandard
      data.list$d = d
      cat(paste("iteration: ", iter, " ended.", "\n"))
      return(data.list)
    }
    if(data.model == "gauss"){
      save(locsettings, result.lorenz.marcin, file = paste0("filtering-lorenz-chv-ekvl/b", b, "_lorenzrunmarcin_N_", n, "_dt_", dt, ".Rdata")) 
    }else{
      save(locsettings, result.lorenz.marcin, file = paste0("filtering-lorenz-chv-ekvl/b", b, "_lorenzrunmarcin_N_", n, "_dt_", dt, "_", data.model, ".Rdata")) 
    }
    
    parallel::stopCluster(cl)
    registerDoSEQ()
    rm(cl)
    cat(paste0("\n\t\ti = ", i))
  }
  
  cat(paste0("\nCompleted Marcin's code for data.model = ", data.model))
}
