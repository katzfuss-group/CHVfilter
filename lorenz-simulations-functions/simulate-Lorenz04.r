rm(list=ls())
library(VEnKF)
library(rootSolve)

    
#N = 32
#Force = 10
#K = 10
#dt = 0.005
#M = 5
#b = 0.1
#Force = 10
#dt = 0.005
#M = 40
#K = 32
N = 768#1280
Force = 10
K = 35#50
dtvec = 0.0005#c(0.005, 0.0005)
M = 30
b = 0.2
Tmax = 80
#seed = 1988

seed.values = 1988 + (1:5 - 1) * 100

for(i in 1:length(dtvec)){

  dt = dtvec[i]
for(seed in seed.values){
fileName.init = paste("lorenz-simulations-functions/init_Lorenz04_N", N, "F", Force, "dt", dt, "K", K, "seed", seed, sep = "_")
fileName.all = paste("lorenz-simulations-functions/Lorenz04_N", N, "F", Force, "dt", dt, "K", K, "seed", seed, sep = "_")





#### generate solutions ####
generateInit = T

if (generateInit) {
  set.seed(seed)
  X0 = rnorm(N)  
} else {
  X0 = scan(fileName.init, quiet = TRUE)  
}
X = Lorenz04M2Sim(X0, Force, K, dt, M, iter = Tmax, burn = 0, order = 4)
#X1 = Lorenz04M2Sim(X0, Force, K, dt, M, iter=Tmax, burn=0, order=1)

#browser()

if ( generateInit ) {
  Xlast = X[,Tmax]
  write(Xlast, file = fileName.init)  
} else {
  write(X, file = fileName.all)
}


generateInit = F

if (generateInit) {
  set.seed(seed)
  X0 = rnorm(N)  
} else {
  X0 = scan(fileName.init, quiet = TRUE)  
}
X = Lorenz04M2Sim(X0, Force, K, dt, M, iter = Tmax, burn = 0, order = 4)
#X1 = Lorenz04M2Sim(X0, Force, K, dt, M, iter=Tmax, burn=0, order=1)

#browser()

if ( generateInit ) {
  Xlast = X[,Tmax]
  write(Xlast, file = fileName.init)  
} else {
  write(X, file = fileName.all)
}
cat(paste("\tSeed = ",seed, ", i = ", i, "\n", sep=""))
}
}