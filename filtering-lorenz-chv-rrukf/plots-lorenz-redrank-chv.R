source("scores.r")
cat("\tPlease make sure to run run-filtering-lorenz-chv-redrank.R before running this file. \n")

#####################################################################################-
############# Average Log Score plotting (reduced rank) #############################-
#####################################################################################-
data.model.list = c("gauss", "gamma"
); N = 768
dtvec = 0.0005#c(0.005, 0.0005); 
b = 0.2
seed.values = 1988 + (1:5 - 1)*100
ScoreFun = calculateLSs

for(j in 1:length(dtvec)){
  
  dt = dtvec[j]
  pdf(paste0("filtering-lorenz-chv-rrukf/b", b, "_LogSc_avg_rrukf_", dt, "_dt", ".pdf"), height = 4, width = 10)
  par(mfrow = c(1, 2), mar = c(4, 3, 1, 1))
  for(data.model in data.model.list){
    
    RMSE.score.avg = array(dim = c(length(seed.values), 40, 3))
    for(i in 1:5){
      
      if(data.model != "gauss"){
        
        load(paste0("filtering-lorenz-chv-rrukf/b", b, "_lorenzrunredrankUKF_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
        load(paste0("filtering-lorenz-chv-rrukf/b", b, "_lorenzrunours_N_", N, "_dt_", dt, "_", data.model, "_adderrors.Rdata"))
      }else{
        
        load(paste0("filtering-lorenz-chv-rrukf/b", b, "_lorenzrunredrankUKF_N_", N, "_dt_", dt, ".Rdata"))
        load(paste0("filtering-lorenz-chv-rrukf/b", b, "_lorenzrunours_N_", N, "_dt_", dt, "_adderrors.Rdata"))
      }
      RMSE.score.avg[i, , ] = ScoreFun(result.lorenz.redrankukf[[i]]$predsMRA, 
                                       result.lorenz.ours[[i]]$predsMRA, 
                                       result.lorenz.ours[[i]]$XY$x)
      # RMSE.score.2 = ScoreFun(result.lorenz.marcin[[i]]$predsMRA, 
      #                               result.lorenz.marcin.2[[i]]$predsMRA, 
      #                               result.lorenz.marcin[[i]]$XY$x)
      # RMSE.score = cbind(RMSE.score, RMSE.score.2[ , 3])
      duplicate = length(result.lorenz.ours[[i]]$predsMRA)
      if(duplicate < length(result.lorenz.redrankukf[[i]]$XY$x)){
        
        RMSE.score.avg[i, duplicate, 3] = NA 
      }
    }
    
    RMSE.score = apply(RMSE.score.avg, c(2, 3), mean)
    matplot(x = RMSE.score[ , 1], y = RMSE.score[ , -1],
            col = c("blue", "red"#, "black"
            ), type = "l", lty = rep(1, 2),
            xlab = "time", ylab = " ")
    #title(main = paste0("data.model = ", data.model))
    if(data.model == "gauss"){
      
      legend("topleft", legend = c("RRUKF", "CHV"#, "Marcin.2"
      ), 
      col = c("blue", "red"#, "black"
      ), lty = rep(1, 2))
      #title(ylab = paste0("dt = ", dt))
    }
  }
  dev.off()
}