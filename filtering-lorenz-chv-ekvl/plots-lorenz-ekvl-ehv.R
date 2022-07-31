source("scores.r")
cat("\tPlease make sure to run run-filtering-lorenz-chv-ekvl.R before running this file. \n")
#state plotting (Lorenz Curve)---------------
#load("lorenzrunstandardnoevolerror.Rdata")
seed.values = 1988 + (1:5 - 1) * 100; N = 768
dtvec = 0.0005#c(0.005, 0.0005)
data.model = "gauss"
for(i in 1:length(dtvec))
{
  
  dt = dtvec[i]
  if(data.model != "gauss"){
    
    #load(paste0("lorenzrunmarcin_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
    #load(paste0("lorenzrunmarcin_mod_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
    load(paste0("filtering-lorenz-chv-ekvl/b0.2_lorenzrunours_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
  }else{
    
    #load(paste0("lorenzrunmarcin_N_", N, "_dt_", dt, ".Rdata"))
    #load(paste0("lorenzrunmarcin_mod_N_", N, "_dt_", dt, ".Rdata"))
    load(paste0("filtering-lorenz-chv-ekvl/b0.2_lorenzrunours_N_", N, "_dt_", dt, ".Rdata"))
  }
  
  attach(locsettings)
  #Tmax = length(result.lorenz.marcin[[1]]$XY$x)
  #list.len.marcin = list.len.marcin.2 = list.len.ours = c(); 
  #max.iter = length(result.lorenz.marcin);
  #for(j in 1:max.iter) {list.len.marcin = append(list.len.marcin, 
  #              length(result.lorenz.marcin[[j]]$predsMRA))}
  #for(j in 1:max.iter) {list.len.marcin.2 = append(list.len.marcin.2, 
  #              length(result.lorenz.marcin.2[[j]]$predsMRA))}
  #for(j in 1:max.iter) {list.len.ours = append(list.len.ours, 
  #              length(result.lorenz.ours[[j]]$predsMRA))}
  #for(j in 1:max.iter){
  #  
  #  Tmax.mod = max(c(list.len.marcin[i], list.len.ours[i]))
  #  if(isTRUE(data.model == "gauss"))
  #  {
  #    pdf(paste0("lorenzStatePlot_N_", N, "_dt_", dt, "_seed_", seed.values[j], 
  #               ".pdf"), height = 8, width = 10)
  #  }else{
  #    pdf(paste0("lorenzStatePlot_N_", N, "_dt_", dt, "_", data.model, "_seed_", seed.values[j], 
  #               ".pdf"), height = 8, width = 10)
  #  }
  #  
  #  par(mfrow = c(5, 4), mar = c(4, 4, 1, 1))
  #  for(t in 1:Tmax.mod){
  #    
  #    xmat = locs
  #    ymat = cbind(result.lorenz.ours[[j]]$XY$x[[t]])
  #    col = "green"; ltype = 1
  #    if(isTRUE(list.len.ours[j] > 0) && t <= list.len.ours[j]){
  #      
  #      ymat = cbind(ymat, result.lorenz.ours[[j]]$predsMRA[[t]]$state)
  #      col = c(col, "red"); ltype = c(ltype, 1)
  #    }
  #    if(isTRUE(list.len.marcin[j] > 0) && t <= list.len.marcin[j]){
  #      
  #      ymat = cbind(ymat, result.lorenz.marcin[[j]]$predsMRA[[t]]$state)
  #      col = c(col, "blue"); ltype = c(ltype, 1)
  #    }
  #    if(isTRUE(list.len.marcin.2[j] > 0) && t <= list.len.marcin.2[j]){
  #      
  #      ymat = cbind(ymat, result.lorenz.marcin.2[[j]]$predsMRA[[t]]$state)
  #      col = c(col, "black"); ltype = c(ltype, 1)
  #    }
  #    ymat = ymat[locsordOrdering, ]
  #    matplot(x = xmat, y = ymat, lty = ltype, col = col, main = "",
  #            xlab = paste("T =", t), ylab = "", type = "l", lwd = 0.2)
  #  }
  #  dev.off()
  #}
  detach(locsettings)
}

num = 1988#which(is.list.vec == TRUE)
Tmax = length(result.lorenz.ours[[1]]$predsMRAnoevolerror)

attach(locsettings); i = 1
#for(i in num){
t = 1
pdf(paste0("filtering-lorenz-chv-ekvl/lornenStatePlot_N_", length(locsettings$locsOrderedPoints), 
           "_T_", t,"_seed_", num,  ".pdf"), height = 8, width = 10)
par(mfrow = c(1, 1), mar = c(3, 3, 1, 1))
#for(t in 1:Tmax){

matplot(x = sort(locsOrderedPoints), 
        y = cbind((result.lorenz.ours[[i]]$XY$x[[t]])[locsettings$locsordOrdering], 
                  #(result.lorenz.marcin[[i]]$predsMRAstandard[[t]]$state)[locsettings$locsordOrdering], 
                  (result.lorenz.ours[[i]]$predsMRA[[t]]$state)[locsettings$locsordOrdering]),
        col = c("green"#, "blue"
                , "red"
        ), xlab = ""#paste("t =", t)
        , ylab = "", type = "l", lty = 1, lwd = 0.3)
#points(x = sort(locsOrderedPoints), 
#       y = (result.lorenz.ours[[i]]$XY$x[[t]])[locsettings$locsordOrdering], 
#       col = "black", pch = 4, cex = 0.5)
#}
dev.off() 


t = 3
pdf(paste0("filtering-lorenz-chv-ekvl/lornenStatePlot_N_", length(locsettings$locsOrderedPoints), 
           "_T_", t,"_seed_", num,  ".pdf"), height = 8, width = 10)
par(mfrow = c(1, 1), mar = c(3, 3, 1, 1))
#for(t in 1:Tmax){

matplot(x = sort(locsOrderedPoints), 
        y = cbind((result.lorenz.ours[[i]]$XY$x[[t]])[locsettings$locsordOrdering], 
                  #(result.lorenz.marcin[[i]]$predsMRAstandard[[t]]$state)[locsettings$locsordOrdering], 
                  (result.lorenz.ours[[i]]$predsMRA[[t]]$state)[locsettings$locsordOrdering]),
        col = c("green"#, "blue"
                , "red"
        ), xlab = ""#paste("t =", t)
        , ylab = "", type = "l", lty = 1, lwd = 0.3)
#points(x = sort(locsOrderedPoints), 
#       y = (result.lorenz.ours[[i]]$XY$x[[t]])[locsettings$locsordOrdering], 
#       col = "black", pch = 4, cex = 0.5)
#}
dev.off()


t = 30
pdf(paste0("filtering-lorenz-chv-ekvl/lornenStatePlot_N_", length(locsettings$locsOrderedPoints), 
           "_T_", t,"_seed_", num,  ".pdf"), height = 8, width = 10)
par(mfrow = c(1, 1), mar = c(3, 3, 1, 1))
#for(t in 1:Tmax){

matplot(x = sort(locsOrderedPoints), 
        y = cbind((result.lorenz.ours[[i]]$XY$x[[t]])[locsettings$locsordOrdering], 
                  #(result.lorenz.marcin[[i]]$predsMRAstandard[[t]]$state)[locsettings$locsordOrdering], 
                  (result.lorenz.ours[[i]]$predsMRA[[t]]$state)[locsettings$locsordOrdering]),
        col = c("green"#, "blue"
                , "red"
        ), xlab = ""#paste("t =", t)
        , ylab = "", type = "l", lty = 1, lwd = 0.3)
#points(x = sort(locsOrderedPoints), 
#       y = (result.lorenz.ours[[i]]$XY$x[[t]])[locsettings$locsordOrdering], 
#       col = "black", pch = 4, cex = 0.5)
#}
dev.off()
#}
detach(locsettings)







#score plotting (Lorenz Curve)---------------

######################################################################-
############# RMSE Score plotting ####################################-
######################################################################-
#data.model.list = c("gauss", "gamma"
#); N = 768; b = 0.2
#dtvec = 0.0005#c(0.005, 0.0005); 
#seed.values = 1988 + (1:5 - 1)*100
#for(i in 1:5){
#  
#  pdf(paste0("b", b, "_RMSE_seed_", seed.values[i], ".pdf"), height = 6, width = 10)
#  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
#  for(j in 1:length(dtvec)){
#    
#    dt = dtvec[j]
#    for(data.model in data.model.list){
#      
#      if(data.model != "gauss"){
#        
#        load(paste0("b", b, "_lorenzrunmarcin_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
#        load(paste0("b", b, "_lorenzrunours_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
#        load(paste0("b", b, "_lorenzrunmarcin_mod_N_", N, "_dt_", dt, ".Rdata"))
#      }else{
#        
#        load(paste0("b", b, "_lorenzrunmarcin_N_", N, "_dt_", dt, ".Rdata"))
#        load(paste0("b", b, "_lorenzrunmarcin_mod_N_", N, "_dt_", dt, ".Rdata"))
#        load(paste0("b", b, "_lorenzrunours_N_", N, "_dt_", dt, ".Rdata"))
#      }
#      RMSE.score = calculateRMSPE(result.lorenz.marcin[[i]]$predsMRA, 
#                                  result.lorenz.ours[[i]]$predsMRA, 
#                                  result.lorenz.marcin[[i]]$XY$x)
#      RMSE.score.2 = calculateRMSPE(result.lorenz.marcin[[i]]$predsMRA, 
#                                    result.lorenz.marcin.2[[i]]$predsMRA, 
#                                    result.lorenz.marcin[[i]]$XY$x)
#      RMSE.score = cbind(RMSE.score, RMSE.score.2[ , 3])
#      if(data.model == "logistic"){
#        temp.pos1 = which(RMSE.score[ , 2] > 10)
#        temp.pos2 = which(RMSE.score[ , 3] > 10)
#        RMSE.score[temp.pos1, 2] = RMSE.score[temp.pos2, 3] = NA
#      }
#      matplot(x = RMSE.score[ , 1], y = RMSE.score[ , -1],
#              col = c("blue", "red", "black"), type = "l", lty = rep(1, 2, 3),
#              xlab = "time", ylab = " ")
#      if(j == 1){title(main = paste0("data.model = ", data.model))}
#      if(data.model == "gauss"){
#        
#        legend("topleft", legend = c("Marcin", "Ours", "Marcin.2"), 
#               col = c("blue", "red", "black"), lty = rep(1, 2))
#        title(ylab = paste0("dt = ", dt))
#      }
#    }
#  }
#  dev.off()
#}

######################################################################-
############# Log Score plotting #####################################-
######################################################################-
#data.model.list = c("gauss", "gamma"
#); N = 768
#dtvec = c(0.005, 0.0005); b = 0.2
#seed.values = 1988 + (1:5 - 1)*100
#ScoreFun = calculateLSs
#for(i in 1:5){
#  
#  pdf(paste0("b", b, "_LogSc_seed_", seed.values[i], ".pdf"), height = 6, width = 10)
#  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(1, 0, 0, 0))
#  for(j in 1:length(dtvec)){
#    
#    dt = dtvec[j]
#    for(data.model in data.model.list){
#      
#      if(data.model != "gauss"){
#        
#        load(paste0("b", b, "_lorenzrunmarcin_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
#        load(paste0("b", b, "_lorenzrunours_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
#        load(paste0("b", b, "_lorenzrunmarcin_mod_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
#      }else{
#        
#        load(paste0("b", b, "_lorenzrunmarcin_N_", N, "_dt_", dt, ".Rdata"))
#        load(paste0("b", b, "_lorenzrunmarcin_mod_N_", N, "_dt_", dt, ".Rdata"))
#        load(paste0("b", b, "_lorenzrunours_N_", N, "_dt_", dt, ".Rdata"))
#      }
#      RMSE.score = ScoreFun(result.lorenz.marcin[[i]]$predsMRA, 
#                            result.lorenz.ours[[i]]$predsMRA, 
#                            result.lorenz.marcin[[i]]$XY$x)
#      # RMSE.score.2 = ScoreFun(result.lorenz.marcin[[i]]$predsMRA, 
#      #                               result.lorenz.marcin.2[[i]]$predsMRA, 
#      #                               result.lorenz.marcin[[i]]$XY$x)
##      # RMSE.score = cbind(RMSE.score, RMSE.score.2[ , 3])
#      duplicate = length(result.lorenz.ours[[i]]$predsMRA)
#      if(duplicate < length(result.lorenz.marcin[[i]]$XY$x)){
#        
#        RMSE.score[duplicate, 3] = NA 
#      }
#      matplot(x = RMSE.score[ , 1], y = RMSE.score[ , -1],
#              col = c("blue", "red"#, "black"
#              ), type = "l", lty = rep(1, 2),
#              xlab = "time", ylab = " ")
#      if(j == 1){title(main = paste0("data.model = ", data.model))}
#      if(data.model == "gauss"){
#        
#        legend("topleft", legend = c("Marcin", "Ours"#, "Marcin.2"
#        ), 
#        col = c("blue", "red"#, "black"
#        ), lty = rep(1, 2))
#        title(ylab = paste0("dt = ", dt))
#      }
#    }
#  }
#  dev.off()
#}

######################################################################-
############# Average Log Score plotting #############################-
######################################################################-
data.model.list = c("gauss", "gamma"
); N = 768
dtvec = 0.0005#c(0.005, 0.0005); 
b = 0.2
seed.values = 1988 + (1:5 - 1)*100
ScoreFun = calculateLSs

for(j in 1:length(dtvec)){
  
  dt = dtvec[j]
  pdf(paste0("filtering-lorenz-chv-ekvl/b", b, "_LogSc_avg_", dt, "_dt", ".pdf"), height = 4, width = 10)
  par(mfrow = c(1, 2), mar = c(4, 3, 1, 1))
  for(data.model in data.model.list){
    
    RMSE.score.avg = array(dim = c(length(seed.values), 40, 3))
    for(i in 1:5){
      
      if(data.model != "gauss"){
        
        load(paste0("filtering-lorenz-chv-ekvl/b", b, "_lorenzrunmarcin_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
        load(paste0("filtering-lorenz-chv-ekvl/b", b, "_lorenzrunours_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
        #load(paste0("b", b, "_lorenzrunmarcin_mod_N_", N, "_dt_", dt, "_", data.model, ".Rdata"))
      }else{
        
        load(paste0("filtering-lorenz-chv-ekvl/b", b, "_lorenzrunmarcin_N_", N, "_dt_", dt, ".Rdata"))
        #load(paste0("b", b, "_lorenzrunmarcin_mod_N_", N, "_dt_", dt, ".Rdata"))
        load(paste0("filtering-lorenz-chv-ekvl/b", b, "_lorenzrunours_N_", N, "_dt_", dt, ".Rdata"))
      }
      RMSE.score.avg[i, , ] = ScoreFun(result.lorenz.marcin[[i]]$predsMRA, 
                                       result.lorenz.ours[[i]]$predsMRA, 
                                       result.lorenz.ours[[i]]$XY$x)
      # RMSE.score.2 = ScoreFun(result.lorenz.marcin[[i]]$predsMRA, 
      #                               result.lorenz.marcin.2[[i]]$predsMRA, 
      #                               result.lorenz.marcin[[i]]$XY$x)
      # RMSE.score = cbind(RMSE.score, RMSE.score.2[ , 3])
      duplicate = length(result.lorenz.ours[[i]]$predsMRA)
      if(duplicate < length(result.lorenz.marcin[[i]]$XY$x)){
        
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
      
      legend("topleft", legend = c("EKVL", "CHV"#, "Marcin.2"
      ), 
      col = c("blue", "red"#, "black"
      ), lty = rep(1, 2))
      #title(ylab = paste0("dt = ", dt))
    }
  }
  dev.off()
}