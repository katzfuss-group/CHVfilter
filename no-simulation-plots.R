library(Rcpp)
Rfast::sourceR(path = "GPvecchia//R//")
source("aux_funs_general.R")

#DAG and matrix sparsity plots (fig 1)-------------------
#install.packages("igraph")
library(igraph); library(ggraph); library(tidyverse); library(ggplot2); library(latex2exp)
library(tidygraph); library(RColorBrewer); library(GPvecchia)
d1 = data.frame(from = "origin", 
                to = paste("group", seq(1,2), sep = ""), edge.cat = rep(1, 2))
d2 = data.frame(from = rep(d1$to, each = 2), 
                to = paste("subgroup", seq(1,4), sep="_"), edge.cat = rep(1, 4))
d3 = data.frame(from = d1$from, to = d2$to, edge.cat = rep(2, 4))
node.col.vals = brewer.pal(7, "Dark2")
graph.data = rbind(d1, d2, d3)
mygraph <- as_tbl_graph(graph.data)

node.label.names <- c("bolditalic(X)^0",
                      "bolditalic(X)[1]",
                      "bolditalic(X)[2]",
                      "bolditalic(X)['1,1']",
                      "bolditalic(X)['1,2']",
                      "bolditalic(X)['2,1']",
                      "bolditalic(X)['2,2']")

graph.plot = ggraph(mygraph, layout = "sugiyama") + theme(
  panel.background = element_rect(fill = "transparent",
                                  colour = NA_character_))
graph.plot = graph.plot +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), start_cap = circle(6, 'mm'),
                 end_cap = circle(8, 'mm'), aes(edge_linetype = factor(edge.cat)),
                 show.legend = FALSE) 
graph.plot = graph.plot +
  geom_node_circle(aes(r = 0.2, fill = factor(1:7)),
                   show.legend = FALSE) + scale_fill_manual(values = node.col.vals)
graph.plot = graph.plot +
  geom_node_text(aes(label = node.label.names), parse = TRUE,show.legend = FALSE, size = 7)
pdf(file = "R_DAG.pdf", width = 8, height = 5)
graph.plot
dev.off()

M = 2; r = rep(1, M + 1); J = rep(2, M)
B = matrix(1, nrow = sum(r * cumprod(c(1, J))), ncol = sum(r))
B[ , 2] = c(0, 2, 3, rep(c(2,3), each = 2))
B[,3] = c(rep(0, 3), 4:7)
temp.mat <- as.matrix(decompressOperator_v4(B, M = M, J = J, r = r, action = "decompress"))
temp.mat[temp.mat == 0] = NA
pdf("sparsity_mat.pdf")
par(mar = c(1, 1, 1, 1))
image(rotateMatrix(as.matrix(temp.mat)), xaxt = "n",
      yaxt = "n", col = node.col.vals)
box()
dev.off()




#plot of points (fig 2 and 7)-------------------
#circular
M = 7; r = c(6, rep(3, M)); J = rep(2, M)
locsettings = generatePoints_v2(M, r, J)
angles = seq(0, 2*pi, length.out = max(locsettings$locsordOrdering) + 1)
angles = angles[ - length(angles)]
angles = angles[locsettings$MatchLocationOrdering]

pdf("locsettings_circ_1D.pdf", width = 8, height = 8)
par(mar = c(0, 1, 0, 1))
plot(x =100, y = 100,
     xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25), 
     type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
#plotrix::draw.circle(0, 0, 1)
cum.points = c(0, cumsum(r * cumprod(c(1, J)))) #no. of points from res. 0 upto res. M
col.select = RColorBrewer::brewer.pal(n = M + 1, "Dark2")
for(m in 0:M)
{
  points.seq = (cum.points[m + 1] + 1):(cum.points[m + 2])
  angles.to.point = angles[points.seq]
  if(m == 0){
    lines.to.draw = c(0, pi)
  }else{
    if(m < M){
      lines.to.draw = zoo::rollapply(angles.to.point, width = r[m + 1], by = r[m + 1], mean)
    }
  }
  points(x = cos(angles.to.point), y = sin(angles.to.point), col = col.select[m+1], pch = 19, cex = 0.6)
  if(m < M){
    for(t in 1:length(lines.to.draw)){
      
      segments(x0 = 0.75 * cos(lines.to.draw[t]), y0 = 0.75 * sin(lines.to.draw[t]),
               x1 = 1.25 * cos(lines.to.draw[t]), y1 = 1.25 * sin(lines.to.draw[t]),
               col = col.select[m+1], lwd = 9 / (m + 1), lty = 1)
    }
  }
}
dev.off()

#linear
M = 2; r = c(1, 1, 1); J = c(1, 1)
locs = seq(0, 1, length.out = 7)
MatchLocationOrdering = c(4, 2, 6, 1, 3, 5, 7)
locsOrderedPoints = locs[MatchLocationOrdering]
locsordOrdering = order(locsOrderedPoints)
col.select = RColorBrewer::brewer.pal(n = length(locs), "Dark2")
pdf("locsettings_lin_1D.pdf", width = 6, height = 2)
par(mfrow = c(1, 3), mar = c(1, 1, 0, 1))
for(i in 0:M)
{
  plot(x = locsOrderedPoints, y = rep(0, length(locsOrderedPoints)), bty = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", lwd = 1, lty = 1, type = "p", 
       xlim = c(0, 1), ylim = c(-0.2, 0.2), cex = 1, pch = 20)
  title(xlab = paste("m =", i), line = 0)
  segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0, col = "black", lwd = 0.5, lty = 1)
  points(x = sort(locsOrderedPoints), y = rep(0, length(locsOrderedPoints)),
         col = "black", cex = 3.5, pch = 20)
  for(m in 0:i){
    
    lines.to.draw = locsOrderedPoints[2^m :(2^(m+1) - 1)]
    if(m < i){
      
      for(t in 1:length(lines.to.draw)){
        segments(x0 = lines.to.draw[t], y0 = -0.2, 
                 x1 = lines.to.draw[t], y1 = 0.2,
                 col = col.select[2^m + t - 1], lwd = 2/(m + 1), lty = 1)
      }
    }
    points(x = locsOrderedPoints[2^m :(2^(m+1) - 1)], y = rep(0, 2^m), 
           cex = 3.5, pch = 20, col = col.select[2^m :(2^(m+1) - 1)])
    #points(x = locsOrderedPoints[2^m :(2^(m+1) - 1)], y = rep(0, 2^m), 
    #       pch = 0, col = "black", cex = 0.8, lwd = 0.5)
  }
}
dev.off()




#defining the locations and priors for Vecchia (1-D)-------------------
locsOrderedPoints = array(dim = 1)
locsOrderedPoints[1] = (0 + 1) / 2
M = 4

for(m in 1:M)
{
  locsOrderedPoints[2^m :(2^(m+1) - 1)] =
    setdiff((seq(0, 1, by = 1/(2^(m+1)))), c(locsOrderedPoints[1:(2^m - 1)], 0, 1))
}
locsOrderedPoints[2^M] = 0; locsOrderedPoints[2^(M+1) - 1] = 1
locs = sort(locsOrderedPoints)
locsordOrdering = order(locsOrderedPoints)
MatchLocationOrdering = match(locsOrderedPoints, locs)#c(0, locs, 1))





M = 4; J = rep(2, M); r = rep(1, M+1); mra.values.list = list(M = M, r = r, J = J)


vecchia_obj = vecchia_specify_modified(matrix(locsOrderedPoints), ordering = 'none', 
                                       conditioning = 'mra', mra.options = mra.values.list) 


PointerDataCompress = extractPointerData(M, J, r, action = "compress")
PointerDataDecompress = extractPointerData(M, J, r, action = "decompress")


covmatrix = fields::Exponential(d = fields::rdist(vecchia_obj$locsord), range = 0.3)                                                                                                                                                                                                  

sig.sel = getMatCov(vecchia_obj, covmatrix)    



inds = Filter(function(i) !is.na(i), as.vector(t(vecchia_obj$U.prep$revNNarray - 1)))

ptrs = c(0, cumsum(apply(vecchia_obj$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))

cov.vals = Filter(function(i) !is.na(i), c(t(sig.sel)))

vals = createUcppM(ptrs, inds, cov.vals)

LMatrix = as.matrix(Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE))



save(mra.values.list, locsOrderedPoints,
     locsordOrdering, MatchLocationOrdering, covmatrix,
     sig.sel, LMatrix, file = "auxoutput_linear.Rdata")



#prior cholesky columns (fig 3)--------------
locsOrderedPoints = array(dim = 1)
locsOrderedPoints[1] = (0 + 1) / 2
M = 4

for(m in 1:M)
{
  locsOrderedPoints[2^m :(2^(m+1) - 1)] =
    setdiff((seq(0, 1, by = 1/(2^(m+1)))), c(locsOrderedPoints[1:(2^m - 1)], 0, 1))
}
locsOrderedPoints[2^M] = 0; locsOrderedPoints[2^(M+1) - 1] = 1
locs = sort(locsOrderedPoints)
locsordOrdering = order(locsOrderedPoints)
MatchLocationOrdering = match(locsOrderedPoints, locs)#c(0, locs, 1))





M = 4; J = rep(2, M); r = rep(1, M+1); mra.values.list = list(M = M, r = r, J = J)


vecchia_obj = vecchia_specify_modified(matrix(locsOrderedPoints), ordering = 'none', 
                                       conditioning = 'mra', mra.options = mra.values.list) 


PointerDataCompress = extractPointerData(M, J, r, action = "compress")
PointerDataDecompress = extractPointerData(M, J, r, action = "decompress")


covmatrix = fields::Exponential(d = fields::rdist(vecchia_obj$locsord), range = 0.3)                                                                                                                                                                                                  

sig.sel = getMatCov(vecchia_obj, covmatrix)    



inds = Filter(function(i) !is.na(i), as.vector(t(vecchia_obj$U.prep$revNNarray - 1)))

ptrs = c(0, cumsum(apply(vecchia_obj$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))

cov.vals = Filter(function(i) !is.na(i), c(t(sig.sel)))

vals = createUcppM(ptrs, inds, cov.vals)

LMatrix = as.matrix(Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE))


pdf("LMatrixPlots.pdf", width = 5, height = 4)
par(mfrow = c(4, 1), mar = c(1, 3, 1, 2), oma = c(3, 1, 0, 1))
cols = viridis::turbo(15)
colIndices.prev = c()
for(m in 0:(M - 1))
{
  if(m == (M - 1)) xaxtype = "s" else xaxtype = "n"
  colIndices = 2^m + (0:(2^(m+1)-2^m -1))
  matplot(x = sort(locsOrderedPoints), y = LMatrix[locsordOrdering ,colIndices], 
          type = "l", lty = rep(1:2, each = length(colIndices)), ylab = "", xlab = "", xaxt = xaxtype,
          col = cols[colIndices])
  
  #title(ylab = paste("resolution (m) =", m), line = 1, outer = TRUE)
  abline(v = locsOrderedPoints[colIndices], lty = 1, col = "gray")
  if(m > 0)
    abline(v = locsOrderedPoints[colIndices.prev], lty = 2, col = "black")
  
  colIndices.prev = append(colIndices.prev, colIndices)
}
mtext("Cholesky factor columns", side = 2, at = 1.5, line = 2, family = "serif")
dev.off()




#plot of normal and compressed matrix (fig 4)-------------
locsOrderedPoints = array(dim = 1)
locsOrderedPoints[1] = (0 + 1) / 2
M = 4

for(m in 1:M)
{
  locsOrderedPoints[2^m :(2^(m+1) - 1)] =
    setdiff((seq(0, 1, by = 1/(2^(m+1)))), c(locsOrderedPoints[1:(2^m - 1)], 0, 1))
}
locsOrderedPoints[2^M] = 0; locsOrderedPoints[2^(M+1) - 1] = 1
locs = sort(locsOrderedPoints)
locsordOrdering = order(locsOrderedPoints)
MatchLocationOrdering = match(locsOrderedPoints, locs)#c(0, locs, 1))





M = 4; J = rep(2, M); r = rep(1, M+1); mra.values.list = list(M = M, r = r, J = J)


vecchia_obj = vecchia_specify_modified(matrix(locsOrderedPoints), ordering = 'none', 
                                       conditioning = 'mra', mra.options = mra.values.list) 


PointerDataCompress = extractPointerData(M, J, r, action = "compress")
PointerDataDecompress = extractPointerData(M, J, r, action = "decompress")


covmatrix = fields::Exponential(d = fields::rdist(vecchia_obj$locsord), range = 0.3)                                                                                                                                                                                                  

sig.sel = getMatCov(vecchia_obj, covmatrix)    



inds = Filter(function(i) !is.na(i), as.vector(t(vecchia_obj$U.prep$revNNarray - 1)))

ptrs = c(0, cumsum(apply(vecchia_obj$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))

cov.vals = Filter(function(i) !is.na(i), c(t(sig.sel)))

vals = createUcppM(ptrs, inds, cov.vals)

LMatrix = as.matrix(Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE))


LMatrixCompressed = decompressOperator_v4(B = LMatrix, M , J, r, action = "compress", 
                                          PointerData = PointerDataCompress) 
pdf("Cholesky_prior_decomp.pdf", width = 6, height = 6)
LMatrix[LMatrix == 0] = NA
par(mar = c(1, 1, 1, 1))
image(rotateMatrix(as.matrix(LMatrix)), 
      xaxt = "n", yaxt = "n", main = "", col = RColorBrewer::brewer.pal(7, "YlGn"))
dev.off()
pdf("Cholesky_prior_comp.pdf", width = 4, height = 6)
LMatrixCompressed[LMatrixCompressed == 0] = NA
par(mar = c(1, 1, 1, 1))
image(rotateMatrix(as.matrix(LMatrixCompressed)), 
      xaxt = "n", yaxt = "n", main = "", col = RColorBrewer::brewer.pal(7, "YlGn"))
dev.off()









#covariance with linear operator(1-D) (fig 5)---------------
locsOrderedPoints = array(dim = 1)
locsOrderedPoints[1] = (0 + 1) / 2
M = 4

for(m in 1:M)
{
  locsOrderedPoints[2^m :(2^(m+1) - 1)] =
    setdiff((seq(0, 1, by = 1/(2^(m+1)))), c(locsOrderedPoints[1:(2^m - 1)], 0, 1))
}
locsOrderedPoints[2^M] = 0; locsOrderedPoints[2^(M+1) - 1] = 1
locs = sort(locsOrderedPoints)
locsordOrdering = order(locsOrderedPoints)
MatchLocationOrdering = match(locsOrderedPoints, locs)#c(0, locs, 1))





M = 4; J = rep(2, M); r = rep(1, M+1); mra.values.list = list(M = M, r = r, J = J)


vecchia_obj = vecchia_specify_modified(matrix(locsOrderedPoints), ordering = 'none', 
                                       conditioning = 'mra', mra.options = mra.values.list) 


PointerDataCompress = extractPointerData(M, J, r, action = "compress")
PointerDataDecompress = extractPointerData(M, J, r, action = "decompress")


covmatrix = fields::Exponential(d = fields::rdist(vecchia_obj$locsord), range = 0.3)                                                                                                                                                                                                  

sig.sel = getMatCov(vecchia_obj, covmatrix)    



inds = Filter(function(i) !is.na(i), as.vector(t(vecchia_obj$U.prep$revNNarray - 1)))

ptrs = c(0, cumsum(apply(vecchia_obj$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))

cov.vals = Filter(function(i) !is.na(i), c(t(sig.sel)))

vals = createUcppM(ptrs, inds, cov.vals)

LMatrix = as.matrix(Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE))


DiffAdvMat1D = (diffAdvOp1DLin(n_G = dim(LMatrix)[1], alphadiffAdv = 0.01, 
                               betadiffAdv = 0.0001))[MatchLocationOrdering, MatchLocationOrdering]

#load("auxoutput_linear.Rdata")
PointerData = list()
PointerData$DataComp = extractPointerData(M = M, J = J, r = r, action = "compress")
PointerData$DataDecomp = extractPointerData(M = M, J = J, r = r, action = "decompress")

EvolMatrix = (tcrossprod(as.matrix(DiffAdvMat1D %*% LMatrix)))[locsordOrdering, locsordOrdering]

EvolMatrixCompDecomp = (tcrossprod(as.matrix(decompressOperator_v4(DiffAdvMat1D %*% decompressOperator_v4(B = LMatrix, M = M, r = r, J = J, 
                                                                                                          action = "compress"),
                                                                   M = M, J = J, r = r, action = "decompress")))
)[locsordOrdering, locsordOrdering]


pdf("covariance_diffadv_operator.pdf", width = 6, height = 6)
DiffAdvMat1D[DiffAdvMat1D == 0] = NA
par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
image(rotateMatrix(as.matrix(DiffAdvMat1D[locsordOrdering, locsordOrdering])),
      xaxt = "n", yaxt = "n", main = "",
      col = RColorBrewer::brewer.pal(7, "YlGn"))
dev.off()
pdf("covariance_diffadv_orig.pdf", width = 6, height = 6)
par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
image(rotateMatrix(EvolMatrix), xaxt = "n", yaxt = "n", main = "",
      col = RColorBrewer::brewer.pal(7, "YlGn"))
dev.off()
pdf("covariance_diffadv_compdecomp.pdf", width = 6, height = 6)
par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
image(rotateMatrix(EvolMatrixCompDecomp), xaxt = "n", yaxt = "n", main = "",
      col = RColorBrewer::brewer.pal(7, "YlGn"))
dev.off()

UKFcovmat = 
  UKFmeanmatredrank.matrix(LMatrix =  LMatrix, 
                          postmean = rep(0, dim(LMatrix)[1]), 
                          evolMat = DiffAdvMat1D, 
                           r = r, RankValue = 5)

UKFcovmat = UKFcovmat$covmat
pdf("covariance_diffadv_ukf.pdf", width = 6, height = 6)
par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
image(rotateMatrix(UKFcovmat[locsordOrdering, locsordOrdering]), xaxt = "n", yaxt = "n", main = "",
      col = RColorBrewer::brewer.pal(7, "YlGn"))
dev.off()




#Basis funs with linear operator (fig 6)---------------
#load("auxoutput_linear.Rdata")
locsOrderedPoints = array(dim = 1)
locsOrderedPoints[1] = (0 + 1) / 2
M = 4

for(m in 1:M)
{
  locsOrderedPoints[2^m :(2^(m+1) - 1)] =
    setdiff((seq(0, 1, by = 1/(2^(m+1)))), c(locsOrderedPoints[1:(2^m - 1)], 0, 1))
}
locsOrderedPoints[2^M] = 0; locsOrderedPoints[2^(M+1) - 1] = 1
locs = sort(locsOrderedPoints)
locsordOrdering = order(locsOrderedPoints)
MatchLocationOrdering = match(locsOrderedPoints, locs)#c(0, locs, 1))





M = 4; J = rep(2, M); r = rep(1, M+1); mra.values.list = list(M = M, r = r, J = J)


vecchia_obj = vecchia_specify_modified(matrix(locsOrderedPoints), ordering = 'none', 
                                       conditioning = 'mra', mra.options = mra.values.list) 


PointerDataCompress = extractPointerData(M, J, r, action = "compress")
PointerDataDecompress = extractPointerData(M, J, r, action = "decompress")


covmatrix = fields::Exponential(d = fields::rdist(vecchia_obj$locsord), range = 0.3)                                                                                                                                                                                                  

sig.sel = getMatCov(vecchia_obj, covmatrix)    



inds = Filter(function(i) !is.na(i), as.vector(t(vecchia_obj$U.prep$revNNarray - 1)))

ptrs = c(0, cumsum(apply(vecchia_obj$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))

cov.vals = Filter(function(i) !is.na(i), c(t(sig.sel)))

vals = createUcppM(ptrs, inds, cov.vals)

LMatrix = as.matrix(Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE))


PointerData = list()
PointerData$DataComp = extractPointerData(M = M, J = J, r = r, action = "compress")
PointerData$DataDecomp = extractPointerData(M = M, J = J, r = r, action = "decompress")
LMatrixCompressed = decompressOperator_v4(B = LMatrix, M , J, r, action = "compress", 
                                          PointerData = PointerData$DataComp) 
DiffAdvMat1D = (diffAdvOp1DLin(n_G = dim(LMatrix)[1], alphadiffAdv = 0.01, 
                               betadiffAdv = 0.0001))[MatchLocationOrdering, MatchLocationOrdering]

EL.orig = as.matrix(DiffAdvMat1D %*% LMatrix)[ , locsordOrdering]

EL.CompDecomp = as.matrix(decompressOperator_v4(DiffAdvMat1D %*% LMatrixCompressed,
                                                M = M, J = J, r = r, action = "decompress"))[ , locsordOrdering]


pdf("basis_diffadv_orig_compdecomp.pdf", width = 5, height = 4)
par(mfrow = c(4, 1), mar = c(1, 3, 1, 2), oma = c(1, 1, 0, 1))
colIndices.prev = c()
for(m in 0:(M - 1))
{
  if(m == (M - 1)) xaxtype = "s" else xaxtype = "n"
  colIndices = 2^m + (0:(2^(m+1)-2^m -1))
  matplot(x = sort(locsOrderedPoints), y = cbind(EL.orig[ , colIndices], EL.CompDecomp[ , colIndices]), 
          type = "l", lty = rep(1:2, each = length(colIndices)), ylab = "n", xlab = "n", xaxt = xaxtype,
          col = rep(c("green", "red"), each = length(colIndices)))
  #title(main = "Plot of L matrix columns", ylab = paste("resolution (m) =", m))
  abline(v = locsOrderedPoints[colIndices], lty = 5, col = "gray")
  
  if(m > 0)
    abline(v = locsOrderedPoints[colIndices.prev], lty = 2, col = "black")
  
  colIndices.prev = append(colIndices.prev, colIndices)
}
dev.off()







