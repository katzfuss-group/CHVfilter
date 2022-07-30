#auxilliary functions needed for experiments-----------------
decompressOperator_v2 = function(A, isZeroOne = TRUE, M = 4, J = 2)
{
  #this code is based on the assumption r = 1, J = 2
  #additional adjustments are needed for further cases
  #this code is also based on assumption that we are using
  #unscented transform on two seperate set of sigma points separately 1:L, (L+1:2L)
  if(isZeroOne){
    
    nrowL = (J^(M+1) + 1)
  }else{
    
    nrowL = (J^(M+1) - 1)
  }
  L = array(0, dim = rep(nrowL, 2))
  
  L[ ,1] = A[ ,1]
  for(m in 1:(M - 1)){
    
    for(j in 1:(J^m)){
      
      #finds out indices from the compressed row
      mTemp = m + 1; J_ct = J; rowindices = J^m + j - 1
      while(mTemp <= M){
        
        newindices = (2^mTemp + (j - 1) * J_ct): (2^mTemp + j * J_ct - 1)
        rowindices = append(rowindices, newindices)
        mTemp = mTemp + 1
        J_ct = J * J_ct
      }
      
      #adjustment for including zero and one, 
      #needs another conditioning if isZeroOne = FALSE
      if(j == 1){
        
        rowindices = append(rowindices, J^(M + 1))
      }else if(j == J^m){
        
        rowindices = append(rowindices, J^(M + 1) + 1)
      }
      L[rowindices, J^m + j - 1] = A[rowindices, (m + 1)]
      #L[rowindices, 2 * J^m + 2 * j - 2] = A[rowindices, 2 * (m + 1) + 1]
    }
  }
  
  #recursion does not work, so seperately writing values for the last resolution
  for(j in 1:J^M){
    
    L[J^M + j - 1, J^M + j - 1] = A[J^M + j - 1 , (M + 1)]
    #L[J^M + j - 1, 2 * J^M + 2 * j - 2] = A[J^M + j - 1 , 2 * (M + 1) + 1]
    if(j == 1){
      
      L[J^(M + 1), J^M + j - 1] = A[J^(M + 1) , M + 1]
    }else if(j == J^M){
      
      L[J^(M + 1) + 1, J^M + j - 1] = A[J^(M + 1) + 1, M + 1]
    }
  }
  
  #inclusion of values when zero and one locations are included 
  if(isZeroOne){
    
    L[J^(M + 1), J^(M + 1)] = A[J^(M + 1), M + 2] 
    L[J^(M + 1) + 1, J^(M + 1) + 1] = A[J^(M + 1) + 1, M + 2] 
  }
  
  return(L)
}

decompressOperator = function(A, M = 4, J = 2)
{
  #this code is based on the assumption r = 1, J = 2
  #additional adjustments are needed for further cases
  #this code is also based on assumption that we are using
  #unscented transform on two seperate set of sigma points separately 1:L, (L+1:2L)
  nrowL = (J^(M+1) - 1)
  L = array(0, dim = c(nrowL, nrowL)) #hardcoded, change for original
  
  L[ ,1] = A[ ,1]
  for(m in 1:(M - 1)){
    
    for(j in 1:(J^m)){
      
      #finds out indices from the compressed row
      mTemp = m + 1; J_ct = J; rowindices = J^m + j - 1
      while(mTemp <= M){
        
        newindices = (2^mTemp + (j - 1) * J_ct): (2^mTemp + j * J_ct - 1)
        rowindices = append(rowindices, newindices)
        mTemp = mTemp + 1
        J_ct = J * J_ct
      }
      
      L[rowindices, J^m + j - 1] = A[rowindices, (m + 1)]
      #L[rowindices, 2 * J^m + 2 * j - 2] = A[rowindices, 2 * (m + 1) + 1]
    }
  }
  
  #recursion does not work, so seperately writing values for the last resolution
  for(j in 1:J^M){
    
    L[J^M + j - 1, J^M + j - 1] = A[J^M + j - 1 , (M + 1)]
    #L[J^M + j - 1, 2 * J^M + 2 * j - 2] = A[J^M + j - 1 , 2 * (M + 1) + 1]
    #if(j == 1){
    
    #  L[J^(M + 1), J^M + j - 1] = A[J^(M + 1) , M + 1]
    #}else if(j == J^M){
    #  
    #  L[J^(M + 1) + 1, J^M + j - 1] = A[J^(M + 1) + 1, M + 1]
    #}
  }
  
  return(L)
}

decompressOperator_v3 = function(B, M = 4, J = rep(2, 4), r = c(2, rep(1, 4)), action = "decompress")
{
  if(action == "compress" && nrow(B) != ncol(B))
  {
    stop("When compressing, the Cholesky factor L must be a square matrix.")
  }
  if(length(J) == 1 || length(r) == 1){
    
    stop("both J and r should be vector valued.")
  }
  if( length(r) != (M + 1)){
    
    stop("r vector should have M + 1 objects, 
         starting from resolution 0 till resolution M")
  }
  
  J = c(1, J); A.cumul.ncols = cumsum(r)
  #this code is based on the assumption r and J are not uniform over the resolutions
  cumul.points = (cumprod(J)) * r
  nrowL = sum(cumul.points)
  complete.row.inds = c()
  complete.col.inds = c()
  complete.vals = c()
  if(action == "decompress")
  {
     dimL = rep(nrowL, 2)
  #  L = array(0, dim = rep(nrowL, 2))
  #  A = B
  }else if(action == "compress"){
     dimL = c(nrowL, sum(r))
  #  A = array(0, dim = c(nrowL, sum(r)))
  #  L = B
  }
  
  for(ctr in 1:(r[1])){
    complete.row.inds = append(complete.row.inds, 1:nrow(B))
    complete.col.inds = append(complete.col.inds, rep(ctr, nrow(B)))
    complete.vals = append(complete.vals, as.vector(B[ ,ctr]))
    #L[ , ctr] = A[ , ctr]
  }
  
  for(m in 1:(M - 1)){
    
    #number of partitions at resolution m
    #noPoints = (r[m + 1]) * prod(J[1:(m + 1)])
    
    noBreaks = prod(J[1:(m + 1)])
    
    cutoff.res.prev = sum(r[1:m] * cumprod(J[1:m]))
    
    for(j in 1: noBreaks){
      
      #finds out indices from the compressed row
      
      #for row m, we start checking for indicators from m+1 onwards
      mTemp = m + 1; 
      
      #the number of division at res. m + 1 is located in J[m + 1 + 1]
      J_ct = J[mTemp + 1]; rowindices = 0#prod(J[1:(m + 1)]) + j - 1
      
      #iterative search for locations
      while(mTemp <= M){
        
        #calculating the cutoff pointer after which the finding for the next resolution starts
        cutoff.res.mTemp = sum((r[1:mTemp]) * cumprod(J[1:mTemp]))
        newindices = (cutoff.res.mTemp + (j - 1) * J_ct * r[mTemp + 1] + 1): 
          (cutoff.res.mTemp + j * J_ct * r[mTemp + 1])
        rowindices = append(rowindices, newindices)
        mTemp = mTemp + 1
        J_ct = J[mTemp + 1] * J_ct ##IMPORTANT: change it to J[mTemp + 1] or something like this
      }
      rowindices = rowindices[-1]
      for(noKnots in 1:r[m + 1]){
        
        KnotIndices  = cutoff.res.prev + 
          (((j - 1) * r[m + 1] + noKnots): ((j * r[m + 1])))
        #"follow this instruction"#rowindices at res m: r[m + 1] * (prod(J[1:(m+1)]) - 1) + noKnots:r[m+1]
        #rowindices = append(KnotIndices, rowindices)
        complete.row.inds = append(complete.row.inds, c(KnotIndices, rowindices))
        if(action == "decompress")
        {
          complete.col.inds = append(complete.col.inds, 
                                     rep(KnotIndices[1],
                                         (length(rowindices) + r[m + 1] - noKnots + 1)))
         # L[rowindices, cutoff.res.prev + (j - 1) * r[m + 1] + noKnots] = 
          complete.vals = append(complete.vals, 
                          as.vector(B[c(KnotIndices, rowindices), (A.cumul.ncols[m] + noKnots)]))  
        }
        else if(action == "compress")
        {
          complete.col.inds = append(complete.col.inds, 
                                     rep(A.cumul.ncols[m] + noKnots, (length(rowindices) + r[m + 1] - noKnots + 1)))
          #A[rowindices, (A.cumul.ncols[m] + noKnots)] = 
          complete.vals = append(complete.vals,
                          as.vector(B[c(KnotIndices, rowindices), cutoff.res.prev + (j - 1) * r[m + 1] + noKnots])) 
        }
      }
      
      #L[rowindices, 2 * J^m + 2 * j - 2] = A[rowindices, 2 * (m + 1) + 1]
    }
  }
  
  #recursion does not work, so seperately writing values for the last resolution
  remove(noBreaks); noBreaks = prod(J)
  remove(cutoff.res.prev); cutoff.res.prev = sum((r[1:M]) * cumprod(J[1:M]))
  for(j in 1:noBreaks){
    
    for (noKnots in 1:r[M + 1]) {
      
      remove(KnotIndices)
      KnotIndices = cutoff.res.prev + 
        (((j - 1) * r[M + 1] + noKnots): ((j * r[M + 1])))
      
      complete.row.inds = append(complete.row.inds, KnotIndices)
      if(action == "decompress")
      {
        complete.col.inds = append(complete.col.inds, 
                                   rep(KnotIndices[1], r[m + 1] - noKnots + 1))
        #L[KnotIndices, KnotIndices] = 
         complete.vals = append(complete.vals, 
                                B[KnotIndices , A.cumul.ncols[M] + noKnots]) 
      }
      else if(action == "compress")
      {
        complete.col.inds = append(complete.col.inds, 
                                   rep(A.cumul.ncols[M] + noKnots, length(KnotIndices)))
        #A[KnotIndices , A.cumul.ncols[M] + noKnots] =
        complete.vals = append(complete.vals, B[KnotIndices, KnotIndices[1]])  
      }
    }
  }
  #complete.row.inds = complete.row.inds[-1]
  #complete.col.inds = complete.col.inds[-1]
  #complete.vals = complete.vals[-1]
  L = Matrix::sparseMatrix(i = complete.row.inds, j = complete.col.inds,
                           x = complete.vals, dims = dimL)
  return(L)
}


extractPointerData = function(M, J, r, action = "decompress")
{
  if(length(J) == 1 || length(r) == 1){
    
    stop("both J and r should be vector valued.")
  }
  if( length(r) != (M + 1)){
    
    stop("r vector should have M + 1 objects, 
         starting from resolution 0 till resolution M")
  }
  
  J = c(1, J); A.cumul.ncols = cumsum(r)
  cumul.points = (cumprod(J)) * r
  nrowL = sum(cumul.points)
  complete.row.inds = c()
  complete.col.inds = c()
  complete.vals.cols = c()
  if(action == "decompress")
  {
    dimL = rep(nrowL, 2)
  }else if(action == "compress"){
    dimL = c(nrowL, sum(r))
    
  }
  
  for(ctr in 1:(r[1])){
    complete.row.inds = append(complete.row.inds, ctr:nrowL)
    complete.col.inds = append(complete.col.inds, rep(ctr, (nrowL - ctr + 1)))
    complete.vals.cols = append(complete.vals.cols, rep(ctr, (nrowL - ctr + 1)))
  }
  
  for(m in 1:(M - 1)){
    
    #number of partitions at resolution m
    
    noBreaks = prod(J[1:(m + 1)])
    
    cutoff.res.prev = sum(r[1:m] * cumprod(J[1:m]))
    
    for(j in 1: noBreaks){
      
      #finds out indices from the compressed row
      
      #for row m, we start checking for indicators from m+1 onwards
      mTemp = m + 1; 
      
      #the number of division at res. m + 1 is located in J[m + 1 + 1]
      J_ct = J[mTemp + 1]; rowindices = 0#prod(J[1:(m + 1)]) + j - 1
      
      #iterative search for locations
      while(mTemp <= M){
        
        #calculating the cutoff pointer after which the finding for the next resolution starts
        cutoff.res.mTemp = sum((r[1:mTemp]) * cumprod(J[1:mTemp]))
        newindices = (cutoff.res.mTemp + (j - 1) * J_ct * r[mTemp + 1] + 1): 
          (cutoff.res.mTemp + j * J_ct * r[mTemp + 1])
        rowindices = append(rowindices, newindices)
        mTemp = mTemp + 1
        J_ct = J[mTemp + 1] * J_ct ##IMPORTANT: change it to J[mTemp + 1] or something like this
      }
      rowindices = rowindices[-1]
      for(noKnots in 1:r[m + 1]){
        
        KnotIndices  = cutoff.res.prev + 
          (((j - 1) * r[m + 1] + noKnots): ((j * r[m + 1])))
        complete.row.inds = append(complete.row.inds, c(KnotIndices, rowindices))
        if(action == "decompress")
        {
          complete.col.inds = append(complete.col.inds, 
                                     rep(KnotIndices[1],
                                         (length(rowindices) + r[m + 1] - noKnots + 1)))
          #complete.vals = append(complete.vals, 
          #                       as.vector(B[c(KnotIndices, rowindices), (A.cumul.ncols[m] + noKnots)])) 
          complete.vals.cols = append(complete.vals.cols, 
                                      rep((A.cumul.ncols[m] + noKnots), (length(rowindices) + r[m + 1] - noKnots + 1)))
        }
        else if(action == "compress")
        {
          complete.col.inds = append(complete.col.inds, 
                                     rep(A.cumul.ncols[m] + noKnots, (length(rowindices) + r[m + 1] - noKnots + 1)))
          #complete.vals = append(complete.vals,
          #                       as.vector(B[c(KnotIndices, rowindices), cutoff.res.prev + (j - 1) * r[m + 1] + noKnots]))
          complete.vals.cols = append(complete.vals.cols,
                                      rep(cutoff.res.prev + (j - 1) * r[m + 1] + noKnots, (length(rowindices) + r[m + 1] - noKnots + 1)))
        }
      }
    }
  }
  
  #recursion does not work, so seperately writing values for the last resolution
  remove(noBreaks); noBreaks = prod(J)
  remove(cutoff.res.prev); cutoff.res.prev = sum((r[1:M]) * cumprod(J[1:M]))
  for(j in 1:noBreaks){
    
    for (noKnots in 1:r[M + 1]) {
      
      remove(KnotIndices)
      KnotIndices = cutoff.res.prev + 
        (((j - 1) * r[M + 1] + noKnots): ((j * r[M + 1])))
      
      complete.row.inds = append(complete.row.inds, KnotIndices)
      if(action == "decompress")
      {
        complete.col.inds = append(complete.col.inds, 
                                   rep(KnotIndices[1], r[m + 1] - noKnots + 1))
        #complete.vals = append(complete.vals, 
        #                       B[KnotIndices , A.cumul.ncols[M] + noKnots]) 
        complete.vals.cols = append(complete.vals.cols,
                                    rep(A.cumul.ncols[M] + noKnots, r[m + 1] - noKnots + 1))
      }
      else if(action == "compress")
      {
        complete.col.inds = append(complete.col.inds, 
                                   rep(A.cumul.ncols[M] + noKnots, length(KnotIndices)))
        #complete.vals = append(complete.vals, B[KnotIndices, KnotIndices[1]]) 
        complete.vals.cols = append(complete.vals.cols, 
                                    rep(KnotIndices[1], length(KnotIndices)))
      }
    }
  }
  DataList = list(complete.row.inds = complete.row.inds,
                  complete.col.inds = complete.col.inds,
                  complete.vals.cols = complete.vals.cols)
  return(DataList)
}



decompressOperator_v4 = function(B, M, J, r, action = "decompress", PointerData = NULL){
  
  dimB = dim(B)
  if(action == "compress"){
    
    if(dimB[1] != dimB[2]){
      
      stop("Make sure that the input matrix is a square matrix.")
    }
    dimL = c((dim(B))[1], sum(r))
  }else if(action == "decompress"){
    
    if(dimB[2] != sum(r))
    {
      stop("The input matrix must have number of columns equal to the total number of knots.")
    }
    dimL = rep(dimB[1], 2)
  }
  if(is.null(PointerData))
  {
    PointerData = extractPointerData(M, J, r, action = action) 
  }
  nvalues = length(PointerData$complete.row.inds)
  valuesvector = rep(NA, nvalues)
  for(i in 1:nvalues){
    
    valuesvector[i] = B[(PointerData$complete.row.inds)[i], 
                        (PointerData$complete.vals.cols)[i]]
  }
  L = Matrix::sparseMatrix(i = PointerData$complete.row.inds,
                           j = PointerData$complete.col.inds,
                           x = valuesvector)
  return(L)
}



filterFmat = function(Fmat = NULL, PointerData = NULL){
  
  if(is.null(Fmat) || is.null(PointerData)){
    
    stop("Both of Fmat and PointerData variables have to be input properly.")
  }
  nvalues = length(PointerData$complete.row.inds)
  valuesvector = rep(NA, nvalues)
  for(i in 1:nvalues){
    
    valuesvector[i] = Fmat[(PointerData$complete.row.inds)[i], 
                        (PointerData$complete.vals.cols)[i]]
  }
  for(i in 1:nvalues){
    
    L = Matrix::sparseMatrix(i = PointerData$complete.row.inds,
                             j = PointerData$complete.vals.cols,
                             x = valuesvector)
  }
  return(L)
}





UKFmeanmatredrank = function(LMatrix, postmean, evolFunction, r, RankValue){
  
  LMatrix = as.matrix(LMatrix)
  if(missing(RankValue)){
    RankValue = floor(0.1 * (dim(LMatrix)[1]))
  }
  meanVector = as.vector(postmean)
  
  
  #re-adjusting mean and variance weights
  LUKFRedRank = RankValue; kappaUKF = 0; alphaUKF = 0.9
  gammaUKFRedRank = (alphaUKF)^2 * (LUKFRedRank + kappaUKF) - LUKFRedRank
  coefUKFRedRank = LUKFRedRank + gammaUKFRedRank #used as denominators in various places of calculating UKF
  betaUKF = 2 #general for Gaussian distribution
  meanUKFweightsRedRank = varUKFweightsRedRank = rep(1/ (2 * coefUKFRedRank), 2 * LUKFRedRank + 1)
  varUKFweightsRedRank[1] = (gammaUKFRedRank / coefUKFRedRank) + (1 - alphaUKF^2 + betaUKF)
  meanUKFweightsRedRank[1] = (gammaUKFRedRank / coefUKFRedRank)
  
  SigmaPointMatrices = matrix(0, nrow = dim(LMatrix)[1], ncol = 2 * RankValue + 1)
  SigmaPointMatrices[ , 1] = meanVector
  for(i in 1:RankValue)
  {
    #the first column entry is eean vector; so ith col of UKFmatrix is (i+1) col of Lcompressed
    SigmaPointMatrices[ , 1 + i] = meanVector + (sqrt(coefUKFRedRank)) * LMatrix[ , i]
    SigmaPointMatrices[ , RankValue + 1 + i] = meanVector - (sqrt(coefUKFRedRank)) * LMatrix[ , i]
  }
  SigmaPointMatrices = apply(SigmaPointMatrices, 2, evolFunction)
  
  
  
  #calculating the UKF matrix
  UKFmeanRedRank = colSums(tcrossprod(meanUKFweightsRedRank, SigmaPointMatrices))
  UKFmatrixRedRank = array(0, dim = dim(LMatrix))
  for(i in 1:(2 * RankValue + 1))
  {
    UKFmatrixRedRank = UKFmatrixRedRank + varUKFweightsRedRank[i] * tcrossprod(SigmaPointMatrices[ ,i] - UKFmeanRedRank)
  }
  values.list = list()
  values.list$mean = UKFmeanRedRank
  values.list$covmat = UKFmatrixRedRank
  return(values.list)
}



UKFmeanmatredrank.matrix = function(LMatrix, postmean, evolMat, r, RankValue){
  
  LMatrix = as.matrix(LMatrix)
  if(missing(RankValue)){
    RankValue = floor(0.1 * (dim(LMatrix)[1]))
  }
  meanVector = as.vector(postmean)
  
  
  #re-adjusting mean and variance weights
  LUKFRedRank = RankValue; kappaUKF = 0; alphaUKF = 0.9
  gammaUKFRedRank = (alphaUKF)^2 * (LUKFRedRank + kappaUKF) - LUKFRedRank
  coefUKFRedRank = LUKFRedRank + gammaUKFRedRank #used as denominators in various places of calculating UKF
  betaUKF = 2 #general for Gaussian distribution
  meanUKFweightsRedRank = varUKFweightsRedRank = rep(1/ (2 * coefUKFRedRank), 2 * LUKFRedRank + 1)
  varUKFweightsRedRank[1] = (gammaUKFRedRank / coefUKFRedRank) + (1 - alphaUKF^2 + betaUKF)
  meanUKFweightsRedRank[1] = (gammaUKFRedRank / coefUKFRedRank)
  
  SigmaPointMatrices = matrix(0, nrow = dim(LMatrix)[1], ncol = 2 * RankValue + 1)
  SigmaPointMatrices[ , 1] = meanVector
  for(i in 1:RankValue)
  {
    #the first column entry is eean vector; so ith col of UKFmatrix is (i+1) col of Lcompressed
    SigmaPointMatrices[ , 1 + i] = meanVector + (sqrt(coefUKFRedRank)) * LMatrix[ , i]
    SigmaPointMatrices[ , RankValue + 1 + i] = meanVector - (sqrt(coefUKFRedRank)) * LMatrix[ , i]
  }
  evolMat = as.matrix(evolMat); evolMat[is.na(evolMat)] = 0
  SigmaPointMatrices = as.matrix(evolMat) %*% SigmaPointMatrices
  
  
  
  #calculating the UKF matrix
  UKFmeanRedRank = colSums(tcrossprod(meanUKFweightsRedRank, SigmaPointMatrices))
  UKFmatrixRedRank = array(0, dim = dim(LMatrix))
  for(i in 1:(2 * RankValue + 1))
  {
    UKFmatrixRedRank = UKFmatrixRedRank + varUKFweightsRedRank[i] * tcrossprod(SigmaPointMatrices[ ,i] - UKFmeanRedRank)
  }
  values.list = list()
  values.list$mean = UKFmeanRedRank
  values.list$covmat = UKFmatrixRedRank
  return(values.list)
}



diffAdvOp1DLin = function(alphadiffAdv = NA, betadiffAdv = NA, n_G){
  
  if(is.na(alphadiffAdv)){
    alphadiffAdv = 0.5/n_G
  }
  if(is.na(betadiffAdv)){
    betadiffAdv = 0.35/(n_G ^2)
  }
  
  c1 = 1 - 2 * betadiffAdv * (n_G ^2)
  c2 = 0.5 * alphadiffAdv * n_G + betadiffAdv * (n_G ^2)
  c3 = -0.5 * alphadiffAdv * n_G + betadiffAdv * (n_G ^2)
  
  A = matrix(0, n_G, n_G)
  
  row.inds = col.inds = vals = c()
  for(i in 2:(n_G - 1))
  {
    row.inds = append(row.inds, rep(i, 3))
    col.inds = append(col.inds, c(i - 1, i, i + 1))
    vals = append(vals, c(c3, c1, c2))
    #A[i, c(i - 1, i, i + 1)] = c(c3, c1, c2)
  }
  row.inds = append(row.inds, rep(1, 2)); col.inds = append(col.inds, c(1, 2))
  vals = append(vals, c(c1, c2))
  #A[1, c(1, 2)] = c(c1, c2)
  row.inds = append(row.inds, rep(n_G, 2)); col.inds = append(col.inds, c(n_G - 1, n_G))
  vals = append(vals, c(c3, c1))
  #A[n_G, c(n_G - 1, n_G)] = c(c3, c1)
  A = Matrix::sparseMatrix(i = row.inds, j = col.inds, x = vals)
  return(A)
}


evolFun.linear = function(alpha = NA, beta = NA, X){
  
  Et = diffAdvOp1DLin(n_G = nrow(X), alpha = alpha, beta = beta)
  newX = matrix(Et %*% as.matrix(X), ncol = 1)
  return(newX)
}

diagOpt = function(x){
  
  n_G = nrow(x)
  returnMat = matrix(diag(rep(2, n_G)) %*% x, ncol = 1)
  return(returnMat)
}

diffAdvOp1DCirc = function(alpha = NA, beta = NA, n_G){
  
  if(is.na(alpha)){
    alphadiffAdv = 0.5/n_G
  }
  if(is.na(beta)){
    betadiffAdv = 0.35/(n_G ^2)
  }
  
  c1 = 1 - 2 * betadiffAdv * (n_G ^2)
  c2 = 0.5 * alphadiffAdv * n_G + betadiffAdv * (n_G ^2)
  c3 = -0.5 * alphadiffAdv * n_G + betadiffAdv * (n_G ^2)
  
  A = matrix(0, n_G, n_G)
  
  for(i in 2:(n_G - 1))
  {
    A[i, c(i - 1, i, i + 1)] = c(c3, c1, c2)
  }
  A[1, c(1, 2, n_G)] = c(c1, c2, c3)
  A[n_G, c(1, n_G - 1, n_G)] = c(c2, c3, c1)
  return(A)
}


rotateMatrix <- function(x) t(apply(x, 2, rev))


  
  
generatePoints = function(M, r, J, domain = "circular"){
  
  if(domain != "circular")
  {
    stop("This code has not been developed for points outside circular domain.")
  }
  if(length(J) != (M) || length(r) != (M + 1))
  {
    stop("Check length of both r and J to match compatibility with M.")
  }
  if(unique(J[2:M] - r[2:M]) != 1)
  {
    stop("This code has been written for a special configuration of r and J.")
  }
  noBreaks = cumprod(c(1, J))
  noPoints = r * noBreaks; cumulnoPoints = cumsum(noPoints)
  
  locsOrderedPoints = array(dim = 1)
  locsOrderedPoints[1: (r[1])] = (seq(0, 1, length.out = r[1] + 1))[ - (r[1] + 1)]
  for(m in 1:M)
  {
    prevlocs = c(sort(locsOrderedPoints), 1)
    newlocs = c()
    for(k in 1:noBreaks[m + 1])
    {
      #this line is written with specific configuration of r and J in mind.
      newlocs = append(newlocs, (seq(prevlocs[k], prevlocs[k + 1], 
                                     length.out = (r[m + 1] + 2)))[-c(1, (r[m+1] + 2))])
      #tempindex = k - cumulnoPoints[m]
      #locsOrderedPoints[k] = 0(prevlocs[tempindex] + prevlocs[tempindex + 1])
    }
    locsOrderedPoints = append(locsOrderedPoints, newlocs)
    
    #locsOrderedPoints[2^m :(2^(m+1) - 1)] =
    #  setdiff((seq(0, 1, by = 1/(2^(m+1)))), c(locsOrderedPoints[1:(2^m - 1)], 0, 1))
  }
  locs = sort(locsOrderedPoints)
  locsordOrdering = order(locsOrderedPoints)
  MatchLocationOrdering = match(locsOrderedPoints, locs)
  data.location = list(locsOrderedPoints = locsOrderedPoints, 
                       locs = locs, 
                       locsordOrdering = locsordOrdering, 
                       MatchLocationOrdering = MatchLocationOrdering)
  return(data.location)
}


generatePoints_v2 = function(M, r, J, domain = "circular", flag = 0){
  
  if(domain != "circular")
  {
    stop("This code has not been developed for points outside circular domain.")
  }
  if(length(J) != (M) || length(r) != (M + 1))
  {
    stop("Check length of both r and J to match compatibility with M.")
  }
  #if(unique(J[2:M] - r[2:M]) != 1)
  #{
  #  stop("This code has been written for a special configuration of r and J.")
  #}
  if(length(unique(r[-1])) != 1 || (r[1] != 2 * r[2]))
  {
    stop("Function needs to be modified for non-uniform r.")
  }
  noBreaks = cumprod(c(1, J))
  noPoints = r * noBreaks; cumulnoPoints = cumsum(noPoints)
  
  totalpoints = sum(noPoints); pointdiff = 1/ (totalpoints)
  
  midpoints = list()
  midpoints[[1]] = c(0, 0.5)
  for(i in 2:(M + 1)){midpoints[[i]] = zoo::rollmean(c(sort(unlist(midpoints)), 1), 2)}
  
  midpointsvec = unlist(midpoints); midpointsveclength = length(midpointsvec)
  
  locsOrderedPoints = c()
  for(i in 1:midpointsveclength)
  {
    locsOrderedPoints = append(locsOrderedPoints,
                               seq(from = midpointsvec[i] - ((r[2] - 1) * pointdiff/2),
                                   to = midpointsvec[i] + ((r[2] - 1) * pointdiff/2),
                                   length.out = r[2]))
    if(i == 1)
    {
      temp.location = which(locsOrderedPoints < 0)
      locsOrderedPoints[temp.location] = 1 + locsOrderedPoints[temp.location]
      locsOrderedPoints = sort(locsOrderedPoints)
    }
  }
  #locsOrderedPoints[1: (r[1])] = (seq(0, 1, length.out = r[1] + 1))[ - (r[1] + 1)]
  #for(m in 1:M)
  #{
  #  prevlocs = c(sort(locsOrderedPoints), 1)
  #  newlocs = c()
  #  for(k in 1:noBreaks[m + 1])
  #  {
  #    #this line is written with specific configuration of r and J in mind.
  #    newlocs = append(newlocs, (seq(prevlocs[k], prevlocs[k + 1], 
  #                                   length.out = (r[m + 1] + 2)))[-c(1, (r[m+1] + 2))])
  #    #tempindex = k - cumulnoPoints[m]
  #    #locsOrderedPoints[k] = 0(prevlocs[tempindex] + prevlocs[tempindex + 1])
  #  }
  #  locsOrderedPoints = append(locsOrderedPoints, newlocs)
    
  #  #locsOrderedPoints[2^m :(2^(m+1) - 1)] =
  #  #  setdiff((seq(0, 1, by = 1/(2^(m+1)))), c(locsOrderedPoints[1:(2^m - 1)], 0, 1))
  #}
  locs = sort(locsOrderedPoints)
  locsordOrdering = order(locsOrderedPoints)
  MatchLocationOrdering = match(locsOrderedPoints, locs)
  if(flag == 0){
    
    data.location = list(locsOrderedPoints = locsOrderedPoints, 
                         locs = locs, 
                         locsordOrdering = locsordOrdering, 
                         MatchLocationOrdering = MatchLocationOrdering)
  }else if(flag == 1){
    
    data.location = list(locsOrderedPoints = locsOrderedPoints, 
                         locs = locs, 
                         locsordOrdering = locsordOrdering, 
                         MatchLocationOrdering = MatchLocationOrdering,
                         midPointsList = midpoints)
  }
  
  return(data.location)
}


rotateLocations = function(LocationData = NULL, M = NULL, r = NULL, J = 
                             NULL, domain = "circular", times = 1){
  
  if(is.null(LocationData))
  {
    LocationData = generatePoints(M, r, J, domain = domain)
  }
  NewMatchLocationOrdering = LocationData$MatchLocationOrdering + times
  if(times > 0){
    
    temp_positions = which(NewMatchLocationOrdering > max(LocationData$MatchLocationOrdering))
    NewMatchLocationOrdering[temp_positions] = 
      NewMatchLocationOrdering[temp_positions] - max(LocationData$MatchLocationOrdering)
  }else{
    
    temp_positions = which(NewMatchLocationOrdering < min(LocationData$MatchLocationOrdering))
    NewMatchLocationOrdering[temp_positions] = 
      MatchLocationOrdering[temp_positions] + max(LocationData$MatchLocationOrdering)
  }
  rotatedLocs = (LocationData$locs)[NewMatchLocationOrdering] 
  #decimalpoints = decimalplaces(1/max(LocationData$locsordOrdering))
  #locdiff = round(LocationData$locs[2] - LocationData$locs[1], decimalpoints + 1)
  #rotatedLocs = round(LocationData$locsOrderedPoints  + 
  #  times * locdiff,
  #  decimalpoints + 1)
  #temp_positions = which(rotatedLocs >= 1)
  #rotatedLocs[temp_positions] = rotatedLocs[temp_positions] - 1
  rotatedOrdering = match(rotatedLocs, LocationData$locsOrderedPoints)
  #NewlocsordOrdering = order(rotatedLocs)
  MatchLocationOrdering = match(rotatedLocs, LocationData$locs)
  data.location = list()
  data.location$locsOrderedPoints = rotatedLocs
  data.location$locsordOrdering = order(rotatedLocs)
  data.location$locs = sort(rotatedLocs)
  data.location$rotatedOrdering = rotatedOrdering
  data.location$MatchLocationOrdering = NewMatchLocationOrdering
  return(data.location)
}



evolFun = function(X) (b*Lorenz04M2Sim(as.numeric(X)/b, Force, K, dt, M, iter = 1, burn = 0))
evolFunMu = function(X) (b*Lorenz04M2Sim(as.numeric(X)/b, Force, K, dt, M, iter = 1, burn = 0))

filter_v2 = function(approx.name, XY, mod, PointerData, DataSetup.init, mu.tt = mu, Sig0 = Sig0){
  
  approx = approximations[[approx.name]]
  PointerDataComp = PointerData$DataComp
  PointerDataDecomp = PointerData$DataDeomp
  locsordOrdering = DataSetup.init$locsordOrdering
  MatchLocationOrdering = DataSetup.init$MatchLocationOrdering
  
  Tmax = length(XY$x)
  preds = list()
  obs.aux = as.numeric(XY$y[[1]])
  
  cat(paste("\tfiltering: t=1\n"))
  
  cat("\t\tcalculating forecast moments\n")
  covmodel = GPvecchia::getMatCov(approx, as.matrix(Sig0))
  #mu.tt = 
  #mu.tt = mu
  
  cat("\t\tcalculating posterior\n")
  inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
  vals = createUcppM(ptrs, inds, cov.vals)
  L.tt =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  
  
  #preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt,
  #                                               likelihood_model = data.model, covmodel = covmodel,
  #                                               covparms = NULL, likparms = lik.params, return_all = TRUE, mod = mod)
  
  #cat("\t\tsaving the moments\n")
  #L.tt  = getLtt(approx, preds.aux)#getL00(approx, covmodel, locs) #####plesase check this
  #mu.tt = matrix(preds.aux$mean, ncol = 1)
  #preds[[1]] = list(state = mu.tt, W = preds.aux$W)#, V = preds.aux$V)
  
  if (Tmax == 1) { 
    return( preds )
  } 
  
  for (t in 1:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    if(mod == "standard")
    {
      cat("\t\tcalculating gradient...\n")
      Et = (Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering]
      Fmat = Et %*% L.tt
      covmodel = GPvecchia::getMatCov(approx, as.matrix(Fmat %*% Matrix::t(Fmat) + Sigt))
    }else if(mod == "noevolerror"){
      
      storetempmat = (apply(decompressOperator_v4(L.tt, action = "compress", PointerData = PointerDataComp, M = M.multires, J = J, r = r)[locsordOrdering, ], 2, evolFun))[MatchLocationOrdering, ]
      L.tt.aux = decompressOperator_v4(storetempmat, action = "decompress", PointerData = PointerDataDecomp, M = M.multires, J = J, r = r)
      covmodel = L.tt.aux
    }else{
      
      stop("Check the `mod' variable.")
    }
    
    
    cat("\t\tcalculating forecast moments\n")
    forecast = (evolFunMu(mu.tt[locsordOrdering]))[MatchLocationOrdering]
    
    
    cat("\t\tcalculating posterior\n")
    preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                        likelihood_model = data.model, covmodel = covmodel,
                                        covparms = covparms, likparms = lik.params, return_all = TRUE, mod = mod)
    cat("\t\tsaving the moments\n")
    
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
  }
  return( preds )
}


filter_v22 = function(approx.name, XY, mod, PointerData, DataSetup.init, mu.tt = mu, Sig0 = Sig0){
  
  approx = approximations[[approx.name]]
  PointerDataComp = PointerData$DataComp
  PointerDataDecomp = PointerData$DataDeomp
  locsordOrdering = DataSetup.init$locsordOrdering
  MatchLocationOrdering = DataSetup.init$MatchLocationOrdering
  
  Tmax = length(XY$x)
  preds = list()
  obs.aux = as.numeric(XY$y[[1]])
  
  cat(paste("\tfiltering: t=1\n"))
  
  cat("\t\tcalculating forecast moments\n")
  covmodel = getMatCov(approx, as.matrix(Sig0))
  #mu.tt = 
  #mu.tt = mu
  
  #cat("\t\tcalculating posterior\n")
  inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
  vals = createUcppM(ptrs, inds, cov.vals)
  L.tt =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  
  
  #preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt,
  #                                               likelihood_model = data.model, covmodel = covmodel,
  #                                               covparms = NULL, likparms = lik.params, return_all = TRUE, mod = mod)
  
  #cat("\t\tsaving the moments\n")
  #L.tt  = getLtt(approx, preds.aux)#getL00(approx, covmodel, locs) #####plesase check this
  #mu.tt = matrix(preds.aux$mean, ncol = 1)
  #preds[[1]] = list(state = mu.tt, W = preds.aux$W)#, V = preds.aux$V)
  
  if (Tmax == 1) { 
    return( preds )
  } 
  
  for (t in 1:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    if(mod == "standard")
    {
      cat("\t\tcalculating gradient...\n")
      Et = tryCatch((Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering],
                    error = function(e) 1)
      if(!is.matrix(Et)){
        
        break
      }
      Fmat = Et %*% L.tt
      covmodel = GPvecchia::getMatCov(approx, as.matrix(Fmat %*% Matrix::t(Fmat)))
    }else if(mod == "noevolerror"){
      
      if(t == 1){
        
        cat("\t\tcalculating gradient...\n")
        Et = (Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering]
        Fmat = Et %*% L.tt
        covmodel = GPvecchia::getMatCov(approx, as.matrix(Fmat %*% Matrix::t(Fmat) + Sigt))
        inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
        ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
        cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
        vals = createUcppM(ptrs, inds, cov.vals)
        covmodel =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
      }else{
        
        storetempmat = (apply(decompressOperator_v4(L.tt, action = "compress", PointerData = PointerDataComp, M = M.multires, J = J, r = r)[locsordOrdering, ], 2, evolFun))[MatchLocationOrdering, ]
        L.tt.aux = decompressOperator_v4(storetempmat, action = "decompress", PointerData = PointerDataDecomp, M = M.multires, J = J, r = r)
        covmodel = L.tt.aux
      }
    }else{
      
      stop("Check the `mod' variable.")
    }
    
    
    cat("\t\tcalculating forecast moments\n")
    forecast = (evolFunMu(mu.tt[locsordOrdering]))[MatchLocationOrdering]
    
    
    cat("\t\tcalculating posterior\n")
    preds.aux = tryCatch(calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                        likelihood_model = data.model, covmodel = covmodel,
                                        covparms = covparms, likparms = lik.params, return_all = TRUE, mod = mod),
                         error = function(e) 1)
    if((!is.list(preds.aux)) || any(is.na(preds.aux$mean)))
    {
      break
    }
    cat("\t\tsaving the moments\n")
    
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
  }
  return( preds )
}


filter_v23 = function(approx.name, XY, mod, PointerData, DataSetup.init, mu.tt = mu, Sig0 = Sig0){
  
  approx = approximations[[approx.name]]
  PointerDataComp = PointerData$DataComp
  PointerDataDecomp = PointerData$DataDeomp
  locsordOrdering = DataSetup.init$locsordOrdering
  MatchLocationOrdering = DataSetup.init$MatchLocationOrdering
  
  Tmax = length(XY$x)
  preds = list()
  obs.aux = as.numeric(XY$y[[1]])
  
  cat(paste("\tfiltering: t=1\n"))
  
  cat("\t\tcalculating forecast moments\n")
  covmodel = getMatCov(approx, as.matrix(Sig0))
  #mu.tt = 
  #mu.tt = mu
  
  cat("\t\tcalculating posterior\n")
  inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
  vals = createUcppM(ptrs, inds, cov.vals)
  L.tt =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  
  
  #preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt,
  #                                               likelihood_model = data.model, covmodel = covmodel,
  #                                               covparms = NULL, likparms = lik.params, return_all = TRUE, mod = mod)
  
  #cat("\t\tsaving the moments\n")
  #L.tt  = getLtt(approx, preds.aux)#getL00(approx, covmodel, locs) #####plesase check this
  #mu.tt = matrix(preds.aux$mean, ncol = 1)
  #preds[[1]] = list(state = mu.tt, W = preds.aux$W)#, V = preds.aux$V)
  
  if (Tmax == 1) { 
    return( preds )
  } 
  
  for (t in 1:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    if(mod == "standard")
    {
      cat("\t\tcalculating gradient...\n")
      Et = tryCatch((Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering],
                    error = function(e) 1)
      if(is.numeric(Et)){
        
        break
      }
      
      Fmat = Et %*% L.tt
      Fmat.mod = filterFmat(Fmat = Fmat, PointerData = PointerDataComp)
      covmodel = getMatCov(approx, as.matrix(Fmat.mod %*% Matrix::t(Fmat.mod)))
    }else if(mod == "noevolerror"){
      
      storetempmat = (apply(decompressOperator_v4(L.tt, action = "compress", PointerData = PointerDataComp, M = M.multires, J = J, r = r)[locsordOrdering, ], 2, evolFun))[MatchLocationOrdering, ]
      L.tt.aux = decompressOperator_v4(storetempmat, action = "decompress", PointerData = PointerDataDecomp, M = M.multires, J = J, r = r)
      covmodel = L.tt.aux
    }else{
      
      stop("Check the `mod' variable.")
    }
    
    
    cat("\t\tcalculating forecast moments\n")
    forecast = (evolFunMu(mu.tt[locsordOrdering]))[MatchLocationOrdering]
    
    
    cat("\t\tcalculating posterior\n")
    preds.aux = tryCatch(calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                        likelihood_model = data.model, covmodel = covmodel,
                                        covparms = covparms, likparms = lik.params, return_all = TRUE, mod = mod),
                         error = function(e) 1)
    if((!is.list(preds.aux)) || any(is.na(preds.aux$mean)))
    {
      break
    }
    cat("\t\tsaving the moments\n")
    
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
  }
  return( preds )
}


filter_v24 = function(approx.name, XY, mod, PointerData, DataSetup.init, mu.tt = mu, Sig0 = Sig0){
  
  approx = approximations[[approx.name]]
  PointerDataComp = PointerData$DataComp
  PointerDataDecomp = PointerData$DataDeomp
  locsordOrdering = DataSetup.init$locsordOrdering
  MatchLocationOrdering = DataSetup.init$MatchLocationOrdering
  mod.init = mod
  
  Tmax = length(XY$x)
  preds = list()
  obs.aux = as.numeric(XY$y[[1]])
  
  cat(paste("\tfiltering: t=1\n"))
  
  cat("\t\tcalculating forecast moments\n")
  covmodel = GPvecchia::getMatCov(approx, as.matrix(Sig0))
  #mu.tt = 
  #mu.tt = mu
  
  cat("\t\tcalculating posterior\n")
  inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
  vals = createUcppM(ptrs, inds, cov.vals)
  L.tt =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  
  
  #preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt,
  #                                               likelihood_model = data.model, covmodel = covmodel,
  #                                               covparms = NULL, likparms = lik.params, return_all = TRUE, mod = mod)
  
  #cat("\t\tsaving the moments\n")
  #L.tt  = getLtt(approx, preds.aux)#getL00(approx, covmodel, locs) #####plesase check this
  #mu.tt = matrix(preds.aux$mean, ncol = 1)
  #preds[[1]] = list(state = mu.tt, W = preds.aux$W)#, V = preds.aux$V)
  
  if (Tmax == 1) { 
    return( preds )
  } 
  
  for (t in 1:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    if(mod.init == "standard")
    {
      cat("\t\tcalculating gradient...\n")
      Et = (Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering]
      Fmat = Et %*% L.tt
      Fmat.mod = filterFmat(Fmat = Fmat, PointerData = PointerDataComp)
      covmodel = Fmat.mod
      covmodel = GPvecchia::getMatCov(approx, as.matrix(Fmat.mod %*% Matrix::t(Fmat.mod) + Sigt))
    }else if(mod.init == "noevolerror"){
      
      storetempmat = (apply(decompressOperator_v4(L.tt, action = "compress", PointerData = PointerDataComp, M = M.multires, J = J, r = r)[locsordOrdering, ], 2, evolFun))[MatchLocationOrdering, ]
      L.tt.aux = decompressOperator_v4(storetempmat, action = "decompress", PointerData = PointerDataDecomp, M = M.multires, J = J, r = r)
      covmodel = GPvecchia::getMatCov(approx, as.matrix(L.tt.aux %*% t(L.tt.aux)))
    }else{
      
      stop("Check the `mod' variable.")
    }
    
    
    cat("\t\tcalculating forecast moments\n")
    forecast = (evolFunMu(mu.tt[locsordOrdering]))[MatchLocationOrdering]
    
    
    cat("\t\tcalculating posterior\n")
    preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                        likelihood_model = data.model, covmodel = covmodel,
                                        covparms = covparms, likparms = lik.params, return_all = TRUE, mod = "standard")
    cat("\t\tsaving the moments\n")
    
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
  }
  return( preds )
}


filter_v3 = function(approx.name, XY, mod, PointerData, DataSetup.init, recompute.iter = NULL, mu.tt = mu, Sig0 = Sig0, rotate.times = NULL, is.rotate = TRUE){
  
  if(rotate.times == 0 || is.null(rotate.times) || mod == "standard") is.rotate = FALSE
  if(!is.rotate) cat("\t\tThe co-ordinates are not being rotated.\n")
  approx = approximations[[approx.name]]
  PointerDataComp = PointerData$DataComp
  PointerDataDecomp = PointerData$DataDeomp
  LocationData.new = DataSetup.init
  locsordOrdering = DataSetup.init$locsordOrdering
  MatchLocationOrdering = DataSetup.init$MatchLocationOrdering
  
  Tmax = length(XY$x)
  preds = list()
  locsordOrderingMat = array(dim = c(length(locsordOrdering), Tmax))
  obs.aux = as.numeric(XY$y[[1]])
  
  #cat(paste("\tfiltering: t=1\n"))
  
  #cat("\t\tcalculating forecast moments\n")
  covmodel = GPvecchia::getMatCov(approx, as.matrix(Sig0))
  
  #cat("\t\tcalculating posterior\n")
  inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
  vals = createUcppM(ptrs, inds, cov.vals)
  L.tt =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  
   
  
  for (t in 1:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    cat("\t\tcalculating forecast moments\n")
    forecast = (evolFun(mu.tt[locsordOrdering]))
    
    if(mod == "standard")
    {
      cat("\t\tcalculating gradient...\n")
      Et = (Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering]
      Fmat = Et %*% L.tt
      covmodel = GPvecchia::getMatCov(approx, as.matrix(Fmat %*% Matrix::t(Fmat) + Sigt))
    }else if(mod == "noevolerror"){
      
      if((is.rotate) && (t %% 4 == 0)){
        
        cat("\t\tcalculating the gradient at rotating location\n")
        Et = (Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], 
                                           K, M, dt, Force)))
        #L.tt.aux = apply(L.tt[locsordOrdering, locsordOrdering], 2, evolFun)
        
        LocationData.new = rotateLocations(LocationData = LocationData.new, times = rotate.times)
        
        MatchLocationOrdering = LocationData.new$MatchLocationOrdering
        locsordOrdering = LocationData.new$locsordOrdering
        #L.tt.aux = L.tt.aux[MatchLocationOrdering, MatchLocationOrdering]
        
        Et = Et[MatchLocationOrdering, MatchLocationOrdering]
        
        
        
        Fmat = Et %*% 
          (L.tt[LocationData.new$rotatedOrdering, LocationData.new$rotatedOrdering])
        covmodel = GPvecchia::getMatCov(approx, #as.matrix(tcrossprod(L.tt.aux)))
                                        as.matrix(tcrossprod(Fmat)))
        
        inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
        ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
        cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
        vals = createUcppM(ptrs, inds, cov.vals)
        covmodel =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
      }
      else{
        
        storetempmat = (apply(decompressOperator_v4(L.tt[locsordOrdering, ], action = "compress", PointerData = PointerDataComp, M = M.multires, J = J, r = r), 2, evolFun))[MatchLocationOrdering, ]
        L.tt.aux = decompressOperator_v4(storetempmat, action = "decompress", PointerData = PointerDataDecomp, M = M.multires, J = J, r = r)
        covmodel = L.tt.aux 
      }
    }else{
      
      stop("Check the `mod' variable.")
    }
    
    forecast = forecast[MatchLocationOrdering]
    
    
    cat("\t\tcalculating posterior\n")
    preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                        likelihood_model = data.model, covmodel = covmodel,
                                        covparms = covparms, likparms = lik.params, return_all = TRUE, mod = mod)
    cat("\t\tsaving the moments\n")
    
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
    locsordOrderingMat[ ,t] = locsordOrdering
  }
  preds[[Tmax + 1]] = locsordOrderingMat
  return( preds )
}




filter_v4 = function(approx.name, XY, mod, PointerData, DataSetup.init, recompute.iter = NULL, mu.tt = mu, Sig0 = Sig0, rotate.times = NULL, is.rotate = TRUE){
  
  if(rotate.times == 0 || is.null(rotate.times) || mod == "standard") is.rotate = FALSE
  if(!is.rotate) cat("\t\tThe co-ordinates are not being rotated.\n")
  approx = approximations[[approx.name]]
  PointerDataComp = PointerData$DataComp
  PointerDataDecomp = PointerData$DataDeomp
  LocationData.new = DataSetup.init
  locsordOrdering = DataSetup.init$locsordOrdering
  MatchLocationOrdering = DataSetup.init$MatchLocationOrdering
  
  Tmax = length(XY$x)
  preds = list()
  locsordOrderingMat = array(dim = c(length(locsordOrdering), Tmax))
  obs.aux = as.numeric(XY$y[[1]])
  
  #cat(paste("\tfiltering: t=1\n"))
  
  #cat("\t\tcalculating forecast moments\n")
  covmodel = GPvecchia::getMatCov(approx, as.matrix(Sig0))
  
  #cat("\t\tcalculating posterior\n")
  inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
  vals = createUcppM(ptrs, inds, cov.vals)
  L.tt =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  
  
  
  for (t in 1:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    cat("\t\tcalculating forecast moments\n")
    forecast = (evolFun(mu.tt[locsordOrdering]))
    
    if(mod == "standard")
    {
      cat("\t\tcalculating gradient...\n")
      Et = (Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering]
      Fmat = Et %*% L.tt
      covmodel = GPvecchia::getMatCov(approx, as.matrix(Fmat %*% Matrix::t(Fmat) + Sigt))
    }else if(mod == "noevolerror"){
      
      if((is.rotate) && (t %% 4 == 0)){
        
        cat("\t\tcalculating the gradient at rotating location\n")
        Et = (Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], 
                                           K, M, dt, Force)))
        #L.tt.aux = apply(L.tt[locsordOrdering, locsordOrdering], 2, evolFun)
        
        LocationData.new = rotateLocations(LocationData = LocationData.new, times = rotate.times)
        
        MatchLocationOrdering = LocationData.new$MatchLocationOrdering
        locsordOrdering = LocationData.new$locsordOrdering
        #L.tt.aux = L.tt.aux[MatchLocationOrdering, MatchLocationOrdering]
        
        Et = Et[MatchLocationOrdering, MatchLocationOrdering]
        
        
        
        Fmat = Et %*% 
          (L.tt[LocationData.new$rotatedOrdering, LocationData.new$rotatedOrdering])
        covmodel = GPvecchia::getMatCov(approx, #as.matrix(tcrossprod(L.tt.aux)))
                                        as.matrix(tcrossprod(Fmat)))
        
        inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
        ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
        cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
        vals = createUcppM(ptrs, inds, cov.vals)
        covmodel =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
      }
      else{
        if(t == 1){
          
          cat("\t\tcalculating gradient...\n")
          Et = (Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering]
          Fmat = Et %*% L.tt
          covmodel = GPvecchia::getMatCov(approx, as.matrix(Fmat %*% Matrix::t(Fmat) + Sigt))
          inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
          ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
          cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
          vals = createUcppM(ptrs, inds, cov.vals)
          covmodel =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
        }else{
          
          storetempmat = (apply(decompressOperator_v4(L.tt[locsordOrdering, ], action = "compress", PointerData = PointerDataComp, M = M.multires, J = J, r = r), 2, evolFun))[MatchLocationOrdering, ]
          L.tt.aux = decompressOperator_v4(storetempmat, action = "decompress", PointerData = PointerDataDecomp, M = M.multires, J = J, r = r)
          covmodel = L.tt.aux
        }
      }
    }else{
      
      stop("Check the `mod' variable.")
    }
    
    forecast = forecast[MatchLocationOrdering]
    
    
    cat("\t\tcalculating posterior\n")
    preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                        likelihood_model = data.model, covmodel = covmodel,
                                        covparms = covparms, likparms = lik.params, return_all = TRUE, mod = mod)
    cat("\t\tsaving the moments\n")
    
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
    locsordOrderingMat[ ,t] = locsordOrdering
  }
  preds[[Tmax + 1]] = locsordOrderingMat
  return( preds )
}


filter_advdifflin = function(approx.name, XY, mod, PointerData, DataSetup.init, mu.tt = mu, Sig0 = Sig0, alpha = NA, beta = NA){
  
  approx = approximations[[approx.name]]
  PointerDataComp = PointerData$DataComp
  PointerDataDecomp = PointerData$DataDeomp
  locsordOrdering = DataSetup.init$locsordOrdering
  MatchLocationOrdering = DataSetup.init$MatchLocationOrdering
  
  Et = #Matrix::sparseMatrix(i = 1:length(locsordOrdering), 
      #                                   j = 1:length(locsordOrdering),
      #                                  x = rep(2, length(locsordOrdering)))
    diffAdvOp1DLin(n_G = length(locsordOrdering), alpha = alpha, beta = beta)
  
  Tmax = length(XY$x)
  preds = list()
  obs.aux = as.numeric(XY$y[[1]])
  
  cat(paste("\tfiltering: t=1\n"))
  
  cat("\t\tcalculating forecast moments\n")
  covmodel = GPvecchia::getMatCov(approx, as.matrix(Sig0))
  #mu.tt = 
  #mu.tt = mu
  
  cat("\t\tcalculating posterior\n")
  inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
  vals = createUcppM(ptrs, inds, cov.vals)
  L.tt =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  
  
  #preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt,
  #                                               likelihood_model = data.model, covmodel = covmodel,
  #                                               covparms = NULL, likparms = lik.params, return_all = TRUE, mod = mod)
  
  #cat("\t\tsaving the moments\n")
  #L.tt  = getLtt(approx, preds.aux)#getL00(approx, covmodel, locs) #####plesase check this
  #mu.tt = matrix(preds.aux$mean, ncol = 1)
  #preds[[1]] = list(state = mu.tt, W = preds.aux$W)#, V = preds.aux$V)
  
  if (Tmax == 1) { 
    return( preds )
  } 
  
  for (t in 1:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    if(mod == "standard")
    {
      cat("\t\tcalculating gradient...\n")
      #Et = (Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering]
      Fmat = Et %*% L.tt
      Fmat.mod = filterFmat(Fmat = Fmat, PointerData = PointerDataComp)
      covmodel = GPvecchia::getMatCov(approx, as.matrix(Fmat.mod %*% Matrix::t(Fmat.mod) #+ Sigt
                                                        ))
    }else if(mod == "noevolerror"){
      
      storetempmat = (Et %*% (decompressOperator_v4(L.tt, action = "compress", PointerData = PointerDataComp, M = M.multires, J = J, r = r)[locsordOrdering, ]))[MatchLocationOrdering, ]
      L.tt.aux = decompressOperator_v4(storetempmat, action = "decompress", PointerData = PointerDataDecomp, M = M.multires, J = J, r = r)
      covmodel = L.tt.aux
    }else{
      
      stop("Check the `mod' variable.")
    }
    
    
    cat("\t\tcalculating forecast moments\n")
    forecast = matrix((Et %*%(mu.tt[locsordOrdering]))[MatchLocationOrdering], ncol = 1)
    
    
    cat("\t\tcalculating posterior\n")
    preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                        likelihood_model = data.model, covmodel = covmodel,
                                        covparms = covparms, likparms = lik.params, return_all = TRUE, mod = mod)
    cat("\t\tsaving the moments\n")
    
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
  }
  return( preds )
}



filterLorenz_v5 = function(approx.name, XY, mod, PointerData, DataSetup.init, mu.tt = mu, Sig0 = Sig0, Sigt = NULL, r = NULL){
  
  approx = approximations[[approx.name]]
  PointerDataComp = PointerData$DataComp
  PointerDataDecomp = PointerData$DataDeomp
  locsordOrdering = DataSetup.init$locsordOrdering
  MatchLocationOrdering = DataSetup.init$MatchLocationOrdering
  
  Tmax = length(XY$x)
  preds = list()
  obs.aux = as.numeric(XY$y[[1]])
  
  cat(paste("\tfiltering: t=1\n"))
  
  cat("\t\tcalculating forecast moments\n")
  covmodel = getMatCov(approx, as.matrix(Sig0))
  #mu.tt = 
  #mu.tt = mu
  
  cat("\t\tcalculating posterior\n")
  #inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
  #ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  #cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
  #vals = createUcppM(ptrs, inds, cov.vals)
  #L.tt =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  
  preds.aux = tryCatch(calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt,
                                               likelihood_model = data.model, covmodel = covmodel,
                                               covparms = covparms, likparms = lik.params, return_all = TRUE, mod = "standard"),
                       error = function(e) 1)
  
  if((!is.list(preds.aux)) || any(is.na(preds.aux$mean)))
  {
    return(preds)
  }
  
  
  cat("\t\tsaving the moments\n")
  L.tt  = getLtt(approx, preds.aux)#getL00(approx, covmodel, locs) #####plesase check this
  mu.tt = matrix(preds.aux$mean, ncol = 1)
  preds[[1]] = list(state = mu.tt, W = preds.aux$W)#, V = preds.aux$V)
  
  if (Tmax == 1) { 
    return( preds )
  } 
  
  for (t in 2:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    if(mod == "standard")
    {
      cat("\t\tcalculating gradient...\n")
      Et = tryCatch((Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering],
                    error = function(e) 1)
      if(isTRUE(any(is.nan(c(as.matrix(Et)))))|| isTRUE(is.numeric(Et))){
        
        break
      }
      
      Fmat = Et %*% L.tt
      Fmat.mod = filterFmat(Fmat = Fmat, PointerData = PointerDataComp)
      covmodel = getMatCov(approx, as.matrix(Fmat.mod %*% Matrix::t(Fmat.mod)))
    }else if(mod == "noevolerror"){
      
      storetempmat = (apply(decompressOperator_v4(L.tt, action = "compress", PointerData = PointerDataComp, M = M.multires, J = J, r = r)[locsordOrdering, ], 2, evolFun))[MatchLocationOrdering, ]
      L.tt.aux = decompressOperator_v4(storetempmat, action = "decompress", PointerData = PointerDataDecomp, M = M.multires, J = J, r = r)
      covmodel = L.tt.aux
    }else if(mod == "redrank"){
     
      object.temp = UKFmeanmatredrank(L.tt[locsordOrdering, locsordOrdering], 
                                      mu.tt[locsordOrdering], 
                                      evolFun, r)
      covmodel = getMatCov(approx, (as.matrix(object.temp$covmat
                                             [MatchLocationOrdering, MatchLocationOrdering]) + 
                             diag(rep(10^-08, length(obs.aux)))))
      forecast = (object.temp$mean)[MatchLocationOrdering]
    }else{
      
      stop("Check the `mod' variable.")
    }
    
    
    cat("\t\tcalculating forecast moments\n")
    if(mod %in% c("standard", "noevolerror")){
      forecast = (evolFunMu(mu.tt[locsordOrdering]))[MatchLocationOrdering] 
    }
    
    
    cat("\t\tcalculating posterior\n")
    if(mod == "redrank"){
      preds.aux = tryCatch(calculate_posterior_VL(obs.aux, approx, prior_mean = forecast,
                                                   likelihood_model = data.model, covmodel = covmodel,
                                                   covparms = covparms, likparms = lik.params, return_all = TRUE, mod = "standard"),
                           error = function(e) 1)
    }else{
      
      preds.aux = tryCatch(calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                                   likelihood_model = data.model, covmodel = covmodel,
                                                   covparms = covparms, likparms = lik.params, return_all = TRUE, mod = mod),
                           error = function(e) 1)
    }
    
    if((!is.list(preds.aux)) || any(is.na(preds.aux$mean)))
    {
      break
    }
    cat("\t\tsaving the moments\n")
    
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
  }
  return( preds )
}


filterLorenz_v1.add.error = function(approx.name, XY, mod, PointerData, DataSetup.init, mu.tt = mu, Sig0 = Sig0, Sigt, r = NULL){
  
  approx = approximations[[approx.name]]
  PointerDataComp = PointerData$DataComp
  PointerDataDecomp = PointerData$DataDeomp
  locsordOrdering = DataSetup.init$locsordOrdering
  MatchLocationOrdering = DataSetup.init$MatchLocationOrdering
  
  R.total = sum(r)
  
  Tmax = length(XY$x)
  preds = list()
  obs.aux = as.numeric(XY$y[[1]])
  
  cat(paste("\tfiltering: t=1\n"))
  
  cat("\t\tcalculating forecast moments\n")
  covmodel = getMatCov(approx, as.matrix(Sig0))
  #mu.tt = 
  #mu.tt = mu
  
  cat("\t\tcalculating posterior\n")
  #inds = Filter(function(i) !is.na(i), as.vector(t(approx$U.prep$revNNarray - 1)))
  #ptrs = c(0, cumsum(apply(approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  #cov.vals = Filter(function(i) !is.na(i), c(t(covmodel)))
  #vals = createUcppM(ptrs, inds, cov.vals)
  #L.tt =  Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  
  preds.aux = tryCatch(calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt,
                                               likelihood_model = data.model, covmodel = covmodel,
                                               covparms = covparms, likparms = lik.params, return_all = TRUE, mod = "standard"),
                       error = function(e) 1)
  
  if((!is.list(preds.aux)) || any(is.na(preds.aux$mean)))
  {
    return(preds)
  }
  
  if(!missing(Sigt)){
    Sigtc = Matrix::chol(as.matrix(Sigt))
  }else{
    stop("Please supply the variable Sigt.")
  }
  
  cat("\t\tsaving the moments\n")
  L.tt  = getLtt(approx, preds.aux)#getL00(approx, covmodel, locs) #####plesase check this
  mu.tt = matrix(preds.aux$mean, ncol = 1)
  preds[[1]] = list(state = mu.tt, W = preds.aux$W)#, V = preds.aux$V)
  
  if (Tmax == 1) { 
    return( preds )
  } 
  
  for (t in 2:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    if(mod == "standard")
    {
      cat("\t\tcalculating gradient...\n")
      Et = tryCatch((Matrix::Matrix(exactGradient(mu.tt[locsordOrdering], K, M, dt, Force)))[MatchLocationOrdering, MatchLocationOrdering],
                    error = function(e) 1)
      if(isTRUE(any(is.nan(c(as.matrix(Et)))))|| isTRUE(is.numeric(Et))){
        
        break
      }
      
      Fmat = Et %*% L.tt
      Fmat.mod = filterFmat(Fmat = Fmat, PointerData = PointerDataComp)
      covmodel = getMatCov(approx, as.matrix(Fmat.mod %*% Matrix::t(Fmat.mod)))
    }else if(mod == "noevolerror"){
      
      errormat = matrix(rnorm(R.total * nrow(Sigt)), ncol = R.total)
      storetempmat = (apply(decompressOperator_v4(L.tt, action = "compress", PointerData = PointerDataComp, M = M.multires, J = J, r = r)[locsordOrdering, ], 2, evolFun))[MatchLocationOrdering, ]
      storetempmat = as.matrix(storetempmat) +t(Sigtc) %*% errormat#crossprod(Sigtc, errormat)
      L.tt.aux = decompressOperator_v4(storetempmat, action = "decompress", PointerData = PointerDataDecomp, M = M.multires, J = J, r = r)
      covmodel = L.tt.aux
    }else if(mod == "redrank"){
      
      object.temp = UKFmeanmatredrank(L.tt[locsordOrdering, locsordOrdering], 
                                      mu.tt[locsordOrdering], 
                                      evolFun, r)
      covmodel = getMatCov(approx, (as.matrix(object.temp$covmat
                                              [MatchLocationOrdering, MatchLocationOrdering]) + 
                                      Sigt))
      forecast = (object.temp$mean)[MatchLocationOrdering]
    }else{
      
      stop("Check the `mod' variable.")
    }
    
    
    cat("\t\tcalculating forecast moments\n")
    if(mod %in% c("standard", "noevolerror")){
      forecast = (evolFunMu(mu.tt[locsordOrdering]))[MatchLocationOrdering] 
    }
    
    
    cat("\t\tcalculating posterior\n")
    if(mod == "redrank"){
      preds.aux = tryCatch(calculate_posterior_VL(obs.aux, approx, prior_mean = forecast,
                                                  likelihood_model = data.model, covmodel = covmodel,
                                                  covparms = covparms, likparms = lik.params, return_all = TRUE, mod = "standard"),
                           error = function(e) 1)
    }else{
      
      preds.aux = tryCatch(calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                                   likelihood_model = data.model, covmodel = covmodel,
                                                   covparms = covparms, likparms = lik.params, return_all = TRUE, mod = mod),
                           error = function(e) 1)
    }
    
    if((!is.list(preds.aux)) || any(is.na(preds.aux$mean)))
    {
      break
    }
    cat("\t\tsaving the moments\n")
    
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
  }
  return( preds )
}


getLtt = function(vecchia.approx, preds){
  n = nrow(vecchia.approx$locsord)
  orig.order=order(vecchia.approx$ord)
  V = preds$V
  L.tt = (Matrix::solve(Matrix::t(V), sparse=TRUE)[seq(n, 1), seq(n, 1)])[orig.order,]
  return(L.tt)
}

## simulate y given x 
simulate.y = function(x, frac.obs, lik.params){
  
  n = nrow(x)
  n.obs = round(n*frac.obs)
  obs.inds = sample(1:n, n.obs, replace = FALSE)
  data.model = lik.params["data.model"]
  # simulate data
  if(data.model=='poisson'){
    y.obs = rpois(n.obs, exp(x[obs.inds]))
  } else if(data.model=='logistic'){
    y.obs = rbinom(n.obs,1,prob = exp(x[obs.inds])/(1+exp(x[obs.inds])))
  } else if(data.model=='gamma'){
    #default_lh_params = list("alpha"=2, "sigma"=sqrt(.1), "beta"=.9, "phi"=1.5)
    #z = rgamma(n.obs, shape = default_lh_params$alpha, rate = default_lh_params$alpha*exp(-y[obs.inds]))
    y.obs = rgamma(n.obs, shape = lik.params[["alpha"]], rate = lik.params[["alpha"]]*exp(-x[obs.inds]))
    
  } else if(data.model == 'gauss'){
    y.obs = rnorm(n.obs, mean = x[obs.inds], sd=lik.params[["sigma"]])
    
    
  } else {
    print('Error: Distribution not implemented yet.')
  }
  y = rep(NA, n)
  y[obs.inds] = y.obs
  return(y)
}



## simulate x
simulate.xy = function(x0, E, Q, frac.obs, lik.params, Tmax, seed=NULL, sig2=1, smooth = 0.5, range = 1, locs = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  n = nrow(x0);
  x = list(); y = list()
  x[[1]] = x0
  y[[1]] = simulate.y(x0, frac.obs, lik.params)
  
  if (Tmax > 1) { 
    
    if (!is.null(Q) && any(Q)) {
      Qc = Matrix::chol(Q)
    } 
    
    for (t in 2:Tmax) {
      if (isTRUE(sig2 > 0 || (!is.null(Q) && sum(abs(Q))>0))) {
        if (!is.null(Q)) {
          w =  t(Qc) %*% matrix(rnorm(n), ncol = 1)
        } else {
          w = matrix(sig2*RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range),
                                                   x = locs[,1], y = locs[,2], spConform = FALSE), ncol=1)
        } 
      } else {
        w = matrix(rep(0, n), ncol=1)
      }
      
      x[[t]] = E(X = x[[t - 1]]) + w
      y[[t]] = simulate.y(x[[t]], frac.obs, lik.params)
    } 
  }
  
  return(list(x = x, y = y))
  
}