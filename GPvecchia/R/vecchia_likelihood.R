#' evaluation of the likelihood
#'
#' @param z the observed data
#' @param vecchia.approx a vecchia object as generated by vecchia_specify()
#' @param covparms covariance parameters as a vector
#' @param nuggets either a single (constant) nugget or a vector of nugget terms for the observations
#' @param covmodel covariance model, 'matern' by default
#'
#' @return (multivariate normal) loglikelihood implied by the Vecchia approximation
#' @examples
#' z=rnorm(5); locs=matrix(1:5,ncol=1); vecchia.approx=vecchia_specify(locs,m=3)
#' vecchia_likelihood(z,vecchia.approx,covparms=c(1,2,.5),nuggets=.2)
#' @export
vecchia_likelihood=function(z,vecchia.approx,covparms,nuggets,covmodel='matern') {

  if(vecchia.approx$cond.yz=='zy')
    warning("cond.yz='zy' will produce a poor likelihood approximation. Use 'SGV' instead.")

  # remove NAs in data and U
  removeNAs()
    
  # create the U matrix
  U.obj=createU(vecchia.approx,covparms,nuggets,covmodel)

  # compute the loglikelihood
  vecchia_likelihood_U(z,U.obj)
}

# ## remove missing data (NA)
# removeNAs.old=function(){ # overwrites z and U.obj
#   p = parent.frame()
#   if(any(is.na(p$z))){
#     #ind.na=(((1:nrow(p$U.obj$U))[!p$U.obj$latent])[p$U.obj$ord.z])[is.na(p$z)]
#     ind.na=(((1:nrow(p$U.obj$U))[!p$U.obj$latent])[p$U.obj$ord.z])[is.na(p$z[p$U.obj$ord.z])]
#     if(any(apply(p$U.obj$U[,ind.na,drop=FALSE],2,Matrix::nnzero)>2)) stop(
#       'NA data is conditioned upon')
#     p$U.obj$U = p$U.obj$U[-ind.na,-ind.na]
#     p$U.obj$latent = p$U.obj$latent[-ind.na]
#     p$U.obj$ord.z = order(order(p$U.obj$ord.z[p$U.obj$ord.z %in% which(!is.na(p$z))]))
#     p$z = p$z[!is.na(p$z)]
#   }
# }


removeNAs=function(){ # overwrites z and U.obj
    p = parent.frame()
    if(any(is.na(p$z))){

        if(length(p$nuggets)<length(p$z)) {
            new.nuggets = rep(0, length(p$z))
            new.nuggets[!is.na(p$z)] = p$nuggets
            p$nuggets = new.nuggets
        }
        
        p$nuggets[is.na(p$z)] = stats::var(p$z,na.rm=TRUE)*1e8
        p$z[is.na(p$z)] = mean(p$z,na.rm=TRUE)
    }
}


## evaluate vecchia likelihood based on U

vecchia_likelihood_U=function(z,U.obj) {

  ### output: loglikelihood (for z)
  U=U.obj$U
  latent=U.obj$latent
  zord=z[U.obj$ord.z]

  # constant
  const=sum(!latent)*log(2*pi)

  # numerator
  z1=Matrix::crossprod(U[!latent,],zord)
  quadform.num=sum(z1^2)
  logdet.num=-2*sum(log(Matrix::diag(U)))

  # denominator
  if(sum(latent)==0){ # no latents -> denominator not needed

    logdet.denom=quadform.denom=0

  } else {  # if latents, need denominator

    U.y=U[latent,]
    z2=as.numeric(U.y%*%z1)
    V.ord=U2V(U.obj)
    z3=Matrix::solve(V.ord,rev(z2),system='L')
    quadform.denom=sum(z3^2)
    logdet.denom=-2*sum(log(Matrix::diag(V.ord)))

  }

  # putting everything together
  neg2loglik=logdet.num-logdet.denom+quadform.num-quadform.denom+const
  loglik=-neg2loglik/2
  return(loglik)

}


## function to reverse-order a matrix
revMat=function(mat) mat[nrow(mat):1,ncol(mat):1,drop=FALSE]
