SV <- function(Tx,Y, nc=length(unique(Tx)), ee=NULL){
  # Supervised Varimax algorithm
  # Inputs:
  #   Tx = vector of treatments for n samples
  #   Y = n by p matrix of p rating scale items
  #   nc = number of components/factors; default is set to number of unique treatments
  #   ee = eigendecomposition of cor(Y) from previous output of SV; useful for minimizing run-time of permutation tests
  # Outputs - a list with:
  #   MR = causal effects from treatments to factors
  #   RtW = causal effects from factors to treatments
  #   tx_effects = (MR) %*% (RtW), or the causal effects from treatments to factors
  #   optimal_outcomes = rotated factors
  #   R = the rotation matrix found by Varimax
  #   sgn = sign flip done to ensure sign determinancy
  #   io = re-ordering of optimal outcomes done to ensure permutation determinancy
  #   eigen = results of eigendecomposition
  #
  #
  # written by Eric V. Strobl, 12/3/2024
  
  
  # center and weight each column equally
  mY = apply(Y,2,mean)
  sdY = apply(Y,2,sd)
  Y = normalizeData(Y)
  
  # convert tx vector to binary assignments
  txs = unique(Tx)
  t_idx = matrix(0,nrow(Y),length(txs))
  coef = matrix(0,length(txs),ncol(Y))
  for (t in 1:length(txs)){
    t_idx[which(Tx == txs[t]),t] = 1
  }
  
  # regress items on tx
  coef = lm.fit(t_idx,Y)$coefficients
  
  # discover unrotated factors
  if (is.null(ee)){
    ee = eigen(cov(Y))
  }
  
  # perform Varimax rotation
  R = my_varimax(coef%*% ee$vectors[,1:nc,drop=FALSE] %*% diag(1/sqrt(ee$values[1:nc])),normalize=FALSE)$rotmat
  
  # compute optimized factors, aka optimal outcomes
  factors = sweep(Y %*% ee$vectors[,1:nc,drop=FALSE],2,1/sqrt(ee$values[1:nc]),FUN="*") %*% R
  
  # causal effects from tx to factors
  MR = sweep(coef%*% ee$vectors[,1:nc,drop=FALSE],2, 1/(sqrt(ee$values[1:nc])),FUN="*") %*% R
  rownames(MR) = txs
  
  # sign determinancy
  sgn = c(sign(cor(factors,rowSums(Y))))

  # permutation determinancy
  io = order(cor(factors,rowSums(Y))^2,decreasing=T) # order by amount of variance explained
  io = cbind(1:nc,io)
  
  # correct for sign and permutation determinancy
  factors = sweep(factors,2,sgn,FUN="*")[,io[,2],drop=FALSE]
  MR = sweep(MR,2,sgn,FUN="*")[,io[,2],drop=FALSE]
  
  # causal effects from factors to Y
  RtW = sweep(t(R),2, sqrt(ee$values[1:nc]),FUN="*") %*% t(ee$vectors[,1:nc,drop=FALSE])
  RtW  = sweep(RtW,1,sgn,FUN="*")[io[,2],]
  
  # causal effects from tx to Y
  tx_effects = MR %*% RtW
  
  items_to_factors = (sweep(ee$vectors[,1:nc,drop=FALSE],2,1/sqrt(ee$values[1:nc]),FUN="*") %*% R)
  # items_to_factors = sweep(items_to_factors,1,sgn,FUN="*")[,io[,2],drop=FALSE]
  items_to_factors = items_to_factors[,io[,2],drop=FALSE]
  
  return(list(MR = MR, RtW = RtW, tx_effects = tx_effects, 
              optimal_outcomes = factors, R = R,
              sgn = sgn, io = io, eigen = ee,
              items_to_factors = items_to_factors, sdY = sdY, mY = mY,
              t_idx = t_idx))
}
