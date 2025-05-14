permutation_testing <- function(mod,Tx,Y,nperms=1e5,nc=length(unique(Tx))){
  # omnibus test, post hoc test of factors, and post hoc test of treatment pairs by permutations
  # Inputs:
  #   mod = original model from SV function
  #   Tx = vector of treatments for n samples
  #   Y = n by p matrix of p rating scale items
  #   nperms = number of permutations (default: 100,000)
  #
  # Outputs:
  #   Omnibus:
  #     abs_sum = absolute sum statistic
  #     pval = p-value
  #   Post-hoc for factors:
  #     abs_sum = absolute sum statistic for each factor
  #     pval = uncorrected p-values
  #     pFWER = family-wise error rate corrected p-values (Holm method)
  #     pFDR = positive false discovery rate corrected p-values (Storey method)
  #   Post-hoc for treatment pairs:
  #     diff = difference statistic for each treatment pair
  #     pval = uncorrected p-values
  #     pFWER = family-wise error rate corrected p-values (Tukey's range / maxT method)
  #
  #
  # written by Eric V. Strobl, 12/3/2024
  
  require(qvalue)
  
  tx = unique(Tx)
  nt = length(tx)
  
  omnibus = list()
  factors = list()
  tx_pairs = list()
  
  omnibus$pval = 0
  omnibus$abs_sum = sum(abs(mod$MR))
  
  factors$pval = rep(0,nc)
  factors$pFDR = factors$pval
  factors$pFWER = factors$pval
  factors$abs_sum = colSums(abs(mod$MR))
  
  tx_pairs$pval = matrix(0,(nt^2-nt)/2,nc) # uncorrected p-values
  tx_pairs$pFWER = tx_pairs$pval
  
  n = nrow(Y)
  
  # compute difference and absolute difference statistics for original samples
  abs_diff = tx_pairs$pval
  tx_pairs$diff = tx_pairs$pval
  # ord = tx_pairs$pval
  for (j in 1:nc){
    abs_diff[,j] = as.numeric(dist(mod$MR[,j]))
    # ord[,j] = order(abs_diff[,j],decreasing=TRUE)
    
    z = outer(mod$MR[,j],mod$MR[,j],'-'); 
    tx_pairs$diff[,j] = z[lower.tri(z)]
  }
  
  
  # perform permutation statistic
  TP = rank(-factors$abs_sum)
  FP = rep(0,nc)
  for (p in 1:nperms){
    cat('\r',p)
    perm = sample(1:n,n,replace=FALSE)
    mod1 = SV(Tx,Y[perm,],ee = mod$eigen,nc=nc)
    
    # omnibus
    omnibus$pval = omnibus$pval + (omnibus$abs_sum <= sum(abs(mod1$MR)))
    
    # post hoc of factors
    factors$pval = factors$pval + (factors$abs_sum <= colSums(abs(mod1$MR)))
    for (j in 1:nc){
      FP[j] = FP[j] + sum(factors$abs_sum[j] <= colSums(abs(mod1$MR)) ) # FDR
    }
    factors$pFWER = factors$pFWER + (factors$abs_sum <= max(colSums(abs(mod1$MR))) ) #FWER
    
    # post hoc of treatment pairs
    # range = apply(mod1$MR,2,max) - apply(mod1$MR,2,min) # range statistic for FWER
    for (j in 1:nc){
      dis1 = as.numeric(dist(mod1$MR[,j]))
      # ranges1 = sort(dis1,decreasing=TRUE) # successive maxima
      # tx_pairs$pFWER[ord[,j],j] = tx_pairs$pFWER[ord[,j],j] + (abs_diff[ord[,j],j] <= ranges1)
      tx_pairs$pFWER[,j] = tx_pairs$pFWER[,j] + (abs_diff[,j] <= max(dis1))
      tx_pairs$pval[,j] = tx_pairs$pval[,j] + (abs_diff[,j] <= dis1)
    }
  }
  omnibus$pval = omnibus$pval/nperms
  
  factors$pval = factors$pval/nperms
  pi0 = (nc-1)/nc # because omnibus was rejected
  factors$pFDR = FP/TP
  factors$pFDR = pmin(pi0*(factors$pFDR/nperms),1)
  factors$pFWER = factors$pFWER/nperms
  
  tx_pairs$pval = tx_pairs$pval/nperms
  tx_pairs$pFWER = tx_pairs$pFWER/nperms
  
  # for (j in 1:nc){ # enforce monotonicity
  #   tx_pairs$pFWER[ord[,j],j] = cummax(tx_pairs$pFWER[ord[,j],j])
  # }
  
  # label tx pairs
  labels = c()
  for (c in 1:(length(tx)-1)){
    for (d in (c+1):length(tx)){
      labels = c(labels, paste(as.character(tx[c]), "-", as.character(tx[d]),sep=""))
    }
  }
  rownames(tx_pairs$pval) = labels
  rownames(tx_pairs$pFWER) = labels
  rownames(tx_pairs$diff) = labels
  
  return( list(omnibus = omnibus, factors = factors, tx_pairs = tx_pairs) )
}