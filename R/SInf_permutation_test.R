SInf_permutation_test <- function(mod,Tx,Y,nc,ncs = 2:5, nperms=1e5){
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
  
  # permutations needed for model selection
  omnibus_stats = matrix(0,nperms,length(ncs))
  
  n = nrow(Y)
  mods = list()
  mod1 = NULL
  for (c in 1:length(ncs)){
    mods[[c]] = SV(Tx,Y,nc=ncs[c])
    mods[[c]]$abs_sum = sum(abs(mods[[c]]$MR))
    for (p in 1:nperms){
      cat('\r',p)
      perm = sample(1:n,n,replace=FALSE)
      mod1 = SV(Tx,Y[perm,],ee = mods[[c]]$eigen,nc=ncs[c])
      
      # omnibus
      omnibus_stats[p,c] = sum(abs(mod1$MR))
    }
  }
  
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
  
  abs_diff = tx_pairs$pval
  tx_pairs$diff = tx_pairs$pval
  for (j in 1:nc){
    abs_diff[,j] = as.numeric(dist(mod$MR[,j]))
    
    z = -outer(mod$MR[,j],mod$MR[,j],'-'); 
    tx_pairs$diff[,j] = z[lower.tri(z)]
  }
  
  TP = rank(-factors$abs_sum)
  FP = rep(0,nc)
  for (p in 1:nperms){ # outer permutation
    cat('\r',p)
    
    ### MODEL SELECTION
    q = 0
    while ( q != nc){
      ps_omnibus = c()
      mod1 = list()
      for (c in 1:length(ncs)){
        perm = sample(1:n,n,replace=FALSE)
        mod1[[c]] = SV(Tx,Y[perm,],ee = mods[[c]]$eigen,nc=ncs[c])
        ps_omnibus = c(ps_omnibus, sum(sum(abs(mod1[[c]]$MR)) < omnibus_stats[,c])/nperms)
      }
      q = ncs[which(ps_omnibus == min(ps_omnibus))[1]] # model selection performed
    }
    
    qi = which(ncs==q) # convert to index
    
    # omnibus
    omnibus$pval = omnibus$pval + (omnibus$abs_sum <= sum(abs(mod1[[qi]]$MR)))
    
    # post hoc of factors
    factors$pval = factors$pval + (factors$abs_sum <= colSums(abs(mod1[[qi]]$MR)) )
    for (j in 1:nc){
      FP[j] = FP[j] + sum(factors$abs_sum[j] <= colSums(abs(mod1[[qi]]$MR)) ) # FDR 5%
    }
    factors$pFWER = factors$pFWER + (factors$abs_sum <= max(colSums(abs(mod1[[qi]]$MR))) ) #FWER
    
    # post hoc of treatment pairs
    # range = apply(mod1[[qi]]$MR,2,max) - apply(mod1[[qi]]$MR,2,min) # range statistic for FWER
    for (j in 1:nc){
      diff1 = as.numeric(dist(mod1[[qi]]$MR[,j]))
      tx_pairs$pFWER[,j] = tx_pairs$pFWER[,j] + (abs_diff[,j] <= max(diff1))
      tx_pairs$pval[,j] = tx_pairs$pval[,j] + (abs_diff[,j] <= diff1)
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