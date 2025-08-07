num_components <- function(Tx,Y,ncs=2:5){

  vars = c()
  for (n in ncs){
    mod = SV(Tx,Y,nc=n)
    values = permutation_testing(mod,Tx,Y,1000,nc=n)$omnibus$pval
    vars = c(vars,values)
  } 
  
  return(ncs[which(vars==min(vars))[1]])


}
