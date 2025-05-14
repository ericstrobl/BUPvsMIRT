generate_synth <- function(nTx = 5, nX=10, nY = 20, nsamps=1000){
  
  Tx = rep(1:nTx,length.out=nsamps)
  X = matrix(rnorm(nsamps*nX),nsamps,nX)
  Y = matrix(rnorm(nsamps*nY),nsamps,nY)
  
  return(list(X=X, Y=Y, Tx = Tx))
  
}