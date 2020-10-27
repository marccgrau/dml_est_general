# Following Setoguchi et al. (2008) with slight adjustments for more sparsity

fD = function(X, beta = beta){
  if (ncol(X) >= 7){
    I_1 = 1*(X[,1]>1)
    I_3 = 1*(X[,3]>1)
    I_5 = 1*(X[,5]>1)
    
    y = log(0.1) + beta[1]*X[,1] + beta[2]*X[,2] + beta[3]*X[,3]*X[,4] + beta[4]*X[,4]*X[,5] + beta[5]*X[,5]*X[,6]
      beta[1]*X[,1]*X[,2]*X[,3] + beta[2]*X[,2]*X[,4]*X[,6] + I_1*X[,2] + I_3*X[,4] + I_5*X[,6]
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}