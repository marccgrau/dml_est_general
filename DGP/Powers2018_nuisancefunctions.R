f1 = function(X){
  y = 0
  return(y)
}

f2 = function(X){
  if (ncol(X) >= 1){
    I_1 <- 1*(X[,1]>1)
    
    y = 5*I_1 - 5
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

f3 = function(X){
  if (ncol(X) >= 1){
    
    y = 2*X[,1] - 4
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

f4 = function(X){
  if(ncol(X) >= 6){
    
    y = X[,2]*X[,4]*X[,6] + 2*X[,2]*X[,4]*(1 - X[,6]) + 3*X[,2]*(1 - X[,4])*X[,6] +
      4*X[,2]*(1 - X[,4])*(1 - X[,6]) + 5*(1 - X[,2])*X[,4]*X[,6] + 6*(1 - X[,2])*X[,4]*(1 - X[,6]) +
      7*(1 - X[,2])*(1 - X[,4])*X[,6] + 8*(1 - X[,2])*(1 - X[,4])*(1 - X[,6])
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

f5 = function(X){
  if(ncol(X) >= 9){
    y = X[,1] + X[,3] + X[,5] + X[,7] + X[,8] + X[,9] - 2
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

f6 = function(X){
  if(ncol(X) >= 9){
    I_1 <- 1*(X[,1]>1)
    I_3 <- 1*(X[,3]>0)
    I_5 <- 1*(X[,5]>1)
    I_7 <- 1*(X[,7]>0)
    
    y = 4*I_1*I_3 + 4*I_5*I_7 + 2*X[,8]*X[,9]
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

f7 = function(X){
  if(ncol(X) >= 9){
    
    y = 0.5*(X[,1]^2 + X[,2] + X[,3]^2 + X[,4] + X[,5]^2 + X[,6] + X[,7]^2 + X[,8] + X[,9]^2 - 11)
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

f8 = function(X){
  if(ncol(X) >= 5){
    
    y = (1/sqrt(2))*(f4(X) + f5(X))
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}