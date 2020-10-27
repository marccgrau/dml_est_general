
pi1 = function(X, w = 0){
  if (ncol(X) >= 5){
    y = w + X[,2] + X[,4] + X[,3]*X[,4] - X[,5]^2 + rnorm(1,0,1)
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

pi2 = function(X, w = 0){
  if (ncol(X) >= 5){
    y = w + X[,3] + X[,2]*X[,3] + X[,1]*X[,4]*X[,6] + X[,2]*X[,5]*X[,7] + rnorm(1,0,1)
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}


pi3 = function(X, w){
  if (ncol(X) >= 1){
    
    y = w + X[,1] + X[,2]^2 + X[,3]^3 + X[,5]^3 - 2*X[,6]^4 + 0.5  + rnorm(1,0,1)
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

test = rep(NA, 10000)
for (i in 1:10000){
  test[i] = mean(ind_func(X, 0, pi3))
}

pi4 = function(X, w){
  if (ncol(X) >= 1){
    
    y = w  + rnorm(1,0,1)
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

test = rep(NA, 10000)
for (i in 1:10000){
  test[i] = mean(propensity_func(X, pi1, -log(9)))
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