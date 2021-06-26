pi1 = function(X, w = 0){
  if (ncol(X) >= 5){
    y = X[,2] + X[,4] + X[,3]*X[,4] - X[,5]^2
    
    y_scale = scale(y, center = TRUE, scale = TRUE) + w
    return(y_scale)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

pi2 = function(X, w = 0){
  if (ncol(X) >= 7){
    y = X[,3] + X[,2]*X[,3] + X[,1]*X[,4]*X[,6] + X[,2]*X[,5]*X[,7]
    
    y_scale = scale(y, center = TRUE, scale = TRUE) + w
    return(y_scale)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}


pi3 = function(X, w){
  if (ncol(X) >= 6){
    
    y = X[,1] + X[,2]^2 + X[,3]^3 + X[,5]^3 - 2*X[,6]^4 + 0.5
    
    y_scale = scale(y, center = TRUE, scale = TRUE) + w
    return(y_scale)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

mu1 = function(X){
  if (ncol(X) >= 4){
    
    y = 2*X[,1] + 3*X[,2] + X[,1]*X[,4] + X[,3]^2 -4
    
    y_scale = scale(y, center = TRUE, scale = TRUE)
    return(y_scale)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

mu2 = function(X){
  if (ncol(X) >= 6){
    
    y = X[,2]*X[,4]*X[,6] + 2*X[,2]*X[,4]*(1 - X[,6]) + 3*X[,2]*(1 - X[,4])*X[,6] +
      4*X[,2]*(1 - X[,4])*(1 - X[,6]) + 5*(1 - X[,2])*X[,4]*X[,6] + 6*(1 - X[,2])*X[,4]*(1 - X[,6]) +
      7*(1 - X[,2])*(1 - X[,4])*X[,6] + 8*(1 - X[,2])*(1 - X[,4])*(1 - X[,6])
    
    y_scale = scale(y, center = TRUE, scale = TRUE)
    return(y_scale)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

mu3 = function(X){
  if (ncol(X) >= 9){
    I_1 <- 1*(X[,1]>1)
    I_5 <- 1*(X[,5]>1)
    
    y = 4*I_1 + 4*I_5 + 2*X[,7]^2 + X[,9]^3
    
    y_scale = scale(y, center = TRUE, scale = TRUE)
    return(y_scale)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

tau1 = function(X){
  if (ncol(X) >= 9){
    
    y = X[,1] + X[,2] + X[,4]*X[,7] + X[,9]^2 - 2 
    
    y_scale = scale(y, center = TRUE, scale = TRUE) + 0.5
    return(y_scale)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

tau2 = function(X){
  if (ncol(X) >= 9){
    
    y = X[,2] + X[,5]*X[,7] + X[,3]*X[,8]*X[,9] + X[,4]*X[,5]*X[,6]*X[,7] 
    
    y_scale = scale(y, center = TRUE, scale = TRUE) + 0.5
    return(y_scale)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

tau3 = function(X){
  if (ncol(X) >= 9){
    
    y = 0.5*(X[,1]^2 + X[,2] + X[,3]^2 + X[,4] + X[,5]^3 + X[,6] + X[,7]^3 + X[,8] + X[,9]^4 - 5) 
    
    y_scale = scale(y, center = TRUE, scale = TRUE) + 0.5
    return(y_scale)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}


# Powers et al functions --------------------------------------------------
# scenarios
f_2 = function(X){
  
  I_1 <- 1*(X[,1]>1)
  
  if (ncol(X) >= 1){
    
    y = 5*I_1 - 5
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

f_3 = function(X){

  if (ncol(X) >= 1){
    
    y = 2*X[,1] - 4
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

f_5 = function(X){
  
  if (ncol(X) >= 9){
    
    y = X[,1] + X[,3] + X[,5] + X[,7] + X[,8] + X[,9] - 2
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

f_6 = function(X){
  
  I_1 <- 1*(X[,1]>1)
  I_3 <- 1*(X[,3]>1)
  I_5 <- 1*(X[,5]>1)
  I_7 <- 1*(X[,7]>1)
  
  if (ncol(X) >= 9){
    
    y = 4*(I_1*I_3) + 4*(I_5*I_7) + 2*X[,8]*X[,9]
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}

f_3 = function(X){
  I_1 <- 1*(X[,1]>1)
  I_5 <- 1*(X[,5]>1)
  
  
  if (ncol(X) >= 9){
    
    y = 0.5*(X[,1]^2 + X[,2] + X[,3]^2 + X[,4] + X[,5]^3 + X[,6] + X[,7]^3 + X[,8] + X[,9]^4 - 5) 
    
    return(y)
  } else {
    print("Insufficient data for this DGP, increase p")
  }
}



