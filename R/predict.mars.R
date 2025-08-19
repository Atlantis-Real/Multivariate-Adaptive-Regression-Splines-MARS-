load("TestFiles/testpredict.RData")

predict.mars <- function(object,newdata) {
  if(missing(newdata) || is.null(newdata)) {
    B <- as.matrix(object$B)
  }
  else {
    tt <- terms(object$formula,data=newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1] # remove intercept
    B <- make_B(X,object$Bfuncs)
  }
  beta <- object$coefficients
  drop(B %*% beta)
}

make_B <- function(X, Bfuncs){
  B <- matrix(1, nrow = nrow(X), ncol = length(Bfuncs)) # Initialize B matrix
  for (m in 2:length(Bfuncs)){ # Loop through each non-constant basis functions
    for (i in 1:nrow(Bfuncs[[m]])){ # Loop through each branch (hinge function)
      Bas <- Bfuncs[[m]][i,]
      s <- Bas["s"]; v <- Bas["v"]; t <- Bas["t"]
      B[,m] <- B[,m] * h(s, X[,v], t)  # Update B[,m] using hinge function and new data on covariates.
    }
  }
  return(B)
}


make_B <- function(X, Bfuncs){
  B <- matrix(1, nrow = nrow(X), ncol = length(Bfuncs)) # Initialize B matrix
  for (m in 2:length(Bfuncs)){ # Loop through each non-constant basis functions
    for (i in 1:nrow(Bfuncs[[m]])){ # Loop through each branch (hinge function)
      Bas <- Bfuncs[[m]][i,]
      s <- Bas["s"]
      v <- Bas["v"]
      t <- Bas["t"]
      B[,m] <- B[,m]*h(X[, v], s, t)  # Update B[,m] using hinge function and new data on covariates.
    }
  }
  return(B)
}



h <- function(x,s,t) {
  return(pmax(0,s*(x-t)))
  # if x>t, s=+1, this return max(0,x-t)
  # if x<t, s=-1, this return max(0,t-x)
}
