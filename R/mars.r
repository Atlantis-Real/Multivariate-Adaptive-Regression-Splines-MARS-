# constructor, validator and helper for class mars.control objects
new_mars.control <- function(control) {
  structure(control, class="mars.control")}

validate_mars.control <- function(control) {
  stopifnot(is.integer(control$Mmax),is.numeric(control$d),
            is.logical(control$trace))
  if(control$Mmax < 2) {
    warning("Mmax must be >= 2; Reset it to 2")
    control$Mmax <- 2}
  if(control$Mmax %% 2 > 0) {
    control$Mmax <- 2*ceiling(control$Mmax/2)
    warning("Mmax should be an even integer. Reset it to ",control$Mmax)}
  control
}

#' Constructor for `mars.control` objects
#'
#' This function constructs a `mars.control` object that specifies
#' parameters used in the model fitting procedure.
#'
#' @param Mmax Maximum number of basis functions. Should be an even integer. Default value is 2.
#' @param d The coefficient in the penalty term of the generalized cross validation measure. Default is 3.
#' @param trace Should we print status information about the fitting? Default is `FALSE`
#'
#' @return a `mars.control` object
#' @export
#'
#' @examples control <- mars.control(Mmax=10)
mars.control <- function(Mmax=2,d=3,trace=FALSE) {
  Mmax <- as.integer(Mmax)
  control <- list(Mmax=Mmax,d=d,trace=trace)
  control <- validate_mars.control(control)
  new_mars.control(control)
}

control <- mars.control()

#-------------------------------------------------------------------------------

mars <- function(formula, data, control) {
  cc <- match.call() # save the call
  mf <- model.frame(formula, data) # returns a data frame with variables needed
  y <- model.response(mf) # returns the response from the model frame
  mt <- attr(mf, "terms") # obtains attributes of the model frame
  x <- model.matrix(mt, mf)[,-1,drop=FALSE] # creates a design matrix and removing the intercept term
  x_names <- colnames(x)
  validate_mars.control(control)
  fwd <- fwd_stepwise(y, x, control)
  bwd <- bwd_stepwise(fwd, control)

  final_model <- lm(y ~ . - 1, data = data.frame(y=y, bwd$B)) # Fit the final model to the selected basis functions

  out <- c(list( # Package output as a list with required components
    call = cc,
    formula = formula,
    y = y,
    B = bwd$B,
    Bfuncs = bwd$Bfuncs,
    x_names = x_names),
    final_model) # Include the final model fit

    class(out) = c("mars", class(final_model)) # Set class as mars inheriting from class of final_model
    out
}

#-------------------------------------------------------------------------------

h <- function(x,s,t) { # hinge function
  # if x>t, s=+1, this return max(0,x-t)
  # if x<t, s=-1, this return max(0,t-x)
  return(pmax(0,s*(x-t)))
}

#-------------------------------------------------------------------------------

LOF <- function(form,data,control) { # obtain the residual sum of squares (RSS), Output: Value of the GCV criterion
  ff <- lm(form,data)
  RSS <- sum(residuals(ff)^2)
  N <- nrow(data)
  M <- length(coef(ff))-1
  Ctilde <- sum(diag(hatvalues(ff))) + control$d*M
  return(RSS * N/(N-Ctilde)^2)
}

init_B <- function(N,Mmax) {
  # Input: N- # of rows; Mmax: # of basis funcs
  # output: a N by (Mmax+1) dataframe
  B <- data.frame( matrix(NA,nrow=N,ncol=(control$Mmax+1)) )
  B[, 1] <- 1 # first column for intercept: B0
  B[, -1] <- 0  # Initialize other columns to 0
  #names(B) <- c("B0",paste0("B",1:Mmax))
  names(B) <- c("B0", paste0("B", 1:(ncol(B)-1)))
  return(B)
}

split_points <- function(xv,Bm) {
  # input: xv: a variable xv to split
  #        Bm: a parent basis func to split
  # output: feasible splitting points
  #
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}

#-------------------------------------------------------------------------------

fwd_stepwise <- function(y,x,control){

  N <- length(y) # sample size
  n <- ncol(x) # number of predictors
  Mmax <- control$Mmax
  B <- init_B(N,Mmax) # Initialize basis matrix
  Bfuncs <- vector("list", length = control$Mmax+1) # Initialize list to store basis functions

  for(i in 1:(Mmax/2)) { # Loop for forward selection
    M <- 2 * i - 1 # Set the value of M from the value of i
    if(control$trace) cat("M",M,"\n")
    lof_best <- Inf
    for(m in 1:M) { # choose a basis function to split
      svars <- setdiff(1:n, Bfuncs[[m]][, "v"]) # store the indices of variables that are not already in basis function m
      if (control$trace) cat("M, m, svars", M, m, svars, "\n")
      for(v in svars){ # select a variable to split on
        tt <- split_points(x[,v],B[,m])
        for(t in tt) { # Loop over split points
          Bnew <- data.frame(B[,1:M], # Generate new basis matrix
                             Btem1 = h(x[, v], +1, t) * B[, m], # Left child basis
                             Btem2 = h(x[, v], -1, t) * B[, m]) # Right child basis
          gdat <- data.frame(y=y,Bnew)
          lof <- LOF(form = y ~ ., data = gdat, control = control) # Compute LOF for new basis matrix
          if(lof < lof_best) { # Update best LOF and split if new LOF is smaller
            lof_best <- lof
            split_best <- c(m=m, v=v, t=t)
          } # end if
        } # end loop over splits
      } # end loop over variables
    } # end loop over basis functions to split
    mstar <- split_best["m"]; vstar <- split_best["v"]; tstar <- split_best["t"]
    cat("[Info] best (m,v,t,lof): (",mstar,vstar,tstar,lof_best,")\n")
    # Update basis matrix and basis functions
    B[, M + 1] <- B[, mstar] * h(x[, vstar], -1, tstar) # Left child basis
    B[, M + 2] <- B[, mstar] * h(x[, vstar], +1, tstar) # Right child basis
    Bfuncs[[M + 1]] = rbind(Bfuncs[[mstar]], c(-1, vstar, tstar)) # Add left child basis function
    Bfuncs[[M + 2]] = rbind(Bfuncs[[mstar]], c(+1, vstar, tstar)) # Update parent basis with the right child basis
    colnames(Bfuncs[[M + 1]]) <- colnames(Bfuncs[[M + 2]]) <- c("s", "v", "t")
  } # end for loop over i
  colnames(B) <- paste0("B",0:(ncol(B)-1))
  return(list(y=y,B=B,Bfuncs=Bfuncs))
}

#-------------------------------------------------------------------------------

bwd_stepwise <- function(fwd, control) {
  #fwd is a list with elements y, B and Bfuncs
  Mmax <- length(fwd$Bfuncs)   # Total number of basis functions
  Jstar <- seq_len(Mmax)
  Kstar <- Jstar
  dat <- data.frame(y = fwd$y, fwd$B[, Jstar]) # Create dataframe with initial selected variables
  lofstar <- LOF(y~.-1,dat,control) # Compute LOF for initial selected variables

  for (M in Mmax:2) { # Loop from Mmax to 2
    b <- Inf
    L <- Kstar
    if(control$trace) cat("L:",L,"\n")

    for (m in L) { # Inner loop for variable removal
      K <- setdiff(L, m)
      dat <- data.frame(y = fwd$y, fwd$B[, K]) # Create dataframe with current set of variables
      lof <- LOF(y~.,dat,control) # Compute LOF for current set of variables
      if(control$trace) cat("M:K:lof",M,":",K,":",lof,"\n")
      if (lof < b) { # Update Kstar if current LOF is smaller
        b <- lof
        Kstar <- K
      }

      if (lof < lofstar) { # Update Jstar and lofstar if current LOF is smaller
        lofstar <- lof
        Jstar <- K
      }
    }
    if(control$trace) cat("M:Jstar:lofstar",M,":",Jstar,":",lofstar,"\n")
  }
  Jstar <- c(1,Jstar) # Ensure the intercept term is included
  return(list(y=fwd$y,B=fwd$B[,Jstar],Bfuncs=fwd$Bfuncs[Jstar]))
}

