#' Plot method for mars objects
#'
#' This function plots the basis functions of a mars model object.
#'
#' @param object A mars object.
#' @param ... Additional arguments to be passed to methods.
#' @return Plot of basis functions.
#' @export
plot.mars <- function(object, ...) {
  n_basis <- length(object$Bfuncs)
  coeffs <- coef(object) # Extract coefficients

  for (i in 1:n_basis) {
    basis <- object$Bfuncs[[i]]
    basis_name <- names(coeffs)[i]  # Use coefficient names as basis names

    if (is.null(basis)) {
      cat(basis_name, ":", sep = "", "\n")
      cat("Intercept detected. Skipping...\n\n")
      next  # Skip to the next iteration
    }
    n_components <- nrow(basis)

    cat(basis_name, ":", sep = "", "\n")

    associated_predictor_indices <- basis[, "v"]
    associated_predictors <- object$x_names[associated_predictor_indices]

    if (n_components == 1) {
      cat("1 component detected.\n")
      plot_1d_basis(basis, associated_predictors)
    } else if (n_components == 2) {
      cat("2 components detected.\n")
      plot_2d_basis(basis, associated_predictors)
    } else if (n_components >= 3) {
      cat("3+ components detected. Unable to plot beyond 3D.\n")
    }

    cat("\n")
  }
}

plot_1d_basis <- function(basis, associated_predictors) {
  if (nrow(basis) > 0) {
    v <- basis[1, "v"]
    t <- basis[1, "t"]
    s <- basis[1, "s"]
    plot_hinge_function(v, t, s, associated_predictors)
  } else {
    cat("No components detected for 1D plot.\n")
  }
}

plot_2d_basis <- function(basis, associated_predictors) {
  if (nrow(basis) > 1) {
    v1 <- basis[1, "v"]
    t1 <- basis[1, "t"]
    s1 <- basis[1, "s"]

    v2 <- basis[2, "v"]
    t2 <- basis[2, "t"]
    s2 <- basis[2, "s"]

    plot_hinge_function_2d(v1, t1, s1, v2, t2, s2, associated_predictors)
  } else {
    cat("Insufficient components detected for 2D plot.\n")
  }
}

plot_hinge_function <- function(v, t, s, associated_predictors) {
  curve(h(x, s, t), -2, 2, main = paste(associated_predictors))
}

plot_hinge_function_2d <- function(v1, t1, s1, v2, t2, s2, associated_predictors) {
  x1 <- seq(-2, 2, length.out = 100)
  x2 <- seq(-2, 2, length.out = 100)
  z <- outer(x1, x2, Vectorize(function(x1, x2) {
    h(x1, s1, t1) * h(x2, s2, t2)
  }))
  persp(x1, x2, z, theta = 30, phi = 30, col = "lightblue",
        main = paste(paste(associated_predictors, collapse = ":")),
        xlab = associated_predictors[1], ylab = associated_predictors[2], zlab = "z")
}

