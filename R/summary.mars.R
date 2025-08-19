#source("R/mars.r")

#' Summary method for mars objects
#'
#' This function provides a summary of the mars model object.
#'
#' @param object A mars object.
#' @param ... Additional arguments to be passed to methods.
#' @return A summary of the mars model object.
#' @export

summary.mars <- function(object, ...) {
  coef_values <- object$coefficients
  coef_names <- names(coef_values)

  cat("Basis functions:\n")
  # Extract predictor names and coefficients
  predictor_names <- object$x_names
  coefficients <- object$coefficients

  predictor_associations <- vector("list", length(coefficients)-1)# Initialize a list to store predictors associated with each coefficient
  for (i in seq_along(coef_values)) { # Iterate over coefficients
    basis_functions <- object$Bfuncs[[i]] # Extract the matrix of predictors associated with this coefficient
    associated_predictor_indices <- basis_functions[, "v"] # Extract predictor indices from the matrix
    associated_predictors <- predictor_names[associated_predictor_indices] # Retrieve the corresponding predictor names
    s_values <- basis_functions[, "s"] # Extract s values
    t_values <- basis_functions[, "t"] # Extract t values
    predictor_associations[[i]] <- associated_predictors # Store predictors associated with this coefficient
    cat(coef_names[i], ":",sep="", "\n") # Print coefficient name

    if (is.null(basis_functions)) {
      cat(" Intercept\n")
    } else if (length(associated_predictors) == 1) {
      cat(" Component 1:", " variable ", associated_predictors,";"," sign ", s_values,";",sep=""," knot at ", t_values, "\n")
      # If only one predictor associated, print it as Component 1
    } else {
      for (j in seq_along(associated_predictors)) {
        cat(" Component", j, ":"," variable ", associated_predictors[j],";"," sign ", s_values[j],";",sep=""," knot at ", t_values[j], "\n")
        # If multiple predictors, print each as Component i
      }
    }
  }

  cat("\n")

  cat("Call:\n")
  print(object$call)

  cat("\nResiduals:\n")
  res_summary <- summary(object$residuals)
  cat("Min:\t","\t","1Q:\t","\t","Median:\t","3Q:\t","\t","Max:\t","\n")
  cat(res_summary[1],"\t",res_summary[2],"\t",res_summary[3],"\t",res_summary[4],"\t",res_summary[5], "\n")

  cat("\nCoefficients:\n")
  std_errors <- sqrt(diag(vcov(object)))
  t_values <- coef_values / std_errors
  p_values <- 2 * pt(abs(t_values), df = object$df.residual, lower.tail = FALSE)

  cat("\t Estimate\t Std. Error\t t value\t Pr(>|t|)\n")
  for (i in seq_along(coef_values)) {
    cat(coef_names[i],"\t",coef_values[i],"\t",std_errors[i],"\t",t_values[i],"\t", p_values[i],"\n")
  }
  cat("---")
  cat("\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n", "\n")

  # Calculate and print residual standard error
  res_std_error <- sd(residuals(object))
  cat("Residual standard error:", res_std_error, "on", object$df.residual, "degrees of freedom\n")

  # Calculate and print Multiple R-squared
  R_squared <- 1 - sum(residuals(object)^2) / sum((object$y - mean(object$y))^2)
  cat("Multiple R-squared:", R_squared, "\n")

  # Calculate and print Adjusted R-squared
  n <- length(object$y)
  p <- length(coef_values)
  adj_R_squared <- 1 - ((1 - R_squared) * (n - 1) / (n - p - 1))
  cat("Adjusted R-squared:", adj_R_squared, "\n")

  # Calculate and print F-statistic
  mean_sq_model <- sum(coef_values^2) / length(coef_values)
  mean_sq_residual <- sum(residuals(object)^2) / object$df.residual
  F_statistic <- mean_sq_model / mean_sq_residual
  df_model <- length(coef_values)
  df_residual <- object$df.residual
  cat("F-statistic:", F_statistic, "on", df_model, "and", df_residual, "DF\n")

  # Calculate and print p-value
  p_value <- pf(F_statistic, df_model, df_residual, lower.tail = FALSE)
  cat("p-value:", p_value, "\n")
}
