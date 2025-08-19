source("R/mars.r")

anova.mars <- function(object, ...) {
  # Check if the model is of class 'mars'
  if (!inherits(object, "mars")) {
    stop("Input must be a 'mars' object")
  }

  cat("Analysis of Variance Table\n\n")

  # Extract basis functions and their corresponding predictor variables
  Bfuncs <- object$Bfuncs

  # Extract coefficients (Bs)
  coef_values <- object$B

  # Initialize a vector to store the degrees of freedom for each basis function
  df <- rep(1, length(Bfuncs))

  # Obtain residual degrees of freedom
  residual_df <- object$df.residual

  # Compute sum of squares for each coefficient
  sum_sq <- sapply(coef_values, function(coef) sum((coef - mean(coef, na.rm = TRUE))^2, na.rm = TRUE))

  # Compute mean squares
  mean_sq <- sum_sq / df

  # Compute the mean square of residuals
  mean_sq_residual <- sum_sq[length(sum_sq)] / residual_df

  # Compute F-values
  f_value <- mean_sq / mean_sq_residual

  # Compute p-values (Pr(>F))
  p_value <- pf(f_value, df, residual_df, lower.tail = FALSE)

  # Calculate sum of squares for residuals
  residual_sum_sq <- sum((object$fitted.values - object$y)^2)

  # Calculate mean square for residuals
  residual_mean_sq <- residual_sum_sq / residual_df

  # Output the ANOVA results using cat
  cat("Response:\n")
  cat("Term \t Df \t Sum Sq \t Mean Sq \t F value \t Pr(>F) \n")
  for (i in seq_along(df)) {
    cat(names(coef_values)[i], "\t", df[i], "\t", sum_sq[i], "\t", mean_sq[i], "\t", f_value[i], "\t", p_value[i], "\n")
  }

  # Output residuals row
  cat("Residual", residual_df, "\t", residual_sum_sq, "\t", residual_mean_sq, "\n")
  cat("---")
  cat("\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n", "\n")
}
