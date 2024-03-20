#  Plotting of the estimated penetrance and the data-generating curve for simulation studies
plot_penetrance_sim <- function(data, prob, max_age, sex, data_gen_alpha, data_gen_beta, data_gen_threshold, data_gen_asymptote_male, data_gen_asymptote_female) {
  if (prob <= 0 || prob >= 1) {
    stop("prob must be between 0 and 1")
  }
  
  x_values <- seq(0, max_age, length.out = max_age + 1)
  asymptotes <- if (sex == "Male") data$asymptote_male_results else data$asymptote_female_results
  
  # Data-generating curve
  data_gen_asymptote <- if (sex == "Male") data_gen_asymptote_male else data_gen_asymptote_female
  data_generating_curve <- pweibull(x_values - data_gen_threshold, shape = data_gen_alpha, scale = data_gen_beta) * data_gen_asymptote
  
  # Calculating mean, lower, and upper CI for the estimated curve
  distributions <- sapply(asymptotes, function(asymptote) pweibull(x_values - data_gen_threshold, shape = data_gen_alpha, scale = data_gen_beta) * asymptote)
  mean_density <- rowMeans(distributions, na.rm = TRUE)
  ci_lower <- apply(distributions, 1, function(x) quantile(x, probs = (1 - prob) / 2, na.rm = TRUE))
  ci_upper <- apply(distributions, 1, function(x) quantile(x, probs = 1 - (1 - prob) / 2, na.rm = TRUE))
  
  # Plot
  plot(x_values, mean_density,
       type = "l", col = if (sex == "Male") "blue" else "red",
       ylim = c(min(ci_lower, data_generating_curve, na.rm = TRUE), max(ci_upper, data_generating_curve, na.rm = TRUE)),
       xlab = "Age", ylab = "Cumulative Penetrance", main = paste("Penetrance Curve with Credible Interval -", sex)
  )
  lines(x_values, data_generating_curve, col = "black", lty = 1)
  lines(x_values, ci_lower, col = "grey", lty = 2)
  lines(x_values, ci_upper, col = "grey", lty = 2)
  polygon(c(x_values, rev(x_values)), c(ci_lower, rev(ci_upper)), col = rgb(0, 0, 1, 0.1), border = NA)
  
  legend("topleft",
         legend = c("Estimated Penetrance", "95% CI", "Data-Generating Curve"),
         col = c(if (sex == "Male") "blue" else "red", "grey", "black"), lty = c(1, 2, 1), cex = 0.8
  )
}


# Calculating the required penetrance data 

calculate_penetrance_data <- function(data, prob, max_age, data_gen_alpha, data_gen_beta, data_gen_threshold, data_gen_asymptote) {
    if (prob <= 0 || prob >= 1) {
        stop("prob must be between 0 and 1")
    }

    params <- calculate_weibull_parameters(
        data$median_results,
        data$first_quartile_results,
        data$threshold_results,
        data$asymptote_results
    )

    alphas <- params$alpha
    betas <- params$beta
    thresholds <- data$threshold_results
    asymptotes <- data$asymptote_results

    x_values <- seq(0, max_age, length.out = max_age + 1)
    distributions <- vector("list", length(alphas))

    for (i in seq_along(alphas)) {
        distributions[[i]] <- pweibull(x_values - thresholds[i], shape = alphas[i], scale = betas[i]) * asymptotes[i]
    }

    distributions_matrix <- do.call(cbind, distributions)
    mean_density <- rowMeans(distributions_matrix, na.rm = TRUE)

    ci_lower <- apply(distributions_matrix, 1, function(x) quantile(x, probs = (1 - prob) / 2, na.rm = TRUE))
    ci_upper <- apply(distributions_matrix, 1, function(x) quantile(x, probs = 1 - (1 - prob) / 2, na.rm = TRUE))

    data_generating_curve <- pweibull(x_values - data_gen_threshold, shape = data_gen_alpha, scale = data_gen_beta) * data_gen_asymptote

    return(list(
        mean_density = mean_density,
        data_generating_curve = data_generating_curve,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        x_values = x_values
    ))
}

# Calculating MSE for evaluation
calculate_mse <- function(estimated_curve, true_curve) {
    mse <- mean((estimated_curve - true_curve)^2)
    return(mse)
}

# Calculate 95% Credible Interval Coverage
calculate_ci_coverage <- function(true_curve, ci_lower, ci_upper) {
    coverage <- mean(true_curve >= ci_lower & true_curve <= ci_upper)
    return(coverage)
}
