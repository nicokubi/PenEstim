#  Plotting of the estimated penetrance and the data-generating curve for simulation studies
plot_penetrance_sim <- function(data, prob, max_age, sex, data_gen_alpha, data_gen_beta, data_gen_threshold, data_gen_asymptote) {
  if (prob <= 0 || prob >= 1) {
    stop("prob must be between 0 and 1")
  }
  
  params_male <- calculate_weibull_parameters(
    data$median_male_results,
    data$first_quartile_male_results,
    data$threshold_male_results
  )
  
  params_female <- calculate_weibull_parameters(
    data$median_female_results,
    data$first_quartile_female_results,
    data$threshold_female_results
  )
  
  alphas_male <- params_male$alpha
  betas_male <- params_male$beta
  thresholds_male <- data$threshold_male_results
  alphas_female <- params_female$alpha
  betas_female <- params_female$beta
  thresholds_female <- data$threshold_female_results
  
  # Depending on the sex, select the corresponding asymptote values or prepare for both
  asymptotes_male <- data$asymptote_male_results
  asymptotes_female <- data$asymptote_female_results
  
  # x-values
  x_values <- seq(0, max_age, length.out = max_age + 1)
  
  # Data-generating curve
  data_generating_curve <- pweibull(x_values - data_gen_threshold, shape = data_gen_alpha, scale = data_gen_beta) * data_gen_asymptote
  
  plot_distribution <- function(alphas, betas, thresholds, asymptotes, x_values, prob, color, add = FALSE) {
    distributions <- mapply(function(alpha, beta, threshold, asymptote) {
      pweibull(x_values - threshold, shape = alpha, scale = beta) * asymptote
    }, alphas, betas, thresholds, asymptotes, SIMPLIFY = FALSE)
    
    distributions_matrix <- matrix(unlist(distributions), nrow = length(x_values), byrow = FALSE)
    mean_density <- rowMeans(distributions_matrix, na.rm = TRUE)
    ci_lower <- apply(distributions_matrix, 1, quantile, probs = (1 - prob) / 2, na.rm = TRUE)
    ci_upper <- apply(distributions_matrix, 1, quantile, probs = 1 - (1 - prob) / 2, na.rm = TRUE)
    
    if (!add) {
      plot(x_values, mean_density,
           type = "l", col = color,
           ylim = c(min(ci_lower, na.rm = TRUE), max(ci_upper, na.rm = TRUE)),
           xlab = "Age", ylab = "Cumulative Penetrance", main = "Penetrance Curve with Credible Interval")
      lines(x_values, data_generating_curve, col = "black", lty = 1)
      
    } else {
      lines(x_values, mean_density, col = color)
      
    }
    lines(x_values, ci_lower, col = color, lty = 2)
    lines(x_values, ci_upper, col = color, lty = 2)
    polygon(c(x_values, rev(x_values)), c(ci_lower, rev(ci_upper)), col = rgb(0, 0, 1, 0.1), border = NA)
  }
  
  if (sex == "Male") {
    plot_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue")
    legend_text <- "Male Estimated Penetrance"
    legend("topleft",
           legend = legend_text,
           col = c("blue"),
           lty = c(1),
           cex = 0.8)
  } else if (sex == "Female") {
    plot_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red")
    legend_text <- "Female Estimated Penetrance"
    legend("topleft",
           legend = legend_text,
           col = c("red"),
           lty = c(1),
           cex = 0.8)
  } else {
    plot_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue")
    plot_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red", add = TRUE)
    legend_text <- c("Male Estimated Penetrance", "Female Estimated Penetrance")
    legend("topleft",
           legend = legend_text,
           col = c("red"),
           lty = c(1),
           cex = 0.8)
  }
  
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

plot_penetrance_sim_female <- function(data, prob, max_age, data_gen_alpha, data_gen_beta, data_gen_threshold, data_gen_asymptote) {
  if (prob <= 0 || prob >= 1) {
    stop("prob must be between 0 and 1")
  }
  
  params_female <- calculate_weibull_parameters(
    data$median_female_results,
    data$first_quartile_female_results,
    data$threshold_female_results
  )
  
  alphas_female <- params_female$alpha
  betas_female <- params_female$beta
  thresholds_female <- data$threshold_female_results
  asymptotes_female <- data$asymptote_female_results
  
  x_values <- seq(0, max_age, length.out = max_age + 1)
  data_generating_curve <- pweibull(x_values - data_gen_threshold, shape = data_gen_alpha, scale = data_gen_beta) * data_gen_asymptote
  
  plot_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red")
  legend_text <- "Female Estimated Penetrance"
  legend("topleft",
         legend = legend_text,
         col = c("red"),
         lty = c(1),
         cex = 0.8)
}

# Plot for females only
plot_penetrance_sim_female <- function(data, prob, max_age, data_gen_alpha, data_gen_beta, data_gen_threshold, data_gen_asymptote) {
  if (prob <= 0 || prob >= 1) {
    stop("prob must be between 0 and 1")
  }
  
  params_female <- calculate_weibull_parameters(
    data$median_female_results,
    data$first_quartile_female_results,
    data$threshold_female_results
  )
  
  alphas_female <- params_female$alpha
  betas_female <- params_female$beta
  thresholds_female <- data$threshold_female_results
  asymptotes_female <- data$asymptote_female_results
  
  x_values <- seq(0, max_age, length.out = max_age + 1)
  data_generating_curve <- pweibull(x_values - data_gen_threshold, shape = data_gen_alpha, scale = data_gen_beta) * data_gen_asymptote
  
  plot_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red")
  legend_text <- "Female Estimated Penetrance"
  legend("topleft",
         legend = legend_text,
         col = c("red"),
         lty = c(1),
         cex = 0.8)
}

# Plot for males only
plot_penetrance_sim_male <- function(data, prob, max_age, data_gen_alpha, data_gen_beta, data_gen_threshold, data_gen_asymptote) {
  if (prob <= 0 || prob >= 1) {
    stop("prob must be between 0 and 1")
  }
  
  params_male <- calculate_weibull_parameters(
    data$median_male_results,
    data$first_quartile_male_results,
    data$threshold_male_results
  )
  
  alphas_male <- params_male$alpha
  betas_male <- params_male$beta
  thresholds_male <- data$threshold_male_results
  asymptotes_male <- data$asymptote_male_results
  
  x_values <- seq(0, max_age, length.out = max_age + 1)
  data_generating_curve <- pweibull(x_values - data_gen_threshold, shape = data_gen_alpha, scale = data_gen_beta) * data_gen_asymptote
  
  plot_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue")
  legend_text <- "Male Estimated Penetrance"
  legend("topleft",
         legend = legend_text,
         col = c("blue"),
         lty = c(1),
         cex = 0.8)
}
