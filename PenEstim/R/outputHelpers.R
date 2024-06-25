#' Combine Chains
#' Function to combine the posterior samples from the multiple chains.
#'
#' @param results A list of MCMC chain results.
#'
#' @return A list with combined results, including median, threshold, first quartile, and asymptote values.
#'
combine_chains <- function(results) {
  list(
    median_male_results = do.call(c, lapply(results, function(x) x$median_male_samples)),
    median_female_results = do.call(c, lapply(results, function(x) x$median_female_samples)),
    threshold_male_results = do.call(c, lapply(results, function(x) x$threshold_male_samples)),
    threshold_female_results = do.call(c, lapply(results, function(x) x$threshold_female_samples)),
    first_quartile_male_results = do.call(c, lapply(results, function(x) x$first_quartile_male_samples)),
    first_quartile_female_results = do.call(c, lapply(results, function(x) x$first_quartile_female_samples)),
    asymptote_male_results = do.call(c, lapply(results, function(x) x$asymptote_male_samples)),
    asymptote_female_results = do.call(c, lapply(results, function(x) x$asymptote_female_samples)),
    loglikelihood_current_results = do.call(c, lapply(results, function(x) x$loglikelihood_current)),
    loglikelihood_proposal_results = do.call(c, lapply(results, function(x) x$loglikelihood_proposal)),
    acceptance_ratio_results = do.call(c, lapply(results, function(x) x$acceptance_ratio)),
    median_male_proposals = do.call(c, lapply(results, function(x) x$median_male_proposals)),
    median_female_proposals = do.call(c, lapply(results, function(x) x$median_female_proposals)),
    threshold_male_proposals = do.call(c, lapply(results, function(x) x$threshold_male_proposals)),
    threshold_female_proposals = do.call(c, lapply(results, function(x) x$threshold_female_proposals)),
    first_quartile_male_proposals = do.call(c, lapply(results, function(x) x$first_quartile_male_proposals)),
    first_quartile_female_proposals = do.call(c, lapply(results, function(x) x$first_quartile_female_proposals)),
    asymptote_male_proposals = do.call(c, lapply(results, function(x) x$asymptote_male_proposals)),
    asymptote_female_proposals = do.call(c, lapply(results, function(x) x$asymptote_female_proposals))
  )
}

#' Generate Summary
#' @description Function to generate summary statistics
#'
#' @param data A list with combined results.
#'
#' @return A summary data frame containing Median, threshold, First Quartile, and Asymptote Value.
#'
generate_summary <- function(data) {
  summary_data <- data.frame(
    Median_Male = data$median_male_results,
    Median_Female = data$median_female_results,
    Threshold_Male = data$threshold_male_results,
    Threshold_Female = data$threshold_female_results,
    First_Quartile_Male = data$first_quartile_male_results,
    First_Quartile_Female = data$first_quartile_female_results,
    Asymptote_Male = data$asymptote_male_results,
    Asymptote_Female = data$asymptote_female_results
  )
  summary(summary_data)
}

#' Generate Posterior Density Plots
#' 
#' Generates histograms of the posterior samples for the different parameters
#'
#' @param data A list with combined results.
#'
generate_density_plots <- function(data) {
  # Set the plotting parameters
  par(mfrow = c(3, 2), las = 1, mar = c(5, 4, 4, 2) + 0.1) # Adjust grid for 5 plots, now 3 rows and 2 columns

  # Define the specific vectors to plot
  plot_names <- c("median_male_results", "first_quartile_male_results", 
                  "asymptote_male_results", "threshold_male_results", 
                  "median_female_results", "first_quartile_female_results",
                  "asymptote_female_results", "threshold_female_results")

  for (name in plot_names) {
    if (is.null(data[[name]]) || length(data[[name]]) == 0) {
      next # Skip this iteration if the data is empty
    }

    mod_name <- gsub("_", " ", name)
    mod_name <- paste0(toupper(substring(mod_name, 1, 1)), substring(mod_name, 2))

    # Set xlim based on the name of the vector
    xlim <- if (name %in% c("median_male_results", "first_quartile_male_results", 
                            "threshold_male_results","median_female_results", 
                            "first_quartile_female_results", 
                            "threshold_female_results")) {
      c(0, 100)
    } else if (name %in% c("asymptote_male_results", "asymptote_female_results")) {
      c(0, 1) # Assuming asymptote values are between 0 and 1
    } else {
      range(data[[name]], na.rm = TRUE) # Default to data range
    }

    # Ensure xlim is finite and valid
    if (any(is.infinite(xlim))) {
      xlim <- c(min(data[[name]], na.rm = TRUE), max(data[[name]], na.rm = TRUE))
    }

    # Create the histogram
    hist(data[[name]],
      main = paste("Density Plot of", mod_name),
      xlab = mod_name,
      freq = FALSE,
      xlim = xlim,
      breaks = "Sturges", # Default break algorithm
      col = "lightblue" # Optional: color for the histogram
    )
  }

  # Reset to default single plot setting after plotting
  par(mfrow = c(1, 1))
}


#' Plot Trace
#' @param results A list of MCMC chain results.
#' @param n_chains The number of chains.
#'
plot_trace <- function(results, n_chains) {
  par(mfrow = c(n_chains, 4)) # Set up a grid for the plots
  for (chain_id in 1:n_chains) {
    # Extract results for the current chain
    median_male_results <- results[[chain_id]]$median_male_samples
    median_female_results <- results[[chain_id]]$median_female_samples
    threshold_male_results <- results[[chain_id]]$threshold_male_samples
    threshold_female_results <- results[[chain_id]]$threshold_female_samples
    first_quartile_male_results <- results[[chain_id]]$first_quartile_male_samples
    first_quartile_female_results <- results[[chain_id]]$first_quartile_female_samples
    asymptote_male_results <- results[[chain_id]]$asymptote_male_samples
    asymptote_female_results <- results[[chain_id]]$asymptote_female_samples

    # Create trace plots for the current chain

    plot(median_male_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Median - Male"), xlab = "Iteration", ylab = "Median")
    plot(median_female_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Median - Female"), xlab = "Iteration", ylab = "Median")
    plot(threshold_male_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Threshold - Male"), xlab = "Iteration", ylab = "Threshold")
    plot(threshold_female_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Threshold - Female"), xlab = "Iteration", ylab = "Threshold")
    plot(first_quartile_male_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of First Quartile - Male"), xlab = "Iteration", ylab = "First Quartile")
    plot(first_quartile_female_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of First Quartile - Female"), xlab = "Iteration", ylab = "First Quartile")
    plot(asymptote_male_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Asymptote Male"), xlab = "Iteration", ylab = "Asymptote")
    plot(asymptote_female_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Asymptote Female"), xlab = "Iteration", ylab = "Asymptote")
  }
}

# Trace for just a single chain
plot_traceSingle <- function(results) {
  par(mfrow = c(2, 2)) # Set up a grid for the plots
  # Extract results for the current chain
  median_male_results <- results$median_male_samples
  median_female_results <- results$median_female_samples
  threshold_male_results <- results$threshold_male_samples
  threshold_female_results <- results$threshold_female_samples
  first_quartile_male_results <- results$first_quartile_male_samples
  first_quartile_female_results <- results$first_quartile_female_samples
  asymptote_male_results <- results$asymptote_male_samples
  asymptote_female_results <- results$asymptote_female_samples

  # Create trace plots for the current chain

  plot(median_male_results, type = "l", main = "Trace plot of Median - Male", xlab = "Iteration", ylab = "Median")
  plot(median_female_results, type = "l", main = "Trace plot of Median - Female", xlab = "Iteration", ylab = "Median")
  plot(threshold_male_results, type = "l", main = "Trace plot of Threshold - Male", xlab = "Iteration", ylab = "Threshold")
  plot(threshold_female_results, type = "l", main = "Trace plot of Threshold - Female", xlab = "Iteration", ylab = "Threshold")
  plot(first_quartile_male_results, type = "l", main = "Trace plot of First Quartile - Male", xlab = "Iteration", ylab = "First Quartile")
  plot(first_quartile_female_results, type = "l", main = "Trace plot of First Quartile - Female", xlab = "Iteration", ylab = "First Quartile")
  plot(asymptote_male_results, type = "l", main = "Trace plot of Asymptote - Male", xlab = "Iteration", ylab = "Asymptote")
  plot(asymptote_female_results, type = "l", main = "Trace plot of Asymptote - Female", xlab = "Iteration", ylab = "Asymptote")
}

# Running mean calculation
running_mean <- function(res) {
  n <- length(res)
  running_mean <- numeric(n)
  for (i in 1:n) {
    running_mean[i] <- mean(res[1:i])
  }
  return(running_mean)
}

# Running var calculation
running_variance <- function(res) {
  n <- length(res)
  running_var <- numeric(n)
  for (i in 1:n) {
    running_var[i] <- var(res[1:i])
  }
  return(running_var)
}

#' Print Rejection Rates
#' @param results A list of MCMC chain results.
#'
printRejectionRates <- function(results) {
  rejection_rates <- sapply(results, function(x) x$rejection_rate)
  cat("Rejection rates: ", rejection_rates, "\n")
}

#' Apply Burn-In
#'
#' @param results A list of MCMC chain results.
#' @param burn_in The fraction roportion of results to discard as burn-in (0 to 1). The default is no burn-in, burn_in=0.
#'
#' @return A list of results with burn-in applied.
#'
apply_burn_in <- function(results, burn_in) {
  # Ensure 'results' is a list and has at least one chain
  if (!is.list(results) || length(results) < 1) {
    stop("results must be a list with at least one chain.")
  }

  # Ensure 'burn_in' is numeric and between 0 and 1
  if (!is.numeric(burn_in) || burn_in <= 0 || burn_in >= 1) {
    stop("burn_in must be a numeric value between 0 and 1.")
  }

  # Function to perform burn-in on a single chain (list of numeric vectors)
  burn_in_chain <- function(chain, burn_in) {
    lapply(chain, function(param_results) {
      n_results <- length(param_results)
      burn_in_count <- round(n_results * burn_in)
      if (burn_in_count >= n_results) {
        stop("burn_in_count must be less than the total number of results (n_results).")
      }
      param_results[(burn_in_count + 1):n_results]
    })
  }

  # Apply burn-in to all results
  lapply(results, function(chain) {
    burn_in_chain(chain, burn_in)
  })
}

#' Apply Thinning
#'
#' @param results A list of MCMC chain results.
#' @param thinning_factor The factor by which to thin the results (positive integer). The default thinning factor is 1, which implies no thinning.
#'
#' @return A list of results with thinning applied.
#'
apply_thinning <- function(results, thinning_factor) {
  # Ensure 'results' is a list and has at least one chain
  if (!is.list(results) || length(results) < 1) {
    stop("results must be a list with at least one chain.")
  }

  # Ensure 'thinning_factor' is a positive integer
  if (!is.numeric(thinning_factor) || thinning_factor <= 0 || !is.integer(thinning_factor)) {
    stop("thinning_factor must be a positive integer.")
  }

  # Function to perform thinning on a single chain (list of numeric vectors)
  thin_chain <- function(chain, thinning_factor) {
    lapply(chain, function(param_results) {
      if (length(param_results) < thinning_factor) {
        stop("Thinning factor is larger than the number of results.")
      }
      param_results[seq(1, length(param_results), by = thinning_factor)]
    })
  }

  # Apply thinning to all chains
  lapply(results, function(chain) {
    thin_chain(chain, thinning_factor)
  })
}

#' Plot Weibull Distribution with Credible Intervals
#'
#' This function plots the Weibull distribution with credible intervals for the given data. 
#' It allows for visualization of penetrance curves for individuals based on their genetic 
#' and demographic information.
#'
#' @param data Data frame, containing individual demographic and genetic information. Must include columns for 'sex', 'age', 'aff' (affection status), and 'geno' (genotype).
#' @param prob Numeric, the probability level for the credible intervals. Must be between 0 and 1.
#' @param max_age Integer, the maximum age considered in the analysis.
#' @param sex Character, specifying the sex of the individuals for the plot ("Male", "Female", or "NA" for not applicable). Default is "NA".
#'
#' @return A plot showing the Weibull distribution with credible intervals.
#'
plot_penetrance <- function(data, prob, max_age, sex = "NA") {
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

  x_values <- seq(0, max_age, length.out = max_age + 1)

  # Depending on the sex, select the corresponding asymptote values or prepare for both
  asymptotes_male <- data$asymptote_male_results
  asymptotes_female <- data$asymptote_female_results

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
        xlab = "Age", ylab = "Cumulative Penetrance", main = "Penetrance Curve with Credible Interval"
      )
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
  } else if (sex == "Female") {
    plot_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red")
    legend_text <- "Female Estimated Penetrance"
  } else {
    plot_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue")
    plot_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red", add = TRUE)
    legend_text <- c("Male Estimated Penetrance", "Female Estimated Penetrance")
  }

  legend("topleft",
    legend = legend_text,
    col = c("blue", "red"),
    lty = c(1, 1),
    cex = 0.8
  )
}