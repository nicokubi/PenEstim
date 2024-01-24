#' Combine Chains
#' Function to combine the posterior samples from the multiple chains.
#'
#' @param results A list of MCMC chain results.
#'
#' @return A list with combined results, including median, shift, first quartile, and asymptote values.
#' @examples
#' combine_results <- combine_chains(mcmc_results)
#'
combine_chains <- function(results) {
  list(
    median_results = do.call(c, lapply(results, function(x) x$median_samples)),
    shift_results = do.call(c, lapply(results, function(x) x$shift_samples)),
    first_quartile_results = do.call(c, lapply(results, function(x) x$first_quartile_samples)),
    asymptote_results = do.call(c, lapply(results, function(x) x$asymptote_samples)),
    loglikelihood_current_results = do.call(c, lapply(results, function(x) x$loglikelihood_current)),
    loglikelihood_proposal_results = do.call(c, lapply(results, function(x) x$loglikelihood_proposal)),
    acceptance_ratio_results = do.call(c, lapply(results, function(x) x$acceptance_ratio))
  )
}

#' Generate Summary
#' @description Function to generate summary statistics
#'
#' @param data A list with combined results.
#'
#' @return A summary data frame containing Median, Shift, First Quartile, and Asymptote Value.
#'
#' @examples
#' summary_stats <- generate_summary(combine_results)
#'
generate_summary <- function(data) {
  summary_data <- data.frame(
    Median = data$median_results,
    Shift = data$shift_results,
    First_Quartile = data$first_quartile_results,
    Asymptote_Value = data$asymptote_results
  )
  summary(summary_data)
}

#' Generate Density Plots
#'
#' @param data A list with combined results.
#'
#' @examples
#' generate_density_plots(combine_results)
#'
generate_density_plots <- function(data) {
  # Set the plotting parameters
  par(mfrow = c(2, 2), las = 1, mar = c(5, 4, 4, 2) + 0.1)

  # Define the specific vectors to plot
  plot_names <- c("median_results", "first_quartile_results", "asymptote_results", "shift_results")

  for (name in plot_names) {
    if (is.null(data[[name]]) || length(data[[name]]) == 0) {
      next # Skip this iteration if the data is empty
    }

    mod_name <- gsub("_", " ", name)
    mod_name <- paste0(toupper(substring(mod_name, 1, 1)), substring(mod_name, 2))

    # Set xlim based on the name of the vector
    xlim <- if (name %in% c("median_results", "first_quartile_results", "shift_results")) {
      c(0, 100)
    } else if (name == "asymptote_results") {
      c(0.15, 1)
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
      xaxp = c(min(xlim), max(xlim), 10)
    )
  }
}

#' Plot Trace
#' @param results A list of MCMC chain results.
#' @param n_chains The number of chains.
#'
#' @examples
#' plot_trace(mcmc_results, n_chains = 4)
#'
plot_trace <- function(results, n_chains) {
  par(mfrow = c(n_chains, 4)) # Set up a grid for the plots
  for (chain_id in 1:n_chains) {
    # Extract results for the current chain
    median_results <- results[chain_id]$median_results
    shift_results <- results[chain_id]$shift_results
    first_quartile_results <- results[chain_id]$first_quartile_results
    asymptote_results <- results[chain_id]$asymptote_results

    # Create trace plots for the current chain

    plot(median_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Median"), xlab = "Iteration", ylab = "Median")
    plot(shift_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Shift"), xlab = "Iteration", ylab = "Shift")
    plot(first_quartile_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of First Quartile"), xlab = "Iteration", ylab = "First Quartile")
    plot(asymptote_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Asymptote"), xlab = "Iteration", ylab = "Asymptote")
  }
}

#' Print Rejection Rates
#' @param results A list of MCMC chain results.
#'
#' @examples
#' printRejectionRates(mcmc_results)
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
#' @examples
#' burned_results <- apply_burn_in(mcmc_results, burn_in = 0.1)
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
#' @examples
#' thinned_results <- apply_thinning(mcmc_results, thinning_factor = 2)
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
#' @param data A list with combined results from MCMC.
#' @param prob The probability for the credible interval (between 0 and 1).
#'
#' @examples
#' plot_weibull_distribution(combine_results, prob = 0.95)
#'

plot_penetrance <- function(data, prob = probCI) {
  # Recover the parameters for plotting the Weibull
  params <- calculate_weibull_parameters(
    data$median_results,
    data$first_quartile_results,
    data$shift_results, 
    data$asymptote_results
  )

  alphas <- params$alpha
  betas <- params$beta
  asymptotes <- data$asymptote_results
  shifts <- data$shift_results

  # Define the range for the distribution
  x_values <- seq(0, 100, length.out = 100)

  # Initialize an empty list for distributions
  distributions <- vector("list", length(alphas))

  for (i in seq_along(alphas)) {
    if (validate_weibull_parameters(
      data$first_quartile_results[i],
      data$median_results[i], data$shift_results[i],
      data$asymptote_results[i]
    )) {
      distributions[[i]] <- pweibull(x_values - shifts[i],
        shape = alphas[i], scale = betas[i]
      ) * asymptotes[i]
    } else {
      distributions[[i]] <- rep(NA, length(x_values))
    }
  }

  # Convert list to matrix and check its structure
  distributions_matrix <- do.call(cbind, distributions)

  # Calculate credible intervals with na.rm = TRUE
  lower_bound <- apply(distributions_matrix, 1, quantile, probs = (1 - 0.95) / 2, na.rm = TRUE)
  upper_bound <- apply(distributions_matrix, 1, quantile, probs = 1 - (1 - 0.95) / 2, na.rm = TRUE)
  mean_distribution <- rowMeans(distributions_matrix, na.rm = TRUE)


  # Plot the average distribution
  par(mfrow = c(1, 1))
  plot(x_values, mean_distribution,
    type = "l", col = "blue", lwd = 2,
    xlab = "Age", ylab = "Cumulative Penetrance",
    main = "Penetrance Curve with Credible Intervals"
  )

  # Add credible intervals
  polygon(c(x_values, rev(x_values)), c(lower_bound, rev(upper_bound)), col = rgb(1, 0, 0, 0.2), border = NA)
  lines(x_values, mean_distribution, col = "blue", lwd = 2)

  legend("topleft",
    legend = c("Mean Distribution", "Credible Interval"),
    col = c("blue", "red"), lty = 1, cex = 0.8, fill = c(NA, rgb(1, 0, 0, 0.2))
  )
}