# Function to perform burn-in
burn_in_chain <- function(samples, burn_in_percentage) {
  n_samples <- length(samples)
  burn_in_count <- round(n_samples * burn_in_percentage)
  samples[(burn_in_count + 1):n_samples]
}

extract_samples <- function(results) {
  list(
    median_samples = do.call(c, lapply(results, function(x) x$median_samples)),
    shift_samples = do.call(c, lapply(results, function(x) x$shift_samples)),
    first_quartile_samples = do.call(c, lapply(results, function(x) x$first_quartile_samples)),
    asymptote_samples = do.call(c, lapply(results, function(x) x$asymptote_samples))
  )
}

# Function to generate summary statistics
generate_summary <- function(data){
  summary_data <- data.frame(
    Median = data$median_samples,
    Shift = data$shift_samples,
    First_Quartile = data$first_quartile_samples,
    Asymptote_Value = data$asymptote_samples
  )
  summary(summary_data)
}

# Function to generate density plots
generate_density_plots <- function(data, title_suffix = "", true_values = c()) {
  par(mfrow = c(2, 2), las = 1, mar = c(5, 4, 4, 2) + 0.1)
  for (name in names(data)) {
    mod_name <- gsub("_", " ", name)
    mod_name <- paste0(toupper(substring(mod_name, 1, 1)), substring(mod_name, 2))
    
    xlim <- if (name %in% c("median_samples", "first_quartile_samples", "shift_samples")) {
      c(0, 100)
    } else if (name == "asymptote_samples") {
      c(0.15, 1)
    } else {
      NULL
    }
    
    hist(data[[name]], 
         main = paste("Density Plot of", mod_name, title_suffix), 
         xlab = mod_name, 
         freq = FALSE,
         xlim = xlim, 
         xaxp = c(min(xlim), max(xlim), 10)  # Specifying 10 intervals
    )
    
    # Adding a vertical line for true values if they exist for the variable
    if (name %in% names(true_values)) {
      abline(v = true_values[name], col = "red", lwd = 2)
    }
  }
}


# Function to generate trace plots for each chain
plot_trace <- function(results, n_chains) {
  par(mfrow = c(n_chains, 4))  # Set up a grid for the plots
  for (chain_id in 1:n_chains) {
    # Extract samples for the current chain
    median_samples <- results[[chain_id]]$median_samples
    shift_samples <- results[[chain_id]]$shift_samples
    first_quartile_samples <- results[[chain_id]]$first_quartile_samples
    asymptote_samples <- results[[chain_id]]$asymptote_samples
    
    # Create trace plots for the current chain
    plot(median_samples, type = "l", main = paste("Chain", chain_id, "- Trace plot of Median"), xlab = "Iteration", ylab = "Median")
    plot(shift_samples, type = "l", main = paste("Chain", chain_id, "- Trace plot of Shift"), xlab = "Iteration", ylab = "Shift")
    plot(first_quartile_samples, type = "l", main = paste("Chain", chain_id, "- Trace plot of First Quartile"), xlab = "Iteration", ylab = "First Quartile")
    plot(asymptote_samples, type = "l", main = paste("Chain", chain_id, "- Trace plot of Asymptote"), xlab = "Iteration", ylab = "Asymptote")
  }
}

# Function to perform burn-in
burn_in_chain <- function(samples, burn_in_percentage) {
  n_samples <- length(samples)
  burn_in_count <- round(n_samples * burn_in_percentage)
  samples[(burn_in_count + 1):n_samples]
}


