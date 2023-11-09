create_distributions <- function(dataframe = NULL, samples_table = NULL) {
# Default values
default_max_age <- 94
default_min_age <- 25
default_median_age <- 60
default_first_quartile <- 52
default_risk_median <- 50
default_risk_quartile <- 25

# Default samples_table
if (is.null(samples_table)) {
  samples_table <- data.frame(
    age = seq(default_min_age, default_max_age),
    at_risk = rep(100, default_max_age - default_min_age + 1) # Using 100 as a placeholder risk value
  )
} else if (!all(c("age", "at_risk") %in% names(samples_table))) {
  stop("samples_table must contain 'age' and 'at_risk' columns.")
}

# If dataframe is NULL or doesn't have the necessary columns, use default values
if (is.null(dataframe) || !all(c("age", "penetrance_prob") %in% names(dataframe))) {
  max_age <- default_max_age
  min_age <- default_min_age
  median_age <- default_median_age
  first_quartile <- default_first_quartile
} else {
  # Extract values from dataframe
  max_age <- max(dataframe$age)
  min_age <- min(dataframe$age)

  # Median penetrance time
  median_index <- which.min(dataframe$penetrance_prob <= 0.5)

  if (length(median_index) == 0) {
    median_age <- default_median_age
  } else {
    median_age <- dataframe$age[median_index]
  }

  first_quartile <- quantile(dataframe$age, 0.25)
}

# Extract risk numbers from the samples_table
risk_median <- ifelse(any(samples_table$age == median_age), samples_table$at_risk[samples_table$age == median_age], default_risk_median)
risk_quartile <- ifelse(any(samples_table$age == first_quartile), samples_table$at_risk[samples_table$age == first_quartile], default_risk_quartile)


  # Define the scaling transformation
  normalize_median <- function(x) {
    return((x - min_age) / (max_age - min_age))
  }

  # Define the scaling transformation for first quartile
  normalize_first_quartile <- function(x) {
    return((x - min_age) / (median_age - min_age))
  }


  # Compute beta distribution parameters for median and first quartile$
  # First normalize the estimate of the empirical median, obtained from the data
  # We center the beta distribution around the empirical median
  # Then compute the variance, using binomial proportions

  compute_parameters_median <- function(stat, at_risk) {
    med <- normalize_median(stat)
    var <- med * (1 - med) / at_risk

    alpha <- med * ((med * (1 - med) / var) - 1)
    beta <- (1 - med) * ((med * (1 - med) / var) - 1)

    return(list(alpha = alpha, beta = beta))
  }

  # Compute beta distribution parameters for first quartile
  compute_parameters_quartile <- function(stat, at_risk) {
    med <- normalize_first_quartile(stat)
    var <- med * (1 - med) / at_risk

    alpha <- med * ((med * (1 - med) / var) - 1)
    beta <- (1 - med) * ((med * (1 - med) / var) - 1)

    return(list(alpha = alpha, beta = beta))
  }

  params_median <- compute_parameters_median(median_age, risk_median)
  params_quartile <- compute_parameters_quartile(first_quartile, risk_quartile)


  # Asymptote distribution
  asymptote_distribution <- function(n) {
    qbeta(runif(n), 2, 2, max_age - 5, max_age + 5)
  }

  # Shift parameter distribution
  shift_distribution <- function(n) {
    runif(n, 0, min_age)
  }

  # Median distribution
  median_distribution <- function(n) {
    qbeta(runif(n), params_median$alpha, params_median$beta)
  }

  # First quartile distribution
  first_quartile_distribution <- function(n) {
    qbeta(runif(n), params_quartile$alpha, params_quartile$beta)
  }


  # Create a list with the distributions
  proposal_distributions <- list(
    asymptote_distribution = asymptote_distribution,
    shift_distribution = shift_distribution,
    median_distribution = median_distribution,
    first_quartile_distribution = first_quartile_distribution
  )

  # Return the proposal_distributions list
  return(proposal_distributions)
}
