create_distributions <- function(dataframe, samples_table) {
  # Default values
  default_max_age <- 94
  default_min_age <- 25
  default_median_age <- 60
  default_first_quartile <- 52
  default_risk_median <- 50
  default_risk_quartile <- 25

  # Check if necessary columns exist in dataframe
  if (!all(c("age", "penetrance_prob") %in% names(dataframe))) {
    stop("Dataframe must contain 'age' and 'penetrance_prob' columns.")
  }

  if (!all(c("age", "at_risk") %in% names(samples_table))) {
    stop("samples_table must contain 'age' and 'at_risk' columns.")
  }

  # Extract values from dataframe
  max_age <- if (nrow(dataframe) > 0) max(dataframe$age) else default_max_age
  min_age <- if (nrow(dataframe) > 0) min(dataframe$age) else default_min_age

  # Median penetrance time
  median_index <- which.min(dataframe$penetrance_prob <= 0.5)

  if (length(median_index) == 0) {
    median_age <- default_median_age
  } else {
    median_age <- dataframe$age[median_index]
  }

  first_quartile <- if (nrow(dataframe) > 0) quantile(dataframe$age, 0.25) else default_first_quartile

  # Extract risk numbers from the samples_table
  risk_median <- if (nrow(samples_table) > 0 && any(samples_table$age == median_age)) {
    samples_table$at_risk[samples_table$age == median_age]
  } else {
    default_risk_median
  }

  risk_quartile <- if (nrow(samples_table) > 0 && any(samples_table$age == first_quartile)) {
    samples_table$at_risk[samples_table$age == first_quartile]
  } else {
    default_risk_quartile
  }

  # Define the scaling transformation
  normalize <- function(x) {
    return((x - min_age) / (max_age - min_age))
  }

  denormalize <- function(x) {
    return(min_age + (max_age - min_age) * x)
  }

  # Define the scaling transformation for first quartile
  normalize_quartile <- function(x) {
    return((x - min_age) / (median_age - min_age))
  }

  denormalize_quartile <- function(x) {
    return(min_age + (median_age - min_age) * x)
  }

  # Compute beta distribution parameters for median and first quartile$
  #  First normalize the estimate of the empirical median, obtained from the data
  # We center the beta distribution around the empirical median
  #  Then compute the variance, using binomial proportions
  compute_parameters <- function(stat, at_risk) {
    med <- normalize(stat)
    var <- med * (1 - med) / at_risk

    alpha <- med * ((med * (1 - med) / var) - 1)
    beta <- (1 - med) * ((med * (1 - med) / var) - 1)

    return(list(alpha = alpha, beta = beta))
  }

  params_median <- compute_parameters(median_age, risk_median)
  params_quartile <- compute_parameters(first_quartile, risk_quartile)

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
    denormalize(qbeta(runif(n), params_median$alpha, params_median$beta))
  }

  # First quartile distribution
  first_quartile_distribution <- function(n) {
    denormalize_quartile(qbeta(runif(n), params_quartile$alpha, params_quartile$beta))
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
