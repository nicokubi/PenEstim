#' Make Priors 
#'
#' This function generates prior distributions based on user input or default parameters. It is designed to aid in the statistical analysis of risk proportions in populations, particularly in the context of cancer research. The distributions are calculated for various statistical metrics such as asymptote, threshold, median, and first quartile.
#'
#' @param data A data frame containing age and risk data. If NULL or contains NA values, default parameters are used.
#' @param sample_size The total sample size used for risk proportion calculations.
#' @param cancer A character string specifying the type of cancer, used in OR/RR ratio calculations.
#' @param ratio The odds ratio (OR) or relative risk (RR) used in asymptote parameter calculations.
#' @param prior_params A list of prior parameters for the beta distributions. If NULL, default parameters are used.
#' @param risk_proportion A data frame with default proportions of people at risk.
#'
#' @details
#' The function includes internal helper functions for normalizing median and first quartile values, and for computing beta distribution parameters. The function handles various settings: using default parameters, applying user inputs, and calculating parameters based on sample size and risk proportions.
#'
#' If the OR/RR ratio is provided, the asymptote parameters are computed based on this ratio, overriding other inputs for the asymptote.
#'
#' The function returns a list of distribution functions for the asymptote, threshold, median, and first quartile, which can be used for further statistical analysis.
#'
#' @return A list of functions representing the prior distributions for asymptote, threshold, median, and first quartile.
#'
#' @examples
#' # Example usage with default parameters
#' distributions <- create_distributions(NULL, 100, "breast_cancer", 1.2, NULL, risk_proportion_default)
#'
#' @seealso \code{\link{qbeta}}, \code{\link{runif}}
#' @export

# Default parameter settings
prior_params_default <- list(
  asymptote = list(g1 = 1, g2 = 1),
  threshold = list(min = 15, max = 25),
  median = list(m1 = 2, m2 = 2),
  first_quartile = list(q1 = 6, q2 = 3)
)

# Dataframe with the default proportions of poeple at risk
risk_proportion_default <- data.frame(
  median = 0.5,
  first_quartile = 0.9,
  max_age = 0.1
)

# Create a data frame with row names
distribution_data_default <- data.frame(
  row.names = c("min", "first_quartile", "median", "max"),
  age = c(NA, NA, NA, NA),
  at_risk = c(NA, NA, NA, NA)
)

#  Function to create the prior distributions
makePriors <- function(
    data,
    sample_size,
    cancer,
    ratio,
    prior_params,
    risk_proportion) {
  # Helper function definitions
  # Define the scaling transformation
  normalize_median <- function(x) {
    return((x - min_age) / (max_age - min_age))
  }

  # Define the scaling transformation for first quartile
  normalize_first_quartile <- function(x) {
    return((x - min_age) / (median_age - min_age))
  }

  # Function: Compute beta distribution parameters for median
  compute_parameters_median <- function(stat, at_risk) {
    median_norm <- normalize_median(stat)
    alpha <- median_norm * at_risk
    beta <- at_risk - alpha
    return(list(m1 = alpha, m2 = beta))
  }

  # Function: Compute beta distribution parameters for first quartile
  compute_parameters_quartile <- function(stat, at_risk) {
    quartile_norm <- normalize_first_quartile(stat)
    alpha <- quartile_norm * at_risk
    beta <- at_risk - alpha
    return(list(q1 = alpha, q2 = beta))
  }

  # Function: Compute beta distribution parameters for asymptote parameter
  compute_parameters_asymptote <- function(stat, at_risk) {
    max_age_norm <- normalize_median(stat)
    alpha <- max_age_norm * at_risk
    beta <- at_risk - alpha
    return(list(g1 = alpha, g2 = beta))
  }
  
  compute_parameters_asymptote <- function(stat, at_risk) {
  max_age_norm <- normalize_median(stat)
  alpha <- max_age_norm * at_risk
  beta <- alpha
  return(list(g1 = alpha, g2 = beta))
  }

  # Main logic for setting the parameters of the prior distribution
  # Setting 1: When there is no user input in the data_distribution, then the default parameter settings are applied.
  # Setting 2: If the user has modified the prior_params_default object, the customized parameter settings will be applied.
  if (is.null(data) || all(is.na(data))) {
    prior_params <- prior_params
  } else {
    # Setting 3: Extracting the user inputs from the data_distribution_default dataframe for the prior elicitation.
    # Check if all age entries are present
    if (any(is.na(data$age)) || any(!sapply(data$age, is.numeric))) {
      stop("Missing or non-numeric age entries in the data. Add numeric ages.")
    }
    # Extract the ages from the user input.
    max_age <- data["max", "age"]
    min_age <- data["min", "age"]
    first_quartile_age <- data["first_quartile", "age"]
    median_age <- data["median", "age"]

    # Setting 3b: If the user did not provide the number of people at risk, the risk proportions are applied to calculate the number of people at risk
    if (!is.null(data) && all(!is.na(data$age)) && all(is.na(data$at_risk)) && !is.null(sample_size)) {
      risk_median <- risk_proportion$median * sample_size
      risk_first_quartile <- risk_proportion$first_quartile * sample_size
      risk_max_age <- risk_proportion$max_age * sample_size
    } else {
      if (any(is.na(data$at_risk)) || any(!sapply(data$at_risk, is.numeric))) {
        stop("Missing or non-numeric risk entries in the data. Add individuals at risk or total sample size.")
      }
      # Setting 3a: If all the data is provided, extract the number of people at risk from the dataframe
      risk_median <- data$at_risk[data$age == median_age]
      risk_first_quartile <- data$at_risk[data$age == first_quartile_age]
      risk_max_age <- data$at_risk[data$age == max_age]
    }

    # Calculate the parameters based on the extracted statistics above
    res_median <- compute_parameters_median(median_age, risk_median)
    res_first_quartile <- compute_parameters_quartile(first_quartile_age, risk_first_quartile)
    res_asymptote <- compute_parameters_asymptote(max_age, risk_max_age)

    # Update the prior_params object with the newly calculated parameters
    prior_params <- list(
      asymptote = list(g1 = res_asymptote$g1, g2 = res_asymptote$g2),
      threshold = list(min = 0, max = min_age),
      median = list(m1 = res_median$m1, m2 = res_median$m2),
      first_quartile = list(q1 = res_first_quartile$q1, q2 = res_first_quartile$q2)
    )
  }
  # Setting 3c: If OR/RR ratio is provided, compute the asymptote parameters based on the ratio
  # This will overwrite any other inputs for the asymptote
  if (!is.null(ratio) && !is.null(cancer)) {
    SERR_baseline <- calculate_lifetime_risk(cancer = cancer, gene = "SEER")
    prior_params$asymptote <- list(g1 = SERR_baseline * ratio, g2 = SERR_baseline * ratio)
  }

  # Asymptote distribution using either custom or default g1 and g2
  asymptote_distribution <- function(n) {
    qbeta(runif(n), prior_params$asymptote$g1, prior_params$asymptote$g2)
  }

  # threshold parameter distribution using either custom or default min and max
  threshold_distribution <- function(n) {
    runif(n, prior_params$threshold$min, prior_params$threshold$max)
  }

  # Median distribution using parameters from compute_parameters_median
  median_distribution <- function(n) {
    qbeta(runif(n), prior_params$median$m1, prior_params$median$m2)
  }

  # First quartile distribution using parameters from compute_parameters_quartile
  first_quartile_distribution <- function(n) {
    qbeta(runif(n), prior_params$first_quartile$q1, prior_params$first_quartile$q2)
  }

  # Create a list with the distributions
  prior_distributions <- list(
    asymptote_distribution = asymptote_distribution,
    threshold_distribution = threshold_distribution,
    median_distribution = median_distribution,
    first_quartile_distribution = first_quartile_distribution
  )

  # Return the prior_distributions list
  return(prior_distributions)
}