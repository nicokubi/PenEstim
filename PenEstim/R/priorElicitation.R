#' Default Parameter Settings for Prior Distributions
#'
#' Default parameters for the prior distributions used in the \code{makePriors} function.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{asymptote}{A list with components \code{g1} and \code{g2}, default values for the asymptote parameters.}
#'   \item{threshold}{A list with components \code{min} and \code{max}, default values for the threshold parameters.}
#'   \item{median}{A list with components \code{m1} and \code{m2}, default values for the median parameters.}
#'   \item{first_quartile}{A list with components \code{q1} and \code{q2}, default values for the first quartile parameters.}
#' }
#' @export
prior_params_default <- list(
  asymptote = list(g1 = 1, g2 = 1),
  threshold = list(min = 15, max = 25),
  median = list(m1 = 2, m2 = 2),
  first_quartile = list(q1 = 6, q2 = 3)
)

#' Default Risk Proportion Data
#'
#' Default proportions of people at risk used in the \code{makePriors} function.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{median}{Proportion of people at risk at the median age.}
#'   \item{first_quartile}{Proportion of people at risk at the first quartile age.}
#'   \item{max_age}{Proportion of people at risk at the maximum age.}
#' }
#' @export
risk_proportion_default <- data.frame(
  median = 0.5,
  first_quartile = 0.9,
  max_age = 0.1
)

#' Default Distribution Data
#'
#' Default data frame structure with row names for use in the \code{makePriors} function.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{age}{Age values (NA for default).}
#'   \item{at_risk}{Proportion of people at risk (NA for default).}
#' }
#' @export
distribution_data_default <- data.frame(
  row.names = c("min", "first_quartile", "median", "max"),
  age = c(NA, NA, NA, NA),
  at_risk = c(NA, NA, NA, NA)
)
#' Make Priors
#'
#' This function generates prior distributions based on user input or default parameters.
#' It is designed to aid in the statistical analysis of risk proportions in populations, particularly in the context of cancer research.
#' The distributions are calculated for various statistical metrics such as asymptote, threshold, median, and first quartile.
#'
#' @param data A data frame containing age and risk data. If NULL or contains NA values, default parameters are used.
#' @param sample_size Numeric, the total sample size used for risk proportion calculations.
#' @param cancer Character string, specifying the type of cancer, used in OR/RR ratio calculations.
#' @param ratio Numeric, the odds ratio (OR) or relative risk (RR) used in asymptote parameter calculations.
#' @param prior_params List, containing prior parameters for the beta distributions. If NULL, default parameters are used.
#' @param risk_proportion Data frame, with default proportions of people at risk.
#' @param baseline_data Data frame with the baseline risk data. 
#'
#' @details
#' The function includes internal helper functions for normalizing median and first quartile values, and for computing beta distribution parameters.
#' The function handles various settings: using default parameters, applying user inputs, and calculating parameters based on sample size and risk proportions.
#'
#' If the OR/RR ratio is provided, the asymptote parameters are computed based on this ratio, overriding other inputs for the asymptote.
#'
#' The function returns a list of distribution functions for the asymptote, threshold, median, and first quartile, which can be used for further statistical analysis.
#'
#' @return A list of functions representing the prior distributions for asymptote, threshold, median, and first quartile.
#'
#' @seealso \code{\link{qbeta}}, \code{\link{runif}}
#' @export
makePriors <- function(data, sample_size, cancer, ratio, prior_params, risk_proportion, baseline_data) {
  # Helper function definitions
  normalize_median <- function(x) {
    return((x - min_age) / (max_age - min_age))
  }
  
  normalize_first_quartile <- function(x) {
    return((x - min_age) / (median_age - min_age))
  }
  
  compute_parameters_median <- function(stat, at_risk) {
    median_norm <- normalize_median(stat)
    alpha <- median_norm * at_risk
    beta <- at_risk - alpha
    return(list(m1 = alpha, m2 = beta))
  }
  
  compute_parameters_quartile <- function(stat, at_risk) {
    quartile_norm <- normalize_first_quartile(stat)
    alpha <- quartile_norm * at_risk
    beta <- at_risk - alpha
    return(list(q1 = alpha, q2 = beta))
  }
  
  compute_parameters_asymptote <- function(stat, at_risk) {
    max_age_norm <- normalize_median(stat)
    alpha <- max_age_norm * at_risk
    beta <- at_risk - alpha
    return(list(g1 = alpha, g2 = beta))
  }
  
  if (is.null(data) || all(is.na(data))) {
    prior_params <- prior_params_default
  } else {
    if (any(is.na(data$age)) || any(!sapply(data$age, is.numeric))) {
      stop("Missing or non-numeric age entries in the data. Add numeric ages.")
    }
    
    max_age <- data["max", "age"]
    min_age <- data["min", "age"]
    first_quartile_age <- data["first_quartile", "age"]
    median_age <- data["median", "age"]
    
    if (!is.null(data) && all(!is.na(data$age)) && all(is.na(data$at_risk)) && !is.null(sample_size)) {
      risk_median <- risk_proportion$median * sample_size
      risk_first_quartile <- risk_proportion$first_quartile * sample_size
      risk_max_age <- risk_proportion$max_age * sample_size
    } else {
      if (any(is.na(data$at_risk)) || any(!sapply(data$at_risk, is.numeric))) {
        stop("Missing or non-numeric risk entries in the data. Add individuals at risk or total sample size.")
      }
      risk_median <- data$at_risk[data$age == median_age]
      risk_first_quartile <- data$at_risk[data$age == first_quartile_age]
      risk_max_age <- data$at_risk[data$age == max_age]
    }
    
    res_median <- compute_parameters_median(median_age, risk_median)
    res_first_quartile <- compute_parameters_quartile(first_quartile_age, risk_first_quartile)
    res_asymptote <- compute_parameters_asymptote(max_age, risk_max_age)
    
    prior_params <- list(
      asymptote = list(g1 = res_asymptote$g1, g2 = res_asymptote$g2),
      threshold = list(min = 0, max = min_age),
      median = list(m1 = res_median$m1, m2 = res_median$m2),
      first_quartile = list(q1 = res_first_quartile$q1, q2 = res_first_quartile$q2)
    )
  }
  
  if (!is.null(ratio) && !is.null(cancer)) {
  # Calculate the minimal value between Male and Female for each age
  SERR_baseline_min <- pmin(baseline_data_default$Female, baseline_data_default$Male)
    
  # Calculate the cumulative sum of the minimal risk and pick the lifetime risk
  SEER_lifetime <- max(cumsum(SERR_baseline_min))
  prior_params$asymptote <- list(g1 = SEER_lifetime * ratio, g2 = SEER_lifetime * ratio)
  }
  
  asymptote_distribution <- function(n) {
    qbeta(runif(n), prior_params$asymptote$g1, prior_params$asymptote$g2)
  }
  
  threshold_distribution <- function(n) {
    runif(n, prior_params$threshold$min, prior_params$threshold$max)
  }
  
  median_distribution <- function(n) {
    qbeta(runif(n), prior_params$median$m1, prior_params$median$m2)
  }
  
  first_quartile_distribution <- function(n) {
    qbeta(runif(n), prior_params$first_quartile$q1, prior_params$first_quartile$q2)
  }
  
  prior_distributions <- list(
    asymptote_distribution = asymptote_distribution,
    threshold_distribution = threshold_distribution,
    median_distribution = median_distribution,
    first_quartile_distribution = first_quartile_distribution,
    prior_params = prior_params
  )
  
  return(prior_distributions)
}
