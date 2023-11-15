proposal_dist <- list(
  asymptote = list(g1 = 1, g2 = 1),
  shift = list(min = 0, max = 25),
  median = list(m1 = 2, m2 = 2),
  quartile = list(q1 = 6, q2 = 3)
)

create_distributions <- function(
    dataframe = NULL,
    samples_table = NULL,
    sample_size = NULL,
    cancer = NULL,
    ratio = NULL,
    custom_params = proposal_dist,
    median_atrisk = 0.5,
    quartile_atrisk = 0.9,
    max_age_atrisk = 0.1) {
  # Helper function definitions
  # Define the scaling transformation
  normalize_median <- function(x) {
    return((x - min_age) / (max_age - min_age))
  }

  # Define the scaling transformation for first quartile
  normalize_first_quartile <- function(x) {
    return((x - min_age) / (median_age - min_age))
  }

  # Compute beta distribution parameters for median
  compute_parameters_median <- function(stat, at_risk) {
    median_norm <- normalize_median(stat)

    alpha <- median_norm * at_risk
    beta <- at_risk - alpha

    return(list(m1 = alpha, m2 = beta))
  }

  # Compute beta distribution parameters for first quartile
  compute_parameters_quartile <- function(stat, at_risk) {
    quartile_norm <- normalize_first_quartile(stat)

    alpha <- quartile_norm * at_risk
    beta <- at_risk - alpha

    return(list(q1 = alpha, q2 = beta))
  }

  # Compute beta distribution parameters for asymptote parameter
  compute_parameters_asymptote <- function(stat, at_risk) {
    max_age_norm <- normalize_median(stat)
    alpha <- max_age_norm * at_risk
    beta <- at_risk - alpha
    return(list(g1 = alpha, g2 = beta))
  }

  # Main logic for setting the parameters of the proposal distribution
  # Check if dataframe and samples table and ratio are NULL, then use proposal_params as it is
  if (is.null(dataframe) && is.null(samples_table)) {
    # Default parameter settings
    proposal_params <- custom_params
  } else {
    # If dataframe is provided, compute the median and first quartile ages
    if (!is.null(dataframe) && all(c("age", "penetrance_prob") %in% names(dataframe))) {
      max_age <- max(dataframe$age)
      min_age <- min(dataframe$age)
      median_index <- which.min(abs(dataframe$penetrance_prob - 0.5))
      median_age <- dataframe$age[median_index]
      first_quartile <- quantile(dataframe$age, 0.25)
    } else {
      stop("The 'dataframe' argument must contain 'age' and 'penetrance_prob' columns.")
    }
    # If samples table is provided, extract the risk numbers
    if (!is.null(samples_table) && all(c("age", "at_risk") %in% names(samples_table))) {
      risk_median <- ifelse(any(samples_table$age == median_age),
        samples_table$at_risk[samples_table$age == median_age], median_atrisk * sample_size
      )
      risk_quartile <- ifelse(any(samples_table$age == first_quartile),
        samples_table$at_risk[samples_table$age == first_quartile], quartile_atrisk * sample_size
      )
      risk_max_age <- ifelse(any(samples_table$age == max_age),
        samples_table$at_risk[samples_table$age == max_age], max_age_atrisk * sample_size
      )
    } else {
      if (!is.null(sample_size)) {
        risk_median <- median_atrisk * sample_size
        risk_quartile <- quartile_atrisk * sample_size
        risk_max_age <- max_age_atrisk * sample_size
      } else {
        stop("Sample size parameter required.")
      }



      res_median <- compute_parameters_median(median_age, risk_median)
      res_quartile <- compute_parameters_quartile(first_quartile, risk_quartile)
      res_asymptote <- compute_parameters_asymptote(max_age, risk_max_age)

      proposal_params <- list(
        asymptote = list(g1 = res_asymptote$g1, g2 = res_asymptote$g2),
        shift = list(min = 0, max = min_age),
        median = list(m1 = res_median$m1, m2 = res_median$m2),
        quartile = list(q1 = res_quartile$q1, q2 = res_quartile$q2)
      )
    }
  }

  # If Ratio is provided, compute the asymptote parameters based on the ratio
  # This will overwrite any other inputs for the asymptote
  if (!is.null(ratio) && !is.null(cancer)) {
    SERR_baseline <- calculate_lifetime_risk(cancer = cancer, gene = "SEER")
    proposal_params$asymptote <- list(g1 = SERR_baseline * ratio, g2 = SERR_baseline * ratio)
  }


  # Asymptote distribution using either custom or default g1 and g2
  asymptote_distribution <- function(n) {
    qbeta(runif(n), proposal_params$asymptote$g1, proposal_params$asymptote$g2)
  }

  # Shift parameter distribution using either custom or default min and max
  shift_distribution <- function(n) {
    runif(n, proposal_params$shift$min, proposal_params$shift$max)
  }

  # Median distribution using parameters from compute_parameters_median
  median_distribution <- function(n) {
    qbeta(runif(n), proposal_params$median$m1, proposal_params$median$m2)
  }

  # First quartile distribution using parameters from compute_parameters_quartile
  first_quartile_distribution <- function(n) {
    qbeta(runif(n), proposal_params$quartile$q1, proposal_params$quartile$q2)
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
