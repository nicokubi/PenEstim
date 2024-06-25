#' Calculate Baseline Risk
#'
#' This function extracts the penetrance data for a specified cancer type, gene, race, and penetrance type
#' from the provided database.
#'
#' @param cancer_type The type of cancer for which the risk is being calculated.
#' @param gene The gene of interest for which the risk is being calculated.
#' @param race The race of the individual.
#' @param type The type of penetrance calculation.
#' @param db The dataset used for the calculation, containing penetrance data.
#'
#' @return A matrix of penetrance data for the specified parameters.
#' @export
calculateBaseline <- function(cancer_type, gene, race, type, db) {
    # Check if dimnames are available and correct
    if (is.null(db$Penetrance) || is.null(attr(db$Penetrance, "dimnames"))) {
        stop("Penetrance data or its dimension names are not properly defined.")
    }

    dim_names <- attr(db$Penetrance, "dimnames")
    required_dims <- c("Cancer", "Gene", "Race", "Age", "PenetType")
    if (!all(required_dims %in% names(dim_names))) {
        stop("One or more required dimensions are missing in Penetrance data.")
    }

    # Function to safely extract index
    get_index <- function(dim_name, value) {
        idx <- which(dim_names[[dim_name]] == value)
        if (length(idx) == 0) {
            stop(paste("Value", value, "not found in dimension", dim_name))
        }
        idx
    }

    # Extracting indices for each dimension except Age
    cancer_index <- get_index("Cancer", cancer_type)
    gene_index <- get_index("Gene", gene)
    race_index <- get_index("Race", race)
    type_index <- get_index("PenetType", type)

    # Subsetting Penetrance data for all ages using indices
    lifetime_risk <- db$Penetrance[cancer_index, gene_index, race_index, , , type_index]
    return(lifetime_risk)
}

#' Calculate Age-Specific Non-Carrier Penetrance
#'
#' This function calculates the age-specific non-carrier penetrance based on SEER baseline
#' data, penetrances for carriers, and allele frequencies. It adjusts penetrance estimates
#' for genetic testing by incorporating the genetic risk attributable to specified alleles.
#'
#' @param SEER_baseline Numeric, the baseline penetrance derived from SEER data for the general population without considering genetic risk factors.
#' @param alpha Numeric, shape parameter for the Weibull distribution used to model carrier risk.
#' @param beta Numeric, scale parameter for the Weibull distribution used to model carrier risk.
#' @param delta Numeric, location parameter for the Weibull distribution used to model carrier risk.
#' @param gamma Numeric, scaling factor applied to the Weibull distribution to adjust carrier risk.
#' @param af Numeric, the allele frequency of the risk allele in the population.
#' @param max_age Integer, the maximum age up to which the calculations are performed.
#'
#' @return A list containing:
#' \item{weightedCarrierRisk}{Numeric vector, the weighted risk for carriers at each age based on allele frequency.}
#' \item{yearlyProb}{Numeric vector, the yearly probability of not getting the disease at each age.}
#' \item{cumulativeProb}{Numeric vector, the cumulative probability of not getting the disease up to each age.}
#'
calculateNCPen <- function(SEER_baseline, alpha, beta, delta, gamma, af, max_age) {
  # Calculate probability weights for carriers based on allele frequencies
  weights <- 2 * af * (1 - af) # Heterozygous carriers only
  
  # Initialize vectors to store the yearly and cumulative probability of not getting the disease
  weightedCarrierRisk <- numeric(max_age)
  yearlyProb <- numeric(max_age) # For single-year probability
  cumulativeProb <- numeric(max_age) # For cumulative probability
  
  # Start with 100% probability of not having the disease
  cumulativeProbability <- 1
  
  for (age in 1:max_age) {
    # Calculate the risk for carriers at this age
    carrierRisk <- dweibull(age - delta, shape = alpha, scale = beta) * gamma
    # Calculate the weighted risk for carriers based on allele frequency
    weightedCarrierRisk[age] <- carrierRisk * weights
    
    # Calculate the single-year probability of not getting the disease
    yearlyProb[age] <- 1 - weightedCarrierRisk[age]
    
    # Update cumulative probability of not getting the disease
    cumulativeProbability <- cumulativeProbability * yearlyProb[age]
    cumulativeProb[age] <- cumulativeProbability
  }
  
  # Return both yearly and cumulative probabilities
  return(list(
    weightedCarrierRisk = weightedCarrierRisk,
    yearlyProb = yearlyProb, cumulativeProb = cumulativeProb
  ))
}

#' Penetrance Function
#'
#' Calculates the penetrance for an individual based on Weibull distribution parameters.
#' This function estimates the probability of developing cancer given the individual's genetic and demographic information.
#'
#' @param i Integer, index of the individual in the data set.
#' @param data Data frame, containing individual demographic and genetic information. Must include columns for 'sex', 'age', 'aff' (affection status), and 'geno' (genotype).
#' @param alpha_male Numeric, Weibull distribution shape parameter for males.
#' @param alpha_female Numeric, Weibull distribution shape parameter for females.
#' @param beta_male Numeric, Weibull distribution scale parameter for males.
#' @param beta_female Numeric, Weibull distribution scale parameter for females.
#' @param delta_male Numeric, shift parameter for the Weibull function for males.
#' @param delta_female Numeric, shift parameter for the Weibull function for females.
#' @param gamma_male Numeric, asymptote parameter for males (only scales the entire distribution).
#' @param gamma_female Numeric, asymptote parameter for females (only scales the entire distribution).
#' @param max_age Integer, maximum age considered in the analysis.
#' @param baselineRisk Numeric matrix, baseline risk for each age by sex. Rows correspond to sex (1 for male, 2 for female) and columns to age.
#' @param BaselineNC Logical, indicates if non-carrier penetrance should be based on SEER data.
#'
#' @return Numeric vector, containing penetrance values for unaffected and affected individuals.
#'
lik.fn <- function(i, data, alpha_male, alpha_female, beta_male, beta_female,
                   delta_male, delta_female, gamma_male, gamma_female, max_age,
                   baselineRisk, BaselineNC) {
  # Map sex to row index: "Female" is 1st row and "Male" is 2nd row
  sex_index <- ifelse(data$sex[i] == 2, "Female", "Male")
  
  # Select parameters based on individual's sex
  alpha <- ifelse(data$sex[i] == 1, alpha_male, alpha_female)
  beta <- ifelse(data$sex[i] == 1, beta_male, beta_female)
  gamma <- ifelse(data$sex[i] == 1, gamma_male, gamma_female)
  delta <- ifelse(data$sex[i] == 1, delta_male, delta_female)
  
  if (data$age[i] == 0 || data$age[i] == 1) {
    lik.i <- c(1, 1) # Assuming people aged 0 or 1 years are all unaffected
  } else {
    # Ensure age is within the valid range
    age_index <- min(max_age, data$age[i])
    
    # Weibull parameters for penetrance, using sex-specific gamma
    survival_prob <- 1 - pweibull(max(age_index - delta, 1), shape = alpha, scale = beta) * gamma
    c.pen <- (pweibull(max(age_index - delta, 1), shape = alpha, scale = beta)
              - pweibull(max(age_index - 1 - delta, 1), shape = alpha, scale = beta)) * gamma
    
    # Extract the corresponding baseline risk for sex and age
    SEER_baseline_max <- baselineRisk[1:age_index, sex_index]
    SEER_baseline_cum <- cumsum(baselineRisk[, sex_index])[age_index]
    SEER_baseline_i <- baselineRisk[age_index, sex_index]
    
    # Calculate cumulative risk for non-carriers based on SEER data or other model
    if (BaselineNC == TRUE) {
      nc.pen <- SEER_baseline_i
      nc.pen.c <- prod(1 - SEER_baseline_i)
    } else {
      nc.pen <- calculateNCPen(
        SEER_baseline = SEER_baseline_max, alpha = alpha,
        beta = beta, delta = delta, gamma = gamma, max_age = max_age
      )$weightedCarrierRisk[age_index]
      nc.pen.c <- calculateNCPen(
        SEER_baseline = SEER_baseline_max, alpha = alpha,
        beta = beta, delta = delta, gamma = gamma, max_age = max_age
      )$cumulativeProb[age_index]
    }
    
    # Penetrance calculations based on genotype and affection status
    lik.i <- c(nc.pen.c, survival_prob) # for censored observations
    if (data$aff[i] == 1) lik.i <- c(nc.pen * nc.pen.c, c.pen) # for affected observations
  }
  
  # Adjustment for observed genotypes
  if (data$geno[i] == "1/1") lik.i[-1] <- 0
  if (data$geno[i] == "1/2") lik.i[-2] <- 0
  
  return(lik.i)
}

#' Calculate Log Likelihood using clipp Package
#'
#' This function calculates the log likelihood using the clipp package for a set of parameters and data.
#'
#' @param paras Numeric vector, the parameters for the Weibull distribution and scaling factors. 
#'        Should contain in order: gamma_male, gamma_female, delta_male, delta_female, 
#'        given_median_male, given_median_female, given_first_quartile_male, given_first_quartile_female.
#' @param families Data frame, containing pedigree information with columns for 'sex', 'age', 'aff' (affection status), and 'geno' (genotype).
#' @param max_age Integer, maximum age considered in the analysis.
#' @param seer_data Numeric matrix, baseline risk data for each age by sex. Rows correspond to sex (1 for male, 2 for female) and columns to age.
#' @param af Numeric, allele frequency of the risk allele in the population.
#' @param BaselineNC Logical, indicates if non-carrier penetrance should be based on the baseline data or the calculated non-carrier penetrance.
#' @param ncores Integer, number of cores to use for parallel computation.
#'
#' @return Numeric, the calculated log likelihood.
#'
#' @references
#' Details about the clipp package and methods can be found in the package documentation.
#'
#' @import clipp
mhLogLikelihood_clipp <- function(paras, families, max_age, baseline_data, af, BaselineNC, ncores) {
  paras <- unlist(paras)
    # Extract parameters
    gamma_male <- paras[1]
    gamma_female <- paras[2]
    delta_male <- paras[3]
    delta_female <- paras[4]
    given_median_male <- paras[5]
    given_median_female <- paras[6]
    given_first_quartile_male <- paras[7]
    given_first_quartile_female <- paras[8]

    # Calculate Weibull parameters
    params_male <- calculate_weibull_parameters(given_median_male, given_first_quartile_male, delta_male)
    alpha_male <- params_male$alpha
    beta_male <- params_male$beta

    params_female <- calculate_weibull_parameters(given_median_female, given_first_quartile_female, delta_female)
    alpha_female <- params_female$alpha
    beta_female <- params_female$beta

    # Initialize the model
    # Here we assume a simple biallelic model where af is the allele frequency of the
    #  pathogenic variant
    geno_freq <- c(1 - af, af)

    # The transition matrix represents the probabilities of transitioning between
    # different genotype states across generations.
    trans <- matrix(
        c(
            1, 0,
            0.5, 0.5,
            0.5, 0.5,
            1 / 3, 2 / 3
        ),
        nrow = 4, ncol = 2, byrow = TRUE
    )

    baselineRisk <- baseline_data

    # Calculate penetrance
    lik <- t(sapply(1:nrow(families), function(i) {
        lik.fn(i, families, alpha_male, alpha_female, beta_male, beta_female, delta_male, delta_female, gamma_male, gamma_female,
            max_age, baselineRisk, BaselineNC
        )
    }))

    # Compute log-likelihood
    loglik <- pedigree_loglikelihood(families, geno_freq, trans, lik, ncores = ncores)

    # Handle -Inf values
    if (is.infinite(loglik) && loglik == -Inf) {
        loglik <- -50000
    }
    return(loglik)
}
