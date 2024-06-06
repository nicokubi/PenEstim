#' Calculate Age-Specific and Cumulative Non-Carrier Penetrance
#'
#' This function calculates the age-specific and cumulative non-carrier penetrance based on SEER baseline
#' data, genetic risk parameters, and allele frequencies. It accommodates the inclusion of homozygous
#' carriers in the calculation.
#'
#' @param SEER_baseline Numeric vector, the baseline risk for the general population by age.
#' @param alpha Numeric, the shape parameter of the Weibull distribution for the genetic risk.
#' @param beta Numeric, the scale parameter of the Weibull distribution for the genetic risk.
#' @param delta Numeric, the shift parameter for the age in the Weibull risk calculation.
#' @param gamma Numeric, the adjustment factor for the genetic risk.
#' @param af Numeric, the allele frequency of the risk allele in the population.
#' @param max_age Integer, the maximum age to calculate the penetrance for.
#' @param homozygous Logical, indicates if homozygous carriers should be included. Defaults to FALSE.
#'
#' @return A list containing two elements: 'yearlyProb' with the yearly probability of not developing the disease,
#'         and 'cumulativeProb' with the cumulative probability up to each age.
#'
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
#' data, penetrances for carriers, allele frequencies, and an option to include homozygous
#' carriers. It is designed to adjust penetrance estimates for genetic testing by incorporating
#' the genetic risk attributable to specified alleles.
#'
#' @param SEER_baseline Numeric, the baseline penetrance derived from SEER data for the general population without considering genetic risk factors.
#' @param carrierPenetrances Numeric vector, the penetrances for carriers at different ages and possibly by sex, assuming that the structure aligns with the age and possibly sex.
#' @param af Numeric, the allele frequency of the risk allele in the population.
#' @param homozygous Logical, whether to include the weights for homozygous carriers in the calculation. Defaults to FALSE.
#'
#' @return Numeric, the calculated age-specific non-carrier penetrance, adjusted for genetic risk.
#'
#' @examples
#' SEER_baseline <- 0.05 # Example baseline penetrance
#' carrierPenetrances <- c(0.1, 0.2, 0.3) # Example carrier penetrances
#' af <- 0.01 # Example allele frequency
#' calculateNCPen(SEER_baseline, carrierPenetrances, af)
#' calculateNCPen(SEER_baseline, carrierPenetrances, af, homozygous = TRUE)
#'
#' @export
#'
calculateNCPen <- function(SEER_baseline, alpha, beta, delta, gamma, af, max_age, homozygote) {
    # Calculate probability weights for carriers based on allele frequencies
    weights1 <- 2 * af * (1 - af) # Heterozygous carriers
    weights2 <- af^2 # Homozygous carriers, if considered
    weights <- if (homozygote) weights1 + weights2 else weights1

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
#' @param data Data frame, containing individual demographic and genetic information.
#' @param alpha Numeric, Weibull distribution shape parameter.
#' @param beta Numeric, Weibull distribution scale parameter.
#' @param delta Numeric, shift parameter for the Weibull function.
#' @param gamma_male Numeric, asymptote parameter, for males (only scales the entire distribution)
#' @param gamma_female Numeric, asymptote parameter, for females (only scales the entire distribution)
#' @param max_age Integer, maximum age considered in the analysis.
#' @param baselineRisk Numeric vector, baseline risk for each age.
#' @param homozygote Logical, indicates if homozygous carriers should be considered.
#' @param SeerNC Logical, indicates if non-carrier penetrance should be based on SEER data.
#' @param sex Character, the sex of the individuals ('Male', 'Female', or 'NA').
#'
#' @return Numeric vector, containing penetrance values for unaffected and affected individuals.
#' @export
#'
lik.fn <- function(i, data, alpha_male, alpha_female, beta_male, beta_female, 
                   delta_male, delta_female, gamma_male, gamma_female, max_age,
                   baselineRisk, homozygote, SeerNC) {
   
  # Map sex to row index: "Female" is 1st row and "Male" is 2nd row
    sex_index <- ifelse(data$sex[i] == 2, 1, 2)

    # Select gamma based on individual's sex
    # in data the 1 indicates male, 2 indicates female
    alpha <- ifelse(data$sex[i] == 1, alpha_male, alpha_female)
    beta <- ifelse(data$sex[i] == 1, beta_male, beta_female)
    gamma <- ifelse(data$sex[i] == 1, gamma_male, gamma_female)
    delta <- ifelse(data$sex[i] == 1, delta_male, delta_female)

    if (data$age[i] == 0 || data$age[i] == 1)  {
        lik.i <- c(1, 1) # Assuming people aged 0 or 1 years are all unaffected
    } else {
        # Ensure age is within the valid range
        age_index <- min(max_age, data$age[i])

        # Weibull parameters for penetrance, using sex-specific gamma
        survival_prob <- 1 - pweibull(max(age_index - delta, 1), shape = alpha, scale = beta) * gamma
        c.pen <- (pweibull(max(age_index - delta, 1), shape = alpha, scale = beta)
        - pweibull(max(age_index - 1 - delta, 1), shape = alpha, scale = beta)) * gamma

        # Extract the corresponding baseline risk for sex and age
        SEER_baseline_max <- baselineRisk[sex_index, 1:age_index]
        SEER_baseline_cum <- cumsum(baselineRisk[sex_index,])[age_index]
        SEER_baseline_i <- baselineRisk[sex_index, age_index]
        
        SEER_baseline_max <- baselineRisk[1, 1:80]
        SEER_baseline_cum <- cumsum(baselineRisk[1,])[80]
        SEER_baseline_i <- baselineRisk[1, 80]

        # Calculate cumuative risk for non-carriers based on SEER data or other model
        if (SeerNC == TRUE) {
            nc.pen <- SEER_baseline_i
            nc.pen.c <- prod(1- SEER_baseline_i)
        } else {
            nc.pen <- calculateNCPen(
                SEER_baseline = SEER_baseline_max, alpha = alpha,
                beta = beta, delta = delta, gamma = gamma, max_age = max_age,
                homozygote = homozygote
            )$weightedCarrierRisk[age_index]
            nc.pen.c <- calculateNCPen(
                SEER_baseline = SEER_baseline_max, alpha = alpha,
                beta = beta, delta = delta, gamma = gamma, max_age = max_age,
                homozygote = homozygote
            )$cumulativeProb[age_index]
        }

        # Penetrance calculations based on genotype and affection status
        lik.i <- c(nc.pen.c, survival_prob) # for censored observations
        if (data$aff[i] == 1) lik.i <- c(nc.pen * nc.pen.c, c.pen) # for affected observations
    }

    # Adjustment for observed genotypes
    if (data$geno[i] == "1/1") lik.i[-1] <- 0
    if (data$geno[i] == "1/2") lik.i[-2] <- 0
    #if (data$geno[i] == "2/2") lik.i[-3] <- 0

    # Setting the third element to 0 if homozygote is FALSE
    #if (!homozygote) {
       # lik.i[3] <- 0
   # }

    return(lik.i)
}

#' Calculate Log Likelihood using clipp Package
#'
#' This function calculates the log likelihood using the clipp package for a set of parameters and data.
#'
#' @param paras A vector of parameters (proposals)
#' @param families The pedigree data for families.
#' @param max_age The maximum age considered.
#' @param cancer_type The type of cancer.
#' @param db The data object.
#' @param af Allele frequency.
#'
#' @return The log likelihood value.
#'
#' @references
#' Details about the clipp package and methods can be found in the package documentation.
#'
#' @examples
#' # Calculate log likelihood
#' log_likelihood <- mhLogLikelihood_clipp(parameters, pedigree_data, 80, "Lung Cancer", data_object, 0.2)
#' 
#' ¨

mhLogLikelihood_clipp <- function(paras, families, max_age, cancer_type, db, af, homozygote, SeerNC) {

  paras<- unlist(paras)
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
    alpha_male <-  params_male$alpha
    beta_male <- params_male$beta

    params_female <- calculate_weibull_parameters(given_median_female, given_first_quartile_female, delta_female)
    alpha_female <- params_female$alpha
    beta_female <- params_female$beta

    # Initialize the model
    
    if (homozygote) {
    geno_freq <- c(1 - af, af)
    trans <- matrix(c(1, 0,
                         0.5, 0.5,
                         0.5, 0.5,
                         1/3, 2/3), 
                       nrow = 4, ncol = 2, byrow = TRUE)
    } else {
      geno_freq <- geno_freq_monogenic(p_alleles = c(1 - af, af))
      trans <- trans_monogenic(n_alleles = 2)
      
    }
    
    baselineRisk <- calculateBaseline(
        cancer_type = cancer_type,
        gene = "SEER", race = "All_Races", type = "Net", db
    )

    # Calculate penetrance
    lik <- t(sapply(1:nrow(families), function(i) {
        lik.fn(i, families, alpha_male, alpha_female, beta_male, beta_female, delta_male, delta_female, gamma_male, gamma_female,
            max_age, baselineRisk,
            homozygote = homozygote, SeerNC = SeerNC
        )
    }))

    # Compute log-likelihood
    loglik <- pedigree_loglikelihood(families, geno_freq, trans, lik, ncores = 6) 

    # Handle -Inf values
   if (is.infinite(loglik) && loglik == -Inf) {
            loglik <- -50000
    } 
    return(loglik)
}
