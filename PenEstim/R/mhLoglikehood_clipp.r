#' Calculate Baseline Penetrance
#'
#' This function calculates the baseline penetrance for a specific cancer type, gene, race, sex, and type.
#' It extracts the necessary information from the provided data object.
#'
#' @param cancer_type The type of cancer.
#' @param gene The gene (default is "SEER").
#' @param race The race (default is "All_Races").
#' @param sex The sex (default is "Male").
#' @param type The type of penetrance (default is "Crude").
#' @param data The data object containing penetrance information.
#'
#' @return A matrix containing the baseline penetrance values.
#'
#' @seealso \code{\link{penet.fn}}
#'
#' @examples
#' # Calculate baseline penetrance
#' baseline <- calculateBaseline("Lung Cancer", data = my_data)
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
#' This function calculates the penetrance for an individual based on a Weibull distribution.
#' Given this simple penetrance function we assume a single age-specific penetrance function.
#'
#' @param i The index of the individual.
#' @param data The data object containing individual information.
#' @param alpha Weibull distribution parameter alpha.
#' @param beta Weibull distribution parameter beta.
#' @param delta The delta parameter.
#' @param gamma The gamma parameter.
#' @param max_age The maximum age considered.
#' @param baselineRisk The baseline penetrance values.
#'
#' @return A vector containing penetrance values for unaffected and affected individuals.
#'
#' @examples
#' # Calculate penetrance for an individual
#' penetrance <- penet.fn(1, individual_data, 1.0, 2.0, 0.5, 0.8, 80, baseline_risk)
#'
# Define the penetrance function
penet.fn <- function(i, data, alpha, beta, delta, gamma, max_age, baselineRisk, homozygote, SeerNC, sex) {
    # Map sex to row index: "Female" is 1st row and "Male" is 2nd row
    sex_index <- ifelse(data$sex[i] == 2, 1, 2)

    if (data$age[i] == 0) {
        penet.i <- c(1, 1, 1) # Assuming people aged 0 years are all unaffected
    } else {
        # Ensure age is within the valid range
        age_index <- min(max_age, data$age[i])

        # Weibull with four parameters for the penetrance
        # For now we just have one penetrance, irrespective of sex
        c.pen <- dweibull(data$age[i] - delta, alpha, beta) * gamma
        # Cumulative penetrance calculation using product of 1-year survival probabilities
        c.pen.c <- exp(-gamma * pweibull(data$age[i] - delta, alpha, beta, lower.tail = FALSE))

        # Extract the corresponding baseline risk for sex and age
        SEER_baseline_max <- baselineRisk[sex_index, 1:max_age]
        SEER_baseline_cum <- baselineRisk[sex_index, 1:age_index]
        SEER_baseline_i <- baselineRisk[sex_index, age_index]

        # Calculate cumulative risk for non-carriers based on SEER data or other model
        if (SeerNC == TRUE) {
            nc.pen <- SEER_baseline_i
            # calculative cumulative product for being cancer-free for non-carriers
            nc.pen.c <- prod(1 - SEER_baseline_cum)
        } else {
            nc.pen <- calculateNCPen(
                SEER_baseline = SEER_baseline_max, alpha = alpha,
                beta = beta, delta = delta, gamma = gamma, , af = af, max_age = max_age,
                homozygote = homozygote
            )$weightedCarrierRisk[age_index]
            # calculative cumulative product for being cancer-free for non-carriers
            nc.pen.c <- calculateNCPen(
                SEER_baseline = SEER_baseline_max, alpha = alpha,
                beta = beta, delta = delta, gamma = gamma, , af = af, max_age = max_age,
                homozygote = homozygote
            )$cumulativeProb[age_index]
        }

        # Penetrance calculations based on genotype
        # Assuming same penetrance for homozygous and heterozygous carriers
        penet.i <- c(1 - nc.pen.c, 1 - c.pen.c, 1 - c.pen.c) # if person is not affected
        if (data$aff[i] == 1) penet.i <- c(nc.pen, c.pen, c.pen)
    }

    if (data$geno[i] == "1/1") penet.i[-1] <- 0
    if (data$geno[i] == "1/2") penet.i[-2] <- 0
    if (data$geno[i] == "2/2") penet.i[-3] <- 0

    # Setting the third element to 0 if homozygote is FALSE
    if (!homozygote) {
        penet.i[3] <- 0
    }

    # Adjusting for gender
    pen <- 1e-28
    if (sex == "Male" && data$sex[i] != 2) penet.i <- c(pen, pen, pen)
    if (sex == "Female" && data$sex[i] != 1) penet.i <- c(pen, pen, pen)
    if (sex == "NA") penet.i <- penet.i

    return(penet.i)
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
mhLogLikelihood_clipp <- function(paras, families, max_age, cancer_type, db, af, homozygote, SeerNC, sex) {
    # Extract parameters
    given_median <- paras[1]
    given_first_quartile <- paras[2]
    delta <- paras[3]
    gamma <- paras[4]

    # Calculate Weibull parameters
    params <- calculate_weibull_parameters(given_median, given_first_quartile, delta, gamma)
    alpha <- params$alpha
    beta <- params$beta

    # Initialize values
    geno_freq <- geno_freq_monogenic(p_alleles = c(1 - af, af))
    trans <- trans_monogenic(n_alleles = 2)
    baselineRisk <- calculateBaseline(
        cancer_type = cancer_type,
        gene = "SEER", race = "All_Races", type = "Net", db
    )

    # Calculate penetrance
    penet <- t(sapply(1:nrow(families), function(i) {
        penet.fn(i, families, alpha, beta, delta, gamma, max_age, baselineRisk, homozygote = homozygote, SeerNC = SeerNC, sex = sex)
    }))

    # Compute log-likelihood
    loglik <- pedigree_loglikelihood(families, geno_freq, trans, penet, ncores = 1)

    # Handle -Inf values
    if (is.infinite(loglik) && loglik == -Inf) {
        penalty <- 1e-28
        loglik <- log(penalty)
    } else {
        cat("Parameters:", given_median, given_first_quartile, alpha, beta, delta, gamma, "\n")
        cat("Log Likelihood:", loglik, "\n")
    }

    return(loglik)
}
