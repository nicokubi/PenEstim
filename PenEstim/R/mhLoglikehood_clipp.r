#' Calculate Baseline Penetrance
#'
#' Calculates the baseline penetrance for a given cancer type, gene, race, sex, and penetrance type.
#' Extracts the necessary information from the provided database object.
#'
#' @param cancer_type String, the type of cancer.
#' @param gene String, the gene, with default "SEER".
#' @param race String, the race, with default "All_Races".
#' @param sex String, the sex, with default "Male".
#' @param type String, the type of penetrance, with default "Crude".
#' @param db Data frame or list, the data object containing penetrance information.
#'
#' @return Numeric matrix, the baseline penetrance values.
#'
#' @seealso \code{\link{penet.fn}}
#' @examples
#' # Calculate baseline penetrance
#' baseline <- calculateBaseline(
#'     cancer_type = "Lung Cancer", gene = "SEER",
#'     race = "All_Races", type = "Crude", db = my_data
#' )
#' @export
calculateBaseline <- function(cancer_type, gene = "SEER", race = "All_Races", sex = "Male", type = "Crude", db) {
    # Validations and index retrieval
    if (is.null(db$Penetrance) || is.null(dimnames(db$Penetrance))) {
        stop("Penetrance data or its dimension names are not properly defined.")
    }

    # Extracting indices for filtering
    indices <- lapply(list(Cancer = cancer_type, Gene = gene, Race = race, Sex = sex, PenetType = type), function(value, dim_name) {
        idx <- which(dimnames(db$Penetrance)[[dim_name]] == value)
        if (length(idx) == 0) stop(paste("Value", value, "not found in dimension", dim_name))
        idx
    }, dim_name = names(dimnames(db$Penetrance)))

    # Subsetting the Penetrance data for specified filters
    lifetime_risk <- db$Penetrance[indices$Cancer, indices$Gene, indices$Race, , indices$Sex, indices$PenetType]
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
#' @param gamma Numeric, adjustment factor for the penetrance.
#' @param max_age Integer, maximum age considered in the analysis.
#' @param baselineRisk Numeric vector, baseline risk for each age.
#' @param homozygote Logical, indicates if homozygous carriers should be considered.
#' @param SeerNC Logical, indicates if non-carrier penetrance should be based on SEER data.
#' @param sex Character, the sex of the individuals ('Male', 'Female', or 'NA').
#'
#' @return Numeric vector, containing penetrance values for unaffected and affected individuals.
#' @examples
#' # Example usage
#' baselineRisk <- calculateBaseline("Breast", "BRCA1", "All_Races", "Female", "Net", db)
#' individual_data <- data.frame(sex = c(2), age = c(30), aff = c(0), geno = c("1/2"))
#' penet.fn(1, individual_data, 1.0, 2.0, 0.5, 0.8, 80, baselineRisk, FALSE, TRUE, "Female")
#' @export
#'
# Define the penetrance function
penet.fn <- function(i, data, alpha, beta, delta, gamma, max_age, baselineRisk, homozygote, SeerNC, sex) {
    # Map sex to row index: "Female" is 1st row and "Male" is 2nd row
    sex_index <- ifelse(data$sex[i] == 2, 1, 2)

    # Early return for age 0 as fully unaffected
    if (data$age[i] == 0) {
        return(c(1, 1, 1)) # Unaffected for all genetic statuses
    }

    age_index <- min(data$age[i], max_age)
    baselineRiskAge <- baselineRisk[sex_index, age_index]

    # Weibull with four parameters for the penetrance
    # For now we just have one penetrance, irrespective of sex
    c.pen <- dweibull(data$age[i] - delta, alpha, beta) * gamma
    # Cumulative penetrance calculation using product of 1-year survival probabilities
    c.pen.c <- pweibull(data$age[i] - delta, alpha, beta)

    # Calculate cumulative risk for non-carriers based on SEER data or other model
    if (SeerNC == TRUE) {
        nc.pen <- baselineRiskAge
        # calculative cumulative product for being cancer-free for non-carriers
        nc.pen.c <- prod(1 - baselineRisk[sex_index, 1:age_index])
    } else {
        calculatedNCRisk <- calculateNCPen(baselineRisk[sex_index, ], alpha, beta, delta, gamma, data$af[i], age_index, homozygote)
        nc.pen <- calculatedNCRisk$yearlyProb[age_index]
        # calculative cumulative product for being cancer-free for non-carriers
        nc.pen.c <- calculatedNCRisk$cumulativeProb[age_index]
    }

    # Penetrance calculations based on genotype
    # Assuming same penetrance for homozygous and heterozygous carriers
    penet.i <- if (data$aff[i] == 1) {
        penet.i <- c(nc.pen, c.pen, c.pen)
    } else {
        c(nc.pen.c, 1 - c.pen.c, 1 - c.pen.c)
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
        penet.fn(i, families, alpha, beta, delta, gamma, max_age,
            baselineRisk,
            homozygote = homozygote, SeerNC = SeerNC, sex = sex
        )
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
