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

calculateNCPen <- function(SEER_baseline, carrierPenetrances, af, homozygous = FALSE) {
    # Calculate probability weights for carriers based on allele frequencies
    weights1 <- 2 * af * (1 - af) # Heterozygous carriers
    weights2 <- af^2 # Homozygous carriers

    # Conditionally include weights for homozygous carriers
    weights <- if (homozygous) {
        weights1 + weights2
    } else {
        weights1
    }
    # Assuming carrierPenetrances is structured to align with the age and possibly sex
    # This might require adjustment based on your data
    weightedCarrierPenetrances <- carrierPenetrances * weights
    weightedSumCarrierPenetrances <- sum(weightedCarrierPenetrances)

    # Calculate age-specific non-carrier penetrance
    nonCarrierPenetrance <- (SEER_baseline - weightedSumCarrierPenetrances) / (1 - sum(weights))

    # Ensure non-negative value
    nonCarrierPenetrance <- max(0, nonCarrierPenetrance)

    return(nonCarrierPenetrance)
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
    if (is.na(data$sex[i])) {
        # Handle NA in sex - using average risk
        sex_indices <- 1:nrow(baselineRisk)
    } else {
        # Map sex to row index: Assuming "Female" is 1st row and "Male" is 2nd row
        sex_index <- ifelse(data$sex[i] == 2, 1, 2)
    }

    if (data$age[i] == 0) {
        penet.i <- c(1, 1, 1) # Assuming people aged 0 years are all unaffected
    } else {
        # Ensure age is within the valid range
        age_index <- min(max_age, data$age[i])

        # Weibull hazard and survival calculation differs for males and females
        if (data$sex[i] == 2) { # Assuming "Male" is coded as 2
            c.pen <- dweibull(data$age[i] - delta, alpha, beta) * gamma
        } else { # Female
            c.pen <- dweibull(data$age[i] - delta, alpha, beta) * gamma
        }

        # Extract the corresponding baseline risk for sex and age
        SEER_baseline_i <- baselineRisk[sex_index, age_index]

        # If Assuming the non-carrier risk is equal to the SEER baseline risk
        if (SeerNC == TRUE) {
            nc.pen <- SEER_baseline_i
        } else {
            nc.pen <- calculateNCPen(SEER_baseline = SEER_baseline_i, carrierPenetrances = c.pen, af = af)
        }

        # Penetrance calculations based on genotype
        # Assuming same penetrance for homozygous and heterozygous carriers
        penet.i <- c(1 - nc.pen, 1 - c.pen, 1 - c.pen) # if person is not affected
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
    browser()
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

    # Compute log-likelihooda
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
