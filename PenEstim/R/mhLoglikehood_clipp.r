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

 calculateBaseline <- function(cancer_type, gene, race, type, data) {
     # Check if dimnames are available and correct
     if (is.null(data$Penetrance) || is.null(attr(data$Penetrance, "dimnames"))) {
         stop("Penetrance data or its dimension names are not properly defined.")
     }

     dim_names <- attr(data$Penetrance, "dimnames")
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
     lifetime_risk <- data$Penetrance[cancer_index, gene_index, race_index, , , type_index]
     return(lifetime_risk)
 }

#' Penetrance Function
#'
#' This function calculates the penetrance for an individual based on their age, genotype, and other factors.
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
 penet.fn <- function(i, data, alpha, beta, delta, gamma, max_age, baselineRisk) {
     if (is.na(data$sex[i])) {
         # Handle NA in sex - using average risk
         sex_indices <- 1:nrow(baselineRisk)
     } else {
         # Map sex to row index: Assuming "Female" is 1st row and "Male" is 2nd row
         sex_index <- ifelse(data$sex[i] == "Female", 1, 2)
         sex_indices <- sex_index
     }

     if (data$age[i] == 0) {
         penet.i <- c(1, 1) # Assuming people aged 0 years are all unaffected
     } else {
         # Ensure age is within the valid range
         age_index <- min(max_age, data$age[i])

         # Extract the corresponding baseline risk for sex and age
         nc.pen <- mean(baselineRisk[sex_indices, age_index])

         # Weibull hazard and survival for carriers
         c.pen <- dweibull(data$age[i] - delta, alpha, beta) * gamma

         # Penetrance calculations based on genotype
         penet.i <- c(1 - nc.pen, 1 - c.pen) # if person is not affected
         if (data$aff[i] == 1) penet.i <- c(nc.pen, c.pen)
     }

     if (data$geno[i] == "1/1") penet.i[-1] <- 0
     if (data$geno[i] == "1/2") penet.i[-2] <- 0

     return(penet.i)
 }

#' Transform Data Frame
#'
#' This function transforms a data frame into the required format for further analysis.
#'
#' @param df The input data frame in the usual PanelPRO format.
#'
#' @return The transformed data frame in the format required for clipp.
#'
#' @examples
#' # Transform a data frame
#' transformed_df <- transformDF(input_df)
#'
 
 transformDF <- function(df) {
     df %>%
         select(
             individual = SubjectID,
             family = PedigreeID,
             mother = MotherID,
             father = FatherID,
             aff = isAffCOL,
             sex = Sex,
             age = CurAge,
             geno = MSH6
         ) %>%
         mutate(
             geno = ifelse(is.na(geno), "", ifelse(geno == 1, "1/2", ifelse(geno == 0,"1/1",geno)))
         )
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
mhLogLikelihood_clipp <- function(paras, families, max_age, cancer_type, db, af) {

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
    geno_freq <- c((1 - af)^2, 2*(1-af)*af)
    trans <- trans_monogenic2(n_alleles = 2, nonviable = TRUE)
    baselineRisk <- calculateBaseline(cancer_type = cancer_type,
     gene = "SEER", race = "All_Races", type = "Crude", data = db)

    # Calculate penetrance
    penet <- t(sapply(1:nrow(families), function(i) {
        penet.fn(i, families, alpha, beta, delta, gamma, max_age, baselineRisk)
    }))

    # Compute log-likelihooda
    loglik <- pedigree_loglikelihood(families, geno_freq, trans, penet, ncores = 1)

    # Handle -Inf values
    if (is.infinite(loglik) && loglik == -Inf) {
        NEGATIVE_INFINITY_PENALTY <- -10000
        loglik <- NEGATIVE_INFINITY_PENALTY
        cat("Negative infinity encountered, penalized log likelihood:", loglik, "\n")
    } else {
        cat("Parameters:", given_median, given_first_quartile, alpha, beta, delta, gamma, "\n")
        cat("Log Likelihood:", loglik, "\n")
    }

    return(loglik)
}