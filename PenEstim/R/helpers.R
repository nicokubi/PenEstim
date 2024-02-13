#' Calculate Lifetime Risk of Cancer
#'
#' This function calculates the lifetime risk of a specific type of cancer based
#' on genetic and demographic factors using data from the PanelPRO Database.
#'
#' @param cancer_type The type of cancer for which the risk is being calculated.
#' @param gene The gene of interest for which the risk is being calculated.
#' @param race The race of the individual, with a default value of "All_Races".
#' @param sex The sex of the individual, with a default value of "Female".
#' @param type The type of calculation, with a default value of "Crude".
#' @param data The dataset used for the calculation, defaulting to PanelPRODatabase.
#'
#' @return A list containing the cumulative risk (`cumulative_risk`) and the total
#' probability (`total_prob`) of developing the specified cancer.
#'
#' @examples
#' calculate_lifetime_risk("Breast_Cancer", "BRCA1")
#'
#' @export
calculate_lifetime_risk <- function(cancer_type, gene, race = "All_Races", sex = "Female", type = "Crude", db) {
    # Find the indices for the respective attributes
    dim_names <- attr(db$Penetrance, "dimnames")
    gene_index <- which(dim_names$Gene == gene)
    cancer_index <- which(dim_names$Cancer == cancer_type)
    race_index <- which(dim_names$Race == race)
    sex_index <- which(dim_names$Sex == sex)
    type_index <- which(dim_names$PenetType == type)

    # Calculate the cumulative risk
    lifetime_risk <- db$Penetrance[cancer_index, gene_index, race_index, sex_index, , type_index]
    cumulative_risk <- cumsum(lifetime_risk)
    total_prob <- sum(lifetime_risk)

    return(list(cumulative_risk = cumulative_risk, total_prob = total_prob))
}

# Function to generate proposals
generate_proposal <- function(distribution_func, args_list) {
    # Use do.call to dynamically call the distribution function with arguments
    proposal <- do.call(distribution_func, args_list)
    return(proposal)
}

#' Calculate Weibull Parameters
#'
#' This function calculates the shape (\code{alpha}) and scale (\code{beta}) parameters
#' of a Weibull distribution given the median, first quartile, delta, and gamma values.
#'
#' @param given_median The median of the data.
#' @param given_first_quartile The first quartile of the data.
#' @param delta A constant offset value.
#' @param gamma A constant shape parameter.
#'
#' @return A list containing the calculated Weibull parameters, \code{alpha} and \code{beta}.
#'
#' @examples
#' # Example usage:
#' result <- calculate_weibull_parameters(50, 25, 0.1, 2)
#' cat("Alpha:", result$alpha, "\n")
#' cat("Beta:", result$beta, "\n")
#'
#' @export

# Ensure the Q function is defined somewhere in your script
quantile.fn <- function(p, beta, alpha, shift) {
    beta * (-log(1 - p))^(1 / alpha) + shift
}

calculate_weibull_parameters <- function(given_median, given_first_quartile, shift, asymptote) {
    # Objective function to minimize
    objective_function <- function(params, given_median, given_first_quartile, shift, asymptote) {
        beta <- params[1] # Scale parameter
        alpha <- params[2] # Shape parameter

        # Calculate the predicted quantiles based on current beta and alpha
        median_pred <- quantile.fn(0.5, beta, alpha, shift)
        first_quartile_pred <- quantile.fn(0.25, beta, alpha, shift)

        # Scale the given parameters
        given_median <- given_median / asymptote
        given_first_quartile <- given_first_quartile / asymptote

        # Sum of squared differences
        sum_of_squares <- (given_median - median_pred)^2 + (given_first_quartile - first_quartile_pred)^2

        return(sum_of_squares)
    }

    # Initial guesses for beta and alpha
    initial_guesses <- c(beta = 5, alpha = 1)

    # Use optim to minimize the objective function
    result <- optim(initial_guesses, objective_function,
        given_median = given_median,
        given_first_quartile = given_first_quartile,
        asymptote = asymptote,
        shift = shift
    ) # Adding method and lower bounds to ensure positive parameters

    return(list(beta = result$par[1], alpha = result$par[2]))
}

calculate_weibull_parameters_vectorized <- function(given_median, given_first_quartile, shift, asymptote) {
    # Initialize lists to store the results
    betas <- numeric(length(given_median))
    alphas <- numeric(length(given_median))

    # Loop over each set of inputs
    for (i in seq_along(given_median)) {
        # Single set of inputs
        single_given_median <- given_median[i]
        single_given_first_quartile <- given_first_quartile[i]
        single_shift <- shift[i]
        single_asymptote <- asymptote[i]

        # Optimization for a single set of inputs
        result <- calculate_weibull_parameters(single_given_median, single_given_first_quartile, single_shift, single_asymptote)

        # Store the results
        betas[i] <- result$beta
        alphas[i] <- result$alpha
    }

    # Return the results as a list of vectors
    return(list(beta = betas, alpha = alphas))
}

# Example usage (you need to define given_median, given_first_quartile, shift, and asymptote)
# calculate_weibull_parameters(given_median, given_first_quartile, shift, asymptote)

validate_weibull_parameters <-
    function(given_first_quartile, given_median, shift, asymptote) {
        # Check for negative or zero values

        # If all checks pass, return TRUE
        return(TRUE)
    }

#' Prepare Ages from Pedigree Data
#'
#' This function processes a list of data frames containing information about various cancer types.
#' For each cancer type, it adjusts the age column and the corresponding affected-status column.
#' It also computes the maximum cancer age for each row and updates the 'CurAge' field accordingly.
#' This function is designed to handle multiple cancer types and modify age-related fields in the data.
#'
#' @param data A list of data frames, where each data frame contains columns for age and affected-status
#'             for different cancer types. The function expects specific column naming conventions like 'AgeBRA',
#'             'isAffBRA' for each cancer type.
#' @return The modified list of data frames, with updated age and affected-status columns for each cancer type,
#'         and updated 'CurAge' values.
#' @examples
#' # Example usage:
#' data <- list(df1, df2) # df1, df2 are data frames with appropriate columns
#' result <- prepAges(data)
#'
prepAges <- function(data, removeProband = FALSE) {
    # Define the list of genes
    genes <- c(
        "APC", "ATM", "BARD1", "BMPR1A", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CDK4",
        "CDKN2A", "CHEK2", "EPCAM", "MLH1", "MSH2", "MSH6", "MUTYH", "NBN", "PALB2",
        "PMS2", "PTEN", "RAD51C", "RAD51D", "STK11", "TP53"
    )

    # Add all cancer types
    cancer_types <- c(
        "BRA", "BC", "CER", "COL", "ENDO", "GAS", "KID", "LEUK", "MELA",
        "OC", "OST", "PANC", "PROS", "SMA", "STS", "THY", "UB", "HEP", "CBC", "BC2"
    )

    for (i in seq_along(data)) {
        # If removeProband is TRUE, set genes, ages, and affections to NA for probands
        if (removeProband) {
            proband_rows <- data[[i]]$Proband == 1
            for (gene in genes) {
                if (gene %in% colnames(data[[i]])) {
                    data[[i]][proband_rows, gene] <- NA
                }
            }
            for (cancer in cancer_types) {
                age_col <- paste0("Age", cancer)
                aff_col <- paste0("isAff", cancer)
                if (age_col %in% colnames(data[[i]])) {
                    data[[i]][proband_rows, age_col] <- NA
                }
                if (aff_col %in% colnames(data[[i]])) {
                    data[[i]][proband_rows, aff_col] <- NA
                }
            }
        }

        for (cancer in cancer_types) {
            age_col <- paste0("Age", cancer)
            aff_col <- paste0("isAff", cancer)

            if (age_col %in% colnames(data[[i]])) {
                na_rows <- is.na(data[[i]][[age_col]])
                data[[i]][[age_col]][na_rows] <- 1
                data[[i]][[aff_col]][na_rows] <- 0
            }
        }

        # Calculate max cancer age only for columns that exist in the data
        cancer_ages_cols <- paste0("Age", cancer_types)
        cancer_ages_cols <- intersect(cancer_ages_cols, colnames(data[[i]]))
        cancer_ages <- data[[i]][cancer_ages_cols]

        if (length(cancer_ages_cols) > 0) {
            max_cancer_age <- apply(cancer_ages, 1, max, na.rm = TRUE)
            data[[i]]$CurAge <- ifelse(is.na(data[[i]]$CurAge), max_cancer_age, pmax(data[[i]]$CurAge, max_cancer_age))
        }
    }
    return(data)
}
