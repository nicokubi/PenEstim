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
calculate_lifetime_risk <- function(cancer_type, gene, race = "All_Races", sex = "Female", type = "Crude", data) {
    # Find the indices for the respective attributes
    dim_names <- attr(data$Penetrance, "dimnames")
    gene_index <- which(dim_names$Gene == gene)
    cancer_index <- which(dim_names$Cancer == cancer_type)
    race_index <- which(dim_names$Race == race)
    sex_index <- which(dim_names$Sex == sex)
    type_index <- which(dim_names$PenetType == type)

    # Calculate the cumulative risk
    lifetime_risk <- data$Penetrance[cancer_index, gene_index, race_index, sex_index, , type_index]
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
calculate_weibull_parameters <-
    function(given_median, given_first_quartile, shift, asymptote) {
        # Calculate alpha
        alpha <- log(log(4 / 3) / log(2)) /
            log((given_first_quartile - shift) / (given_median - shift))

        # Calculate beta using the median (M)
        beta <- (given_median - shift) / (log(2)^(1 / alpha))

        return(list(alpha = alpha, beta = beta))
    }

validate_weibull_parameters <-
    function(given_first_quartile, given_median, shift, asymptote) {
        # Check for negative or zero values
        if (given_median <= 0 || given_first_quartile <= 0 || shift < 0) {
            return(FALSE)
        }

        # Check if asymptote (gamma) is within the valid range (0,1)
        if (asymptote <= 0 || asymptote >= 1) {
            return(FALSE)
        }

        # Check if the logarithmic calculations will be valid
        if (given_first_quartile <= shift || given_median <= shift) {
            return(FALSE)
        }

        # Check if the denominator in the alpha calculation would be zero
        if ((given_first_quartile - shift) == (given_median - shift)) {
            return(FALSE)
        }

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
prepAges <- function(data, removeProband) {
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


na_rows <- is.na(simFamilies_C_1000_nocen_selected[[1]][["AgeBC"]])
simFamilies_C_1000_nocen_selected[[2]][["AgeBC"]][na_rows]
na_rows
