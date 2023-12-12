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
calculate_lifetime_risk <- function(cancer_type, gene, race = "All_Races", sex = "Female", type = "Crude", data = PanelPRODatabase) {
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
calculate_weibull_parameters <- function(given_median, given_first_quartile, delta, gamma) {
    # Calculate alpha
    alpha <- log(-log((gamma - 0.25) / gamma) / -log((gamma - 0.5) / gamma)) /
        log((given_first_quartile - delta) / (given_median - delta))

    # Calculate beta using the median (M)
    beta <- (given_median - delta) / (log(2)^(1 / alpha))

    return(list(alpha = alpha, beta = beta))
}
