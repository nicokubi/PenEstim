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
