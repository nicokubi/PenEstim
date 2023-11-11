# Helper function to calculate lifetime risk
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
    lifetime_risk_cum <- cumsum(lifetime_risk)
    total_prob <- sum(lifetime_risk)

    return(list(cumulative_risk = lifetime_risk_cum, total_probability = total_prob))
}
calculate_lifetime_risk("SEER","Breast")

