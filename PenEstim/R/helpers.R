#' Calculate Lifetime Risk of Cancer
#'
#' This function calculates the lifetime risk of a specific type of cancer based
#' on genetic and demographic factors using data in the format of the PanelPRO Database.
#'
#' @param cancer_type The type of cancer for which the risk is being calculated.
#' @param gene The gene of interest for which the risk is being calculated.
#' @param race The race of the individual.
#' @param sex The sex of the individual, Female or Male.
#' @param type The type of penetrance calculation, Net or Crude.
#' @param db The dataset used for the calculation.
#'
#' @return A list containing the cumulative risk (`cumulative_risk`) and the total
#' probability (`total_prob`) of developing the specified cancer.
#'
#' @export
calculate_lifetime_risk <- function(cancer_type, gene, race, sex, type, db) {
  # Find the indices for the respective attributes
  dim_names <- attr(db$Penetrance, "dimnames")
  gene_index <- which(dim_names$Gene == gene)
  cancer_index <- which(dim_names$Cancer == cancer_type)
  race_index <- which(dim_names$Race == race)
  type_index <- which(dim_names$PenetType == type)
  
  cumulative_risk <- list(
    female = numeric(94),
    male = numeric(94),
    joint = numeric(93)
  )
  
  lifetime_risk <- list(
    female = numeric(94),
    male = numeric(94),
    joint = numeric(93)
  )
  
  # Calculate the cumulative risk
  risk <- db$Penetrance[cancer_index, gene_index, race_index, , , type_index]
  cumulative_risk$female <- cumsum(risk[1, ])
  cumulative_risk$male <- cumsum(risk[2, ])
  cumulative_risk$joint <- cumsum(colMeans(risk))
  
  lifetime_risk$female <- sum(risk[1, ])
  lifetime_risk$male <- sum(risk[2, ])
  lifetime_risk$joint <- sum(colMeans(risk))
  
  return(list(risk = risk, cumulative_risk = cumulative_risk, lifetime_risk = lifetime_risk))
}

#' Calculate Weibull Parameters
#'
#' This function calculates the shape (\code{alpha}) and scale (\code{beta}) parameters
#' of a Weibull distribution given the median, first quartile, and delta values.
#'
#' @param given_median The median of the data.
#' @param given_first_quartile The first quartile of the data.
#' @param delta A constant offset value.
#'
#' @return A list containing the calculated Weibull parameters, \code{alpha} and \code{beta}.
#
#' @export
calculate_weibull_parameters <- function(given_median, given_first_quartile, delta) {
  # Calculate alpha
  alpha <- log(-log(0.5) / -log(0.75)) / log((given_median - delta) / (given_first_quartile - delta))
  
  # Calculate beta using the median (M)
  beta <- (given_median - delta) / (-log(0.5))^(1 / alpha)
  
  return(list(alpha = alpha, beta = beta))
}

#' Validate Weibull Parameters
#'
#' This function validates the given parameters for calculating Weibull distribution.
#'
#' @param given_first_quartile The first quartile of the data.
#' @param given_median The median of the data.
#' @param threshold A constant threshold value.
#' @param asymptote A constant asymptote value (gamma).
#'
#' @return Boolean indicating whether the parameters are valid (TRUE) or not (FALSE).
#'
#' @export
validate_weibull_parameters <- function(given_first_quartile, given_median, threshold, asymptote) {
  # Check for negative or zero values
  if (given_median <= 0 || given_first_quartile <= 0 || threshold < 0) {
    return(FALSE)
  }
  
  # Check if asymptote (gamma) is within the valid range (0,1)
  if (asymptote <= 0 || asymptote >= 1) {
    return(FALSE)
  }
  
  # Check if the logarithmic calculations will be valid
  if (given_first_quartile <= threshold || given_median <= threshold) {
    return(FALSE)
  }
  
  # Check if the denominator in the alpha calculation would be zero
  if ((given_first_quartile - threshold) == (given_median - threshold)) {
    return(FALSE)
  }
  
  # If all checks pass, return TRUE
  return(TRUE)
}

#' Transform Data Frame
#'
#' This function transforms a data frame from the standard format used in PanelPRO
#' into the required format which conforms to the requirements of PenEstim (and clipp).
#'
#' @param df The input data frame in the usual PanelPRO format.
#'
#' @return The transformed data frame in the format required for clipp.
#'
#' @export
transformDF <- function(df) {
  # Rename and transform columns
  df$individual <- df$ID
  df$isProband <- df$isProband
  df$family <- df$PedigreeID
  df$mother <- df$MotherID
  df$father <- df$FatherID
  df$aff <- df$isAff
  df$sex <- ifelse(df$Sex == 0, 2, df$Sex) # Convert 0s to 2s in sex, keep 1s as is
  
  # Apply row-wise logic to assign age based on aff column
  df$age <- apply(df, 1, function(row) {
    aff <- as.numeric(row[["isAff"]])
    if (aff == 1) {
      return(as.numeric(row[["Age"]]))
    } else {
      return(as.numeric(row[["CurAge"]]))
    }
  })
  
  # Process 'geno' column
  df$geno <- ifelse(is.na(df$geno), "", ifelse(df$geno == 1, "1/2", ifelse(df$geno == 0, "1/1", df$geno)))
  
  # Select only the necessary columns
  df <- df[c("individual", "isProband", "family", "mother", "father", "aff", "sex", "age", "geno")]
  
  return(df)
}
