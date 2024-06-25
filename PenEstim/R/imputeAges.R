#' Impute Ages Based on Affection Status and Sex
#'
#' This function imputes ages for individuals in a dataset based on their affection
#' status and sex using either Weibull, baseline, or empirical distribution.
#'
#' @param data A data frame containing the data.
#' @param na_indices A vector of indices indicating the rows in the data where ages need to be imputed.
#' @param baseline_male A data frame containing baseline data for males.
#' @param baseline_female A data frame containing baseline data for females.
#' @param alpha_male Shape parameter for the Weibull distribution for males.
#' @param beta_male Scale parameter for the Weibull distribution for males.
#' @param delta_male Location parameter for the Weibull distribution for males.
#' @param alpha_female Shape parameter for the Weibull distribution for females.
#' @param beta_female Scale parameter for the Weibull distribution for females.
#' @param delta_female Location parameter for the Weibull distribution for females.
#' @param empirical_density A density object containing the empirical density of ages.
#' @return The data frame with imputed ages.
#' @export
imputeAges <- function(data, na_indices, baseline_male, baseline_female, alpha_male, beta_male, delta_male,
                       alpha_female, beta_female, delta_female, empirical_density) {
  for (i in na_indices) {
    if (is.na(data$age[i])) {
      u <- runif(1)
      relationship_prob <- as.numeric(data$degree_of_relationship[i])

      if (data$aff[i] == 1) {
        # Use Weibull distribution for carriers
        if (runif(1) < relationship_prob) {
          if (data$sex[i] == 1) { # Male
            age <- delta_male + beta_male * (-log(1 - u))^(1 / alpha_male)
          } else if (data$sex[i] == 2) { # Female
            age <- delta_female + beta_female * (-log(1 - u))^(1 / alpha_female)
          }
        } else {
          # Use baseline distribution for non-carriers
          age <- ifelse(data$sex[i] == 1, drawBaseline(baseline_male), drawBaseline(baseline_female))
        }
      } else {
        # Use empirical distribution for unaffected individuals
        age <- drawEmpirical(empirical_density)
      }
      data$age[i] <- max(1, round(age))
    }
  }
  return(data)
}

#' Initialize Ages Using a Uniform Distribution
#'
#' This function initializes ages for individuals in a dataset using a uniform distribution.
#'
#' @param data A data frame containing the data.
#' @param threshold The minimum age value for the uniform distribution.
#' @param max_age The maximum age value for the uniform distribution.
#' @return A list containing the updated data frame and the indices of the rows where ages were initialized.
#' @export
imputeAgesInit <- function(data, threshold, max_age) {
  na_indices <- which(is.na(data$age))
  data$age[na_indices] <- runif(length(na_indices), threshold, max_age)
  return(list(data = data, na_indices = na_indices))
}

#' Calculate Degree of Relationship in Pedigree
#'
#' This function calculates the degree of relationship for individuals in a dataset based on their pedigree.
#'
#' @param data A data frame containing the data with columns 'individual', 'father', 'mother', 'sex', 'aff', 'family', and 'isProband'.
#' @return The data frame with an additional column 'degree_of_relationship' indicating the degree of relationship for each individual.
#' @export
calcPedDegree <- function(data) {
  # Create a copy of the data to avoid modifying the original data directly
  data_copy <- data

  # Replace 0 with NA in the mother and father columns in the copy
  data_copy$mother[data_copy$mother == 0] <- NA
  data_copy$father[data_copy$father == 0] <- NA

  # Ensure both parents are NA if one is missing
  data_copy$father[is.na(data_copy$father) != is.na(data_copy$mother)] <- NA
  data_copy$mother[is.na(data_copy$father) != is.na(data_copy$mother)] <- NA

  # Create the pedigree object
  ped <- pedigree(
    id = data_copy$individual,
    dadid = data_copy$father,
    momid = data_copy$mother,
    sex = data_copy$sex,
    affected = data_copy$aff,
    famid = data_copy$family
  )

  # Calculate the kinship matrix
  kin_matrix <- kinship(ped)

  # Function to calculate the degree of relationship
  calculate_degree <- function(proband_id, family_kin_matrix) {
    degrees <- rep(NA, nrow(family_kin_matrix))
    names(degrees) <- rownames(family_kin_matrix)
    for (i in 1:nrow(family_kin_matrix)) {
      if (i == proband_id) {
        degrees[i] <- 0
      } else {
        kin_value <- family_kin_matrix[proband_id, i]
        degrees[i] <- kin_value * 2
      }
    }
    return(degrees)
  }

  # Initialize a column for degrees of relationship in the copy
  data_copy$degree_of_relationship <- NA

  # Iterate through each family to calculate the degree of relationship
  families <- unique(data_copy$family)
  for (fam in families) {
    family_data <- data_copy[data_copy$family == fam, ]
    proband_actual_id <- family_data$individual[family_data$isProband == 1]

    if (length(proband_actual_id) == 1) {
      # Subset the kinship matrix for the current family
      family_ids <- family_data$individual
      family_kin_matrix <- kin_matrix[family_ids, family_ids]

      # Find the row index of the proband in the family kinship matrix
      proband_index <- which(family_ids == proband_actual_id)

      # Calculate degrees of relationship for this family
      degrees_of_relationship <- calculate_degree(proband_index, family_kin_matrix)

      # Update the copy of the main data frame
      data_copy$degree_of_relationship[data_copy$family == fam] <- degrees_of_relationship
    }
  }

  # Add the new column to the original data
  data$degree_of_relationship <- data_copy$degree_of_relationship

  return(data)
}

#' Calculate Empirical Density for Non-Affected Individuals
#'
#' This function calculates the empirical density for ages of non-affected individuals in a dataset.
#'
#' @param data A data frame containing the data.
#' @param aff_column The name of the column indicating affection status (default is "aff").
#' @param age_column The name of the column indicating ages (default is "age").
#' @param n_points The number of points to use in the density estimation (default is 10000).
#' @return A density object representing the empirical density of ages.
#' @export
calculateEmpiricalDensity <- function(data, aff_column = "aff", age_column = "age", n_points = 10000) {
  # Filter the data to include only non-affected individuals (aff == 0)
  non_affected_data <- subset(data, data[[aff_column]] == 0)

  # Remove NA values from the age column of the filtered data
  cleaned_non_affected_ages <- na.omit(non_affected_data[[age_column]])

  # Estimate the empirical density of the age data
  age_density <- density(cleaned_non_affected_ages, n = n_points)

  return(age_density)
}

#' Draw Ages Using the Inverse CDF Method from the baseline data
#'
#' This function draws ages using the inverse CDF method from baseline data.
#'
#' @param baseline_data A data frame containing baseline data with columns 'cum_prob' and 'age'.
#' @return A single age value drawn from the baseline data.
#' @export
drawBaseline <- function(baseline_data) {
  u <- runif(1)
  age <- approx(baseline_data$cum_prob, baseline_data$age, xout = u)$y
  return(age)
}

#' Draw Ages Using the Inverse CDF Method from Empirical Density
#'
#' This function draws ages using the inverse CDF method from empirical density data.
#'
#' @param empirical_density A density object containing the empirical density of ages.
#' @return A single age value drawn from the empirical density data.
#' @export
drawEmpirical <- function(empirical_density) {
  u <- runif(1)
  age <- approx(cumsum(empirical_density$y) / sum(empirical_density$y), empirical_density$x, xout = u)$y
  return(age)
}
