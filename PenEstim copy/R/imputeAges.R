# Function to draw ages using the inverse CDF method from SEER data
drawSeer <- function(seer_data) {
  u <- runif(1)
  age <- approx(seer_data$cum_prob, seer_data$age, xout = u)$y
  return(age)
}

# Function to draw ages using the inverse CDF method from empirical density
drawEmpirical <- function(empirical_density) {
  u <- runif(1)
  age <- approx(cumsum(empirical_density$y) / sum(empirical_density$y), empirical_density$x, xout = u)$y
  return(age)
}

# Function to impute ages based on aff status and sex using either Weibull, SEER, or empirical distribution
imputeAges <- function(data, na_indices, SEER_male, SEER_female, alpha_male, beta_male, delta_male,
                       alpha_female, beta_female, delta_female, empirical_density) {
  
  for (i in na_indices) {
    if (is.na(data$age[i])) {
      u <- runif(1)
      relationship_prob <- as.numeric(data$degree_of_relationship[i])
      
      if (data$aff[i] == 1) {
        # Use Weibull distribution for carriers
        if (runif(1) < relationship_prob) {
          if (data$sex[i] == 1) {  # Male
            age <- delta_male + beta_male * (-log(1 - u))^(1 / alpha_male)
          } else if (data$sex[i] == 2) {  # Female
            age <- delta_female + beta_female * (-log(1 - u))^(1 / alpha_female)
          }
        } else {
          # Use SEER distribution for non-carriers
          age <- ifelse(data$sex[i] == 1, drawSeer(SEER_male), drawSeer(SEER_female))
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

# Function to initialize ages using a uniform distribution
imputeAgesInit <- function(data, threshold, max_age) {
  na_indices <- which(is.na(data$age))
  data$age[na_indices] <- runif(length(na_indices), threshold, max_age)
  return(list(data = data, na_indices = na_indices))
}

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
  ped <- pedigree(id = data_copy$individual,
                  dadid = data_copy$father,
                  momid = data_copy$mother,
                  sex = data_copy$sex,
                  affected = data_copy$aff,
                  famid = data_copy$family)
  
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

# Function to calculate empirical density for non-affected individuals
calculateEmpiricalDensity <- function(data, aff_column = "aff", age_column = "age", n_points = 10000) {
  # Filter the data to include only non-affected individuals (aff == 0)
  non_affected_data <- subset(data, data[[aff_column]] == 0)
  
  # Remove NA values from the age column of the filtered data
  cleaned_non_affected_ages <- na.omit(non_affected_data[[age_column]])
  
  # Estimate the empirical density of the age data
  age_density <- density(cleaned_non_affected_ages, n = n_points)
  
  return(age_density)
}
