create_distributions <- function(dataframe, samples_table) {
  
  # Extract values from dataframe
   # Extract values from dataframe
  max_age <- max(dataframe$age)
  min_age <- min(dataframe$age)
  
  # Median penetrance time
  median_age <- dataframe$age[which.min(dataframe$penetrance_prob <= 0.5)]
  first_quartile <- quantile(dataframe$age, 0.25)
  
  print(paste("Max Age:", max_age))
  print(paste("Min Age:", min_age))
  print(paste("Median Age:", median_age))
  print(paste("First Quartile:", first_quartile))
  # Define the scaling transformation
  normalize <- function(x) {
    return((x - min_age) / (max_age - min_age))
  }
  
  denormalize <- function(x) {
    return(min_age + (max_age - min_age) * x)
  }
  
  # Extract risk numbers from the samples_table
  risk_median <- samples_table$at_risk[samples_table$age == median_age]
  risk_quartile <- samples_table$at_risk[samples_table$age == first_quartile]
  
  print(paste("Risk for Median Age:", risk_median))
  print(paste("Risk for First Quartile:", risk_quartile))
  
  # Compute beta distribution parameters for median and first quartile
  compute_parameters <- function(stat, at_risk) {
    med <- normalize(stat)
    var <- med * (1 - med)/at_risk
    
    alpha <- med * ((med*(1-med)/var) - 1)
    beta <- (1-med) * ((med*(1-med)/var) - 1)
    
    return(list(alpha = alpha , beta = beta))
  }
  
  params_median <- compute_parameters(median_age,risk_median)
  params_quartile <- compute_parameters(first_quartile,risk_quartile)

  print(params_median)
  print(params_quartile)
  
  # Asymptote distribution
  asymptote_distribution <- function(n) {
    qbeta(runif(n), 2, 2, max_age-1, max_age + 1)
  }
  
  # Shift parameter distribution
  shift_distribution <- function(n) {
    runif(n, 0, min_age)
  }
  
  # Median distribution
  median_distribution <- function(n) {
    denormalize(qbeta(runif(n), params_median$alpha, params_median$beta))
  }
  
  # First quartile distribution
  first_quartile_distribution <- function(n) {
    denormalize(qbeta(runif(n), params_quartile$alpha, params_quartile$beta))
  }
  
  # Return a list of the distributions
  list(
    asymptote_distribution = asymptote_distribution,
    shift_distribution = shift_distribution,
    median_distribution = median_distribution,
    first_quartile_distribution = first_quartile_distribution
  )
}

