#' Draw Initial Parameters Based on Kaplan-Meier Statistics
#'
#' This function calculates Kaplan-Meier statistics for male and female affected individuals in a dataset. 
#' It filters the data by sex and affected status, and computes various statistics such as the median age,
#' first quartile age, and asymptotic survival probability.
#'
#' @param data A data frame containing the columns 'sex', 'aff', 'geno', and 'age'. 
#'   - `sex` should be a numeric code where 1 indicates male and 2 indicates female.
#'   - `aff` should indicate affected status where 1 represents affected.
#'   - `geno` should represent genotype and the function filters for "1/2".
#'   - `age` should be the age at which the event (affected status) is observed.
#'
#' @return A list with the Kaplan-Meier statistics for both male and female:
#'   - `asymptote_male` and `asymptote_female`: The asymptotic survival probability estimates.
#'   - `threshold_male` and `threshold_female`: The minimum age at which an event is observed.
#'   - `median_age_male` and `median_age_female`: The median ages at event occurrence.
#'   - `first_quartile_age_male` and `first_quartile_age_female`: The age by the first quartile of the survival function.
#'
#' @examples
#' # Assuming 'data' is a dataframe with appropriate structure:
#' results <- draw_initial_params_v2(data)
#'
#' @importFrom survival survfit
#' @importFrom survival Surv
#' @export
draw_initial_params_v2 <- function(data) {
  # Filter data by sex and affected status
  data_male_affected <- data[data$sex == 1 & data$aff == 1 & data$geno == "1/2",]
  data_female_affected <- data[data$sex == 2 & data$aff == 1 & data$geno == "1/2",]
  
  # Helper function to compute Kaplan-Meier statistics
  compute_km_stats <- function(subset) {
    if (nrow(subset) == 0) {
      return(list(threshold = NA, median_age = NA, first_quartile_age = NA, asymptote = NA))
    }
    
    km_fit <- survfit(Surv(subset$age, subset$aff) ~ subset$geno, data = subset)
    km_summary <- summary(km_fit)
    
    # Threshold (minimum age at which an event is observed)
    threshold <- min(km_fit$time)
    
    # First quartile age and median age
    surv_probs <- km_summary$surv
    times <- km_summary$time
    first_quartile_age <- times[which.min(abs(surv_probs - 0.75))]
    median_age <- times[which.min(abs(surv_probs - 0.50))]
    
    # Attempt to find a 'flattening' of the survival curve for asymptote estimation
    last_surv_prob <- surv_probs[length(surv_probs)-1]
    asymptote <- 1 - last_surv_prob
    
    return(list(threshold = threshold, median_age = median_age, first_quartile_age = first_quartile_age, asymptote = asymptote))
  }
  
  # Compute KM stats for male and female subsets
  stats_male <- compute_km_stats(data_male_affected)
  stats_female <- compute_km_stats(data_female_affected)
  
  return(list(
    asymptote_male = stats_male$asymptote,
    asymptote_female = stats_female$asymptote,
    threshold_male = stats_male$threshold,
    threshold_female = stats_female$threshold,
    median_age_male = stats_male$median_age,
    median_age_female = stats_female$median_age,
    first_quartile_age_male = stats_male$first_quartile_age,
    first_quartile_age_female = stats_female$first_quartile_age
  ))
}

