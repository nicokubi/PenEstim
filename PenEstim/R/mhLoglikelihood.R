#' Log-Likelihood Calculation for Metropolis-Hastings
#'
#' @param paras Vector of parameters for which likelihood is to be computed.
#' @param data List of families data.
#' @return Log-likelihood value.
#' @importFrom PPP PPP


mhLogLikelihood <- function(paras, families, max_age, PanelPRODatabase) {
  # set age, same as in DB
  age <- seq(1, max_age, 1)
  # use same BRCA1 frequency (not estimated)
  BRCA1freq <- 0.9
  PanelPRODatabase$AlleleFrequency[["BRCA1_anyPV",1]] <- BRCA1freq
  PanelPRODatabase$AlleleFrequency[["BRCA1_anyPV",2]] <- BRCA1freq
  PanelPRODatabase$AlleleFrequency[["BRCA1_anyPV",3]] <- BRCA1freq
  
  # Parameters, drawn from prior/proposal
  given_median <- paras[1]
  given_first_quartile <- paras[2]
  gamma <- paras[3]
  delta <- paras[4]
  
  calculate_weibull_parameters <- function(given_median, given_first_quartile, delta, gamma) {
    # Calculate alpha
    alpha <- log(-log((gamma-0.25)/gamma) / -log((gamma-0.5)/gamma)) /
      log((given_first_quartile - delta) / (given_median - delta))
    
    # Calculate beta using the median (M)
    beta <- (given_median - delta) / (log(2)^(1 / alpha))
    
    return(list(alpha = alpha, beta = beta))
  }
  
  params <- calculate_weibull_parameters(given_median, given_first_quartile, delta,gamma)
  alpha <- params$alpha
  beta <- params$beta
  
  # Now use alpha and beta in your simulation
  penetrance.mod.f <- dweibull(age - delta, alpha, beta) * gamma
  
  # For now focus on just one vector of penetrance estimates
  gene <- "BRCA1_hetero_anyPV"
  cancer <- "Breast"
  race <- "All_Races"
  female <- "Female"
  male <- "Male"
  type <- "Net"
  
  # Find the indices for the resp. attributes
  dim_names <- attr(PanelPRODatabase$Penetrance, "dimnames")
  gene_index <- which(dim_names$Gene == gene)
  cancer_index <- which(dim_names$Cancer == cancer)
  race_index <- which(dim_names$Race == race)
  sex_index <- which(dim_names$Sex == female)
  type_index <- which(dim_names$PenetType == type)
  
  
  # overwrite the penetrance in the PanelPro Database
  PanelPRODatabase$Penetrance[cancer_index, gene_index, race_index, sex_index, ,] <- penetrance.mod.f
  
  # Male Penetrance
  sex_index <- which(dim_names$Sex == male)
  
  # Male penetrance funciton
  penetrance.mod.m <- 0
  # overwrite the penetrance in the PanelPro Database
  PanelPRODatabase$Penetrance[cancer_index, gene_index, race_index, sex_index, ,] <- penetrance.mod.m
  
  # Storing the estimates
  log_likelihood <- 0
  
  for (i in 1:length(data)) {
    data <- families[[i]]  # Get the data for the current family
    
    # Access the posterior probabilities (not normalized) and estimates for BRCA1 gene in Breast cancer
    postprobs <- PPP(pedigree=data, genes = c("BRCA1"), cancers = "Breast", database= PanelPRODatabase)$posterior.prob[[1]]  # check bug
    estimate <- postprobs[postprobs$genes=="BRCA1_hetero_anyPV","estimate"]
    
    if (is.nan(estimate) || estimate <=0) {
      # Handle NaN or zero probabilities by adding a small value and/or penalizing
      estimate <- 1e-28
      ll <- log(estimate)
      log_likelihood <- log_likelihood + ll
      cat("NaN or zero encountered:", alpha, beta, gamma, delta, ll, log_likelihood, estimate, "\n")
    } else {
      ll <- log(estimate)
      log_likelihood <- log_likelihood + ll
      cat(given_median, given_first_quartile, alpha, beta, gamma, delta, ll, log_likelihood, estimate, "\n")
    }
  }
  
  # return the log-likelihood
  return(log_likelihood)
}
