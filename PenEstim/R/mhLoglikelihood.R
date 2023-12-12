#' Log-Likelihood Calculation for Metropolis-Hastings
#'
#' @param paras Vector of parameters for which likelihood is to be computed.
#' @param families List of families data.
#' @param max_age Maximum age to be considered.
#' @param cancer_type The type of cancer for which to estimate penetrance.
#' @param gene_input The gene for which to estimate penetrance.
#' @param PanelPRODatabase List containing penetrance data for different cancer types.
#' @return Log-likelihood value.
#' @importFrom PPP PPP

mhLogLikelihood <- function(paras, families, max_age, PanelPRODatabase, cancer_type, gene_input) {
  # set age, same as in DB
  age <- seq(1, max_age, 1)

  # Parameters, drawn from prior/proposal
  given_median <- paras[1]
  given_first_quartile <- paras[2]
  gamma <- paras[3]
  delta <- paras[4]

  # Recalculate the parameters 
  params <- calculate_weibull_parameters(given_median, given_first_quartile, delta, gamma)
  alpha <- params$alpha
  beta <- params$beta

  # Now use alpha and beta in your simulation
  penetrance.mod.f <- dweibull(age - delta, alpha, beta) * gamma

  # For now focus on just one vector of penetrance estimates
  gene_adj <- paste(gene_input, "_hetero_anyPV", sep = "")
  gene <- gene_adj
  cancer <- cancer_type
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

  # overwrite the penetrance in the PanelPro Database
  PanelPRODatabase$Penetrance[cancer_index, gene_index, race_index, sex_index, , ] <- penetrance.mod.f

  #  Male Penetrance
  sex_index <- which(dim_names$Sex == male)

  # Male penetrance function
  penetrance.mod.m <- 0
  # overwrite the penetrance in the PanelPro Database
  PanelPRODatabase$Penetrance[cancer_index, gene_index, race_index, sex_index, , ] <- penetrance.mod.m

  # Storing the estimates
  log_likelihood <- 0

  for (i in 1:length(families)) {
    data <- families[[i]] # Get the data for the current family

    # Access the posterior probabilities (not normalized) and estimates for the specified gene and cancer type
    postprobs <- PPP(pedigree = data, genes = c(gene_input), cancers = cancer_type, database = PanelPRODatabase, impute.missing.ages = FALSE)$posterior.prob[[1]]
    estimate <- postprobs[postprobs$genes == gene_adj, "estimate"]

    # Check for NA or NaN before proceeding to if condition
       if (is.nan(estimate) || estimate <= 0) {
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

