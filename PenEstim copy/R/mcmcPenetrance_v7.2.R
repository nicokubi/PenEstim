#' Execution of a Single Chain in Metropolis-Hastings for Cancer Risk Estimation
#'
#' Performs a single chain execution in the Metropolis-Hastings algorithm for Bayesian inference,
#' specifically tailored for cancer risk estimation. It estimates parameters related to cancer penetrance
#' based on family data, genetic information, and SEER database estimates.
#'
#' @param seed Seed value for random number generation.
#' @param n_iter Number of iterations for the chain.
#' @param burn_in Fraction of iterations to discard as burn-in (0 to 1).
#' @param chain_id Identifier for the chain.
#' @param data List of families data.
#' @param max_age Maximum age to be considered.
#' @param db Database containing baseline risk estimates.
#' @param prior_distributions List of the parameters for the distributions of the proposal, including asymptote, threshold, median, and first quartile distributions.
#' @param cancer_type Type of cancer for which risk is being estimated.
#' @param gene_input Gene information for risk estimation.
#' @param af Allele frequency for the risk allele.
#' @param median_max Boolean indicating whether to use SEER median or max_age as an upper bound for the median proposal. Defaults to TRUE, i.e., using the SEER median.
#' @param max_penetrance Maximum penetrance considered for analysis.
#' @param SeerNC Boolean indicating if non-carrier penetrance is assumed to be the SEER penetrance.
#' @param sex Specifies the sex for which the estimation is performed; can be "NA" (default), "Female", or "Male".
#' @param var Variance-covariance matrix used in the proposal distribution.
#' @param ageImputation Boolean indicating whether to impute ages for missing data.
#' @return A list containing samples, log likelihoods, acceptance ratio, and rejection rate for each iteration.
#' @importFrom stats set.seed
#' @importFrom parallel makeCluster stopCluster parLapply
#' @examples
#' result <- mhChain_v7.2(
#'   seed = 123, n_iter = 1000, burn_in = 0.1, chain_id = 1, data = familyData,
#'   max_age = 90, db = database,
#'   prior_distributions = propDist, cancer_type = "breast",
#'   gene_input = "BRCA1", af = 0.0001, median_max = TRUE,
#'   max_penetrance = 1, SeerNC = TRUE, sex = "NA",
#'   var = c(0.1, 0.1, 2, 2, 5, 5, 5, 5), ageImputation = FALSE
#' )
#' @export
mhChain_v7.2 <- function(seed, n_iter, burn_in, chain_id, ncores, data, max_age, db,
                         prior_distributions, cancer_type, gene_input, af,
                         median_max, max_penetrance, SeerNC, var,
                         ageImputation) {
  
  # Set seed
  set.seed(seed)
  # Calculate Empirical density
  age_density <- calculateEmpiricalDensity(data, aff_column = "aff", age_column = "age")
  
  if (ageImputation) {
    data <- calcPedDegree(data)
    # Initialize ages
    threshold <- prior_distributions$prior_params$threshold$min
    init_result <- imputeAgesInit(data, threshold, max_age)
    data <- init_result$data
    na_indices <- init_result$na_indices
  }

  # Calculate SEER baseline and midpoint
  SEER_baseline <- calculate_lifetime_risk(
    cancer = cancer_type, gene = "SEER",
    race = "All_Races", type = "Net", db = db
  )
  # Normalize CDF for males and females
  SEER_male <- data.frame(
    age = as.numeric(names(SEER_baseline$cumulative_risk$male)),
    cum_prob = SEER_baseline$cumulative_risk$male / max(SEER_baseline$cumulative_risk$male)
  )
  SEER_female <- data.frame(
    age = as.numeric(names(SEER_baseline$cumulative_risk$female)),
    cum_prob = SEER_baseline$cumulative_risk$female / max(SEER_baseline$cumulative_risk$female)
  )
  midpoint_prob_male <- SEER_baseline$lifetime_risk$male / 2
  midpoint_prob_female <- SEER_baseline$lifetime_risk$female / 2
  midpoint_index_male <-
    which(SEER_baseline$cumulative_risk$male >= midpoint_prob_male)[1]
  midpoint_index_female <-
    which(SEER_baseline$cumulative_risk$female >= midpoint_prob_female)[1]
  baseline_mid_male <-
    as.numeric(names(SEER_baseline$cumulative_risk$male)[midpoint_index_male])
  baseline_mid_female <-
    as.numeric(names(SEER_baseline$cumulative_risk$female)[midpoint_index_female])
  
  # Function to initialize the Weibull parameters using empirical data
  draw_initial_params <- function(data, prior_distributions) {
    # Filter data by sex and affected status
    data_male_affected <- data[data$sex == 1 & data$aff == 1, ]
    data_female_affected <- data[data$sex == 2 & data$aff == 1, ]
    
    # Calculate threshold (minimum age), median, and first quartile by sex among affected individuals
    lower_bound <- prior_distributions$prior_params$threshold$min
    upper_bound <- prior_distributions$prior_params$threshold$max
    
    threshold_male <- ifelse(length(data_male_affected$age) > 0,
                             min(data_male_affected$age, na.rm = TRUE), NA
    )
    threshold_female <- ifelse(length(data_female_affected$age) > 0,
                               min(data_female_affected$age, na.rm = TRUE), NA
    )
    
    threshold_male <- pmax(pmin(threshold_male, upper_bound, na.rm = TRUE), lower_bound, na.rm = TRUE)
    threshold_female <- pmax(pmin(threshold_female, upper_bound, na.rm = TRUE), lower_bound, na.rm = TRUE)
    
    median_male <- ifelse(length(data_male_affected$age) > 0,
                          median(data_male_affected$age, na.rm = TRUE), NA
    )
    median_female <- ifelse(length(data_female_affected$age) > 0,
                            median(data_female_affected$age, na.rm = TRUE), NA
    )
    
    first_quartile_male <- ifelse(length(data_male_affected$age) > 0,
                                  quantile(data_male_affected$age, probs = 0.25, na.rm = TRUE), NA
    )
    first_quartile_female <- ifelse(length(data_female_affected$age) > 0,
                                    quantile(data_female_affected$age, probs = 0.25, na.rm = TRUE), NA
    )
    
    asymptote_male <- runif(1,SEER_baseline$cumulative_risk$male[length(SEER_baseline$cumulative_risk$male)],1)
    asymptote_female <- runif(1,SEER_baseline$cumulative_risk$female[length(SEER_baseline$cumulative_risk$female)],1)
    
    return(list(
      asymptote_male = asymptote_male,
      asymptote_female = asymptote_female,
      threshold_male = threshold_male,
      threshold_female = threshold_female,
      median_male = median_male,
      median_female = median_female,
      first_quartile_male = first_quartile_male,
      first_quartile_female = first_quartile_female
    ))
  }
  
  # Initialize Parameters
  initial_params <- draw_initial_params(data = data, prior_distributions = prior_distributions)
  params_current <- initial_params
  current_states <- list()
  
  num_pars <- 8
  C <- diag(var)
  sd <- 2.38^2 / num_pars
  eps <- 0.01
  
  out <- list(
    asymptote_male_samples = numeric(n_iter),
    asymptote_female_samples = numeric(n_iter),
    asymptote_male_proposals = numeric(n_iter),
    asymptote_female_proposals = numeric(n_iter),
    threshold_male_samples = numeric(n_iter),
    threshold_female_samples = numeric(n_iter),
    threshold_male_proposals = numeric(n_iter),
    threshold_female_proposals = numeric(n_iter),
    median_male_samples = numeric(n_iter),
    median_female_samples = numeric(n_iter),
    median_male_proposals = numeric(n_iter),
    median_female_proposals = numeric(n_iter),
    first_quartile_male_samples = numeric(n_iter),
    first_quartile_female_samples = numeric(n_iter),
    first_quartile_male_proposals = numeric(n_iter),
    first_quartile_female_proposals = numeric(n_iter),
    loglikelihood_current = numeric(n_iter),
    loglikelihood_proposal = numeric(n_iter),
    logprior_current = numeric(n_iter),
    logprior_proposal = numeric(n_iter),
    acceptance_ratio = numeric(n_iter),
    rejection_rate = numeric(n_iter),
    C = vector("list", n_iter)
  )
  
  num_rejections <- 0
  cat("Starting Chain", chain_id, "\n")
  
  calculate_log_prior <- function(params, prior_distributions, max_age) {
    prior_params <- prior_distributions$prior_params

    scaled_asymptote_male <- params$asymptote_male
    scaled_asymptote_female <- params$asymptote_female

    scaled_threshold_male <- params$threshold_male
    scaled_threshold_female <- params$threshold_female

    scaled_median_male <- (params$median_male - params$threshold_male) / (max_age - params$threshold_male)
    scaled_median_female <- (params$median_female - params$threshold_female) / (max_age - params$threshold_female)

    scaled_first_quartile_male <- (params$first_quartile_male - params$threshold_male) /
      (params$median_male - params$threshold_male)
    scaled_first_quartile_female <- (params$first_quartile_female - params$threshold_female) /
      (params$median_female - params$threshold_female)

    log_prior_asymptote_male <- dbeta(scaled_asymptote_male, prior_params$asymptote$g1, prior_params$asymptote$g2, log = TRUE)
    log_prior_asymptote_female <- dbeta(scaled_asymptote_female, prior_params$asymptote$g1, prior_params$asymptote$g2, log = TRUE)

    log_prior_threshold_male <- dunif(scaled_threshold_male, prior_params$threshold$min, prior_params$threshold$max, log = TRUE)
    log_prior_threshold_female <- dunif(scaled_threshold_female, prior_params$threshold$min, prior_params$threshold$max, log = TRUE)

    log_prior_median_male <- dbeta(scaled_median_male, prior_params$median$m1, prior_params$median$m2, log = TRUE)
    log_prior_median_female <- dbeta(scaled_median_female, prior_params$median$m1, prior_params$median$m2, log = TRUE)

    log_prior_first_quartile_male <- dbeta(scaled_first_quartile_male, prior_params$first_quartile$q1, prior_params$first_quartile$q2, log = TRUE)
    log_prior_first_quartile_female <- dbeta(scaled_first_quartile_female, prior_params$first_quartile$q1, prior_params$first_quartile$q2, log = TRUE)

    log_prior_total <- log_prior_asymptote_male + log_prior_asymptote_female +
      log_prior_threshold_male + log_prior_threshold_female +
      log_prior_median_male + log_prior_median_female +
      log_prior_first_quartile_male + log_prior_first_quartile_female

    return(log_prior_total)
  }

  
  for (i in 1:n_iter) {
    weibull_params_male <- calculate_weibull_parameters(params_current$median_male, params_current$first_quartile_male, params_current$threshold_male)
    alpha_male <-  weibull_params_male$alpha
    beta_male <- weibull_params_male$beta
    delta_male <- params_current$threshold_male
    
    weibull_params_female <- calculate_weibull_parameters(params_current$median_female, params_current$first_quartile_female, params_current$threshold_female)
    alpha_female <- weibull_params_female$alpha
    beta_female <- weibull_params_female$beta
    delta_female <- params_current$threshold_female
    
    # Impute ages 
    if (ageImputation) {
    data <- imputeAges(data, na_indices, SEER_male, SEER_female, alpha_male, beta_male, delta_male,
                       alpha_female, beta_female, delta_female)
    }
    
    params_vector <- c(
      params_current$asymptote_male, params_current$asymptote_female,
      params_current$threshold_male, params_current$threshold_female,
      params_current$median_male, params_current$median_female,
      params_current$first_quartile_male, params_current$first_quartile_female
    )
    
    proposal_vector <- mvrnorm(1, mu = params_vector, Sigma = C)
    
    proposal_vector[1] <- ifelse(proposal_vector[1] < 0, -proposal_vector[1], 
                                 ifelse(proposal_vector[1] > 1, 2 - proposal_vector[1], proposal_vector[1]))
    proposal_vector[2] <- ifelse(proposal_vector[2] < 0, -proposal_vector[2], 
                                 ifelse(proposal_vector[2] > 1, 2 - proposal_vector[2], proposal_vector[2]))
    
    out$asymptote_male_proposals[i] <- proposal_vector[1]
    out$asymptote_female_proposals[i] <- proposal_vector[2]
    out$threshold_male_proposals[i] <- proposal_vector[3]
    out$threshold_female_proposals[i] <- proposal_vector[4]
    out$median_male_proposals[i] <- proposal_vector[5]
    out$median_female_proposals[i] <- proposal_vector[6]
    out$first_quartile_male_proposals[i] <- proposal_vector[7]
    out$first_quartile_female_proposals[i] <- proposal_vector[8]
    
    params_proposal <- list(
      asymptote_male = proposal_vector[1],
      asymptote_female = proposal_vector[2],
      threshold_male = proposal_vector[3],
      threshold_female = proposal_vector[4],
      median_male = proposal_vector[5],
      median_female = proposal_vector[6],
      first_quartile_male = proposal_vector[7],
      first_quartile_female = proposal_vector[8]
    )
    
    loglikelihood_current <- mhLogLikelihood_clipp(
      params_current, data, max_age,
      cancer_type, db, af, SeerNC, ncores
    )
    
    logprior_current <- calculate_log_prior(params_current, prior_distributions, max_age)
    
    out$loglikelihood_current[i] <- loglikelihood_current
    out$logprior_current[i] <- logprior_current
    
    out$loglikelihood_proposal[i] <- NA
    out$logprior_proposal[i] <- NA
    out$acceptance_ratio[i] <- NA
    
    is_rejected <- FALSE
    if (
      proposal_vector[1] < 0 || proposal_vector[1] > 1 ||
      proposal_vector[2] < 0 || proposal_vector[2] > 1 ||
      proposal_vector[3] < 0 || proposal_vector[3] > 100 ||
      proposal_vector[4] < 0 || proposal_vector[4] > 100 ||
      proposal_vector[5] < proposal_vector[7] || 
      (median_max && proposal_vector[5] > baseline_mid_male) || 
      (!median_max && proposal_vector[5] > max_age) || 
      proposal_vector[6] < proposal_vector[8] || 
      (median_max && proposal_vector[6] > baseline_mid_female) || 
      (!median_max && proposal_vector[6] > max_age) || 
      proposal_vector[7] < proposal_vector[3] || proposal_vector[7] > proposal_vector[5] ||
      proposal_vector[8] < proposal_vector[4] || proposal_vector[8] > proposal_vector[6]) {
      is_rejected <- TRUE
      num_rejections <- num_rejections + 1
    } else {
      loglikelihood_proposal <- mhLogLikelihood_clipp(
        params_proposal, data, max_age,
        cancer_type, db, af, SeerNC, ncores
      )
      logprior_proposal <- calculate_log_prior(params_proposal, prior_distributions, max_age)
      
      log_acceptance_ratio <- (loglikelihood_proposal + logprior_proposal) - (loglikelihood_current + logprior_current)
      
      if (log(runif(1)) < log_acceptance_ratio) {
        params_current <- params_proposal
      } else {
        num_rejections <- num_rejections + 1
      }
      
      out$loglikelihood_proposal[i] <- loglikelihood_proposal
      out$logprior_proposal[i] <- logprior_proposal
      out$acceptance_ratio[i] <- log_acceptance_ratio
    }
    
    current_states[[i]] <- c(
      params_current$asymptote_male, params_current$asymptote_female,
      params_current$threshold_male, params_current$threshold_female,
      params_current$median_male, params_current$median_female,
      params_current$first_quartile_male, params_current$first_quartile_female
    )
    
    if (i > max(burn_in * n_iter, 3)) {
      C <- sd * cov(do.call(rbind, current_states)) + eps * sd * diag(num_pars)
    }
    
    out$asymptote_male_samples[i] <- params_current$asymptote_male
    out$asymptote_female_samples[i] <- params_current$asymptote_female
    out$threshold_male_samples[i] <- params_current$threshold_male
    out$threshold_female_samples[i] <- params_current$threshold_female
    out$median_male_samples[i] <- params_current$median_male
    out$median_female_samples[i] <- params_current$median_female
    out$first_quartile_male_samples[i] <- params_current$first_quartile_male
    out$first_quartile_female_samples[i] <- params_current$first_quartile_female
    out$C[[i]] <- C
  }
  
  out$rejection_rate <- num_rejections / n_iter
  
  return(out)
}

#' Bayesian Inference using Independent Metropolis-Hastings for Penetrance Estimation
#'
#' This function implements the Independent Metropolis-Hastings algorithm for Bayesian
#' penetrance estimation of cancer risk. It utilizes parallel computing to run multiple
#' chains and provides various options for analyzing and visualizing the results.
#'
#' @param data The pedigree data in the required format. A data frame of
#' family history information with the following columns.
#' * `ID`: A numeric value; ID for each individual. There should not be any duplicated entries.
#' * `Sex`: A numeric value; `0` for female and `1` for male. Missing entries are not currently supported.
#' * `MotherID`: A numeric value; unique ID for someone's mother.
#' * `FatherID`: A numeric value; unique ID for someone's father.
#' * `isProband`: A numeric value; `1` if someone is a proband, `0` otherwise.
#' * `CurAge`: A numeric value; the age of censoring (current age if the person is alive or age of death if the person is dead). Ages ranging from `1` to `94` are allowed.
#' * `isAffX`: A numeric value; the affection status of cancer `X`, where `X` is a `short` cancer code (see Details). Affection status should be encoded as `1` if the individual was diagnosed, `0` otherwise. Missing entries are not currently supported.
#' * `AgeX`: A numeric value; the age of diagnosis for cancer `X`, where `X` is a `short` cancer code (see Details). Ages ranging from `1` to `94` are allowed. If the individual was not diagnosed for a given cancer, their affection age should be encoded as `NA`.
#' * `isDead`: A numeric value; `1` if someone is dead, `0` otherwise. Missing entries are assumed to be `0`.
#' * Columns for germline testing results (e.g., `BRCA1`, `MLH1`) or tumor marker testing results. Positive results should be coded as `1`, negative results should be coded as `0`, and unknown results should be coded as `NA`.
#' * Optional: `race`: A character string; expected values are `"All_Races"`, `"AIAN"`, `"Asian"`, `"Black"`, `"White"`, `"Hispanic"`, `"WH"`, and `"WNH"`.
#' * Optional: `Ancestry`: A character string; expected values are `"AJ"`, `"nonAJ"`, and `"Italian"`.
#' * `Twins`: A numeric value; `0` for non-identical/single births, `1` for the first set of identical twins/multiple births in the family, `2` for the second set, etc.
#' @param cancer_type The type of cancer for which to estimate penetrance.
#' @param gene_input Gene information used for risk estimation.
#' @param n_chains Number of chains for parallel computation.
#' @param n_iter_per_chain Number of iterations for each chain.
#' @param db Database for the baseline risk estimates.
#' @param max_age Maximum age considered for analysis, default is 94.
#' @param removeProband Logical, indicating whether to remove probands from the analysis (default is FALSE).
#' @param median_max Boolean indicating whether to use SEER median age or max_age as an upper bound for the median proposal. Defaults to TRUE.
#' @param sex Option to run the estimation for each sex separately. The default is estimation of one penetrance curve for both sexes (sex = "NA").
#' @param SeerNC Boolean indicating that the non-carrier penetrance is assumed to be the SEER penetrance. Default is TRUE.
#' @param burn_in Fraction of results to discard as burn-in (0 to 1). Default is 0 (no burn-in).
#' @param thinning_factor Factor by which to thin the results, default is 1 (no thinning).
#' @param distribution_data Data for generating prior distributions.
#' @param af Allele frequency for risk allele, default is 0.0001.
#' @param max_penetrance Maximum penetrance considered for analysis, default is 1.
#' @param sample_size Optional sample size for distribution generation.
#' @param ratio Optional ratio parameter for distribution generation.
#' @param prior_params Parameters for prior distributions.
#' @param risk_proportion Proportion of risk for distribution generation.
#' @param summary_stats Boolean to include summary statistics in the output.
#' @param rejection_rates Boolean to include rejection rates in the output.
#' @param density_plots Boolean to include density plots in the output.
#' @param penetrance_plot Boolean to include penetrance plots in the output.
#' @param probCI Probability level for confidence intervals in penetrance plots (default is 0.95).
#' @return A list containing combined results from all chains, along with optional statistics and plots.
#' @importFrom stats rbeta runif
#' @examples
#' result <- PenEstim_v7.2(
#'   data = familyData, cancer_type = "Breast", gene_input = "BRCA1",
#'   n_chains = 4, n_iter_per_chain = 1000, max_age = 90,
#'   burn_in = 0.1, thinning_factor = 2, summary_stats = TRUE,
#'   rejection_rates = TRUE, density_plots = TRUE, penetrance_plot = TRUE
#' )
#' @export
PenEstim_v7.2 <- function(data, cancer_type, gene_input, n_chains = 1,
                          n_iter_per_chain = 10000,
                          db = PPP::PanelPRODatabase,
                          ncores = 6,
                          sex = "NA",
                          max_age = 94,
                          removeProband = FALSE,
                          ageImputation = FALSE,
                          median_max = TRUE,
                          SeerNC = TRUE,
                          var = c(0.1, 0.1, 2, 2, 5, 5, 5, 5),
                          burn_in = 0,
                          thinning_factor = 1,
                          distribution_data = distribution_data_default,
                          af = PPP::PanelPRODatabase$AlleleFrequency[paste0(gene_input, "_anyPV"), "nonAJ"],
                          max_penetrance = 1,
                          sample_size = NULL,
                          ratio = NULL,
                          prior_params = prior_params_default,
                          risk_proportion = risk_proportion_default,
                          summary_stats = TRUE,
                          rejection_rates = TRUE,
                          density_plots = TRUE,
                          penetrance_plot = TRUE,
                          probCI = 0.95) {
  # Validate inputs
  if (missing(data)) {
    stop("Error: 'data' parameter is missing. Please provide a valid list of pedigrees.")
  }
  if (!(cancer_type %in% CANCER_TYPES)) {
    stop(paste("Error: Cancer type", shQuote(cancer_type), "is not supported. Please choose from the supported list."))
  }
  if (!(gene_input %in% GENE_TYPES)) {
    stop(paste("Error: Gene type", shQuote(gene_input), "is not supported. Please choose from the supported list."))
  }
  if (missing(n_chains) || !is.numeric(n_chains) || n_chains <= 0) {
    stop("Error: 'n_chains' parameter is missing or invalid. Please specify a positive integer.")
  }
  if (missing(n_iter_per_chain) || !is.numeric(n_iter_per_chain) || n_iter_per_chain <= 0) {
    stop("Error: 'n_iter_per_chain' parameter is missing or invalid. It must be a positive integer.")
  }
  if (n_chains > parallel::detectCores()) {
    stop("Error: 'n_chains' exceeds the number of available CPU cores.")
  }
  
  # Create the seeds for the individual chains
  seeds <- sample.int(1000, n_chains)
  
  # Apply the prepAges function to treat missing data
  data <- prepAges(data)
  
  # Apply the transformation to adjust the format for the clipp package
  data <- do.call(rbind, lapply(data, transformDF,
                                cancer_type = cancer_type,
                                gene = gene_input
  ))
  
  if (removeProband) {
    data <- data[data$isProband != 1, ]
  }
  
  # Create the prior distributions
  prop <- makePriors(
    data = distribution_data,
    sample_size = sample_size,
    cancer = cancer_type,
    ratio = ratio,
    prior_params = prior_params,
    risk_proportion = risk_proportion
  )
  
  cores <- parallel::detectCores()
  
  if (n_chains > cores) {
    stop("Error: 'n_chains exceeds the number of available CPU cores.")
  }
  cl <- parallel::makeCluster(n_chains)
  
  # Load required packages to the clusters
  parallel::clusterEvalQ(cl, {
    library(PPP)
    library(clipp)
    library(stats4)
    library(MASS)
    library(dplyr)
    library(parallel)
  })
  
  parallel::clusterExport(cl, c(
    "mhChain_v7.2", "mhLogLikelihood_clipp", "calculate_lifetime_risk", "calculateNCPen",
    "calculate_weibull_parameters", "validate_weibull_parameters", "calculateBaseline", "prior_params",
    "transformDF", "makePriors", "lik.fn", "mvrnorm", "var", "calculateEmpiricalDensity",
    "seeds", "n_iter_per_chain", "sex", "burn_in", "imputeAges", "imputeAgesInit", "draw_age_from_seer",
    "data", "prop", "af", "max_age", "SeerNC", "median_max", "ncores",
    "PanelPRODatabase", "cancer_type", "gene_input", "CANCER_TYPES",
    "GENE_TYPES", "CANCER_NAME_MAP"
  ), envir = environment())
  
  results <- parallel::parLapply(cl, 1:n_chains, function(i) {
    mhChain_v7.2(seeds[i],
                 n_iter = n_iter_per_chain,
                 burn_in = burn_in,
                 chain_id = i,
                 data = data,
                 db = db,
                 ncores = ncores,
                 prior_distributions = prop,
                 max_age = max_age,
                 cancer_type = cancer_type,
                 gene_input = gene_input,
                 af = af,
                 max_penetrance = max_penetrance,
                 median_max = median_max,
                 SeerNC = SeerNC,
                 var = var,
                 ageImputation = ageImputation
    )
  })
  
  # Check rejection rates and issue a warning if they are all above 90%
  all_high_rejections <- all(sapply(results, function(x) x$rejection_rate > 0.9))
  if (all_high_rejections) {
    warning("Low acceptance rate. Please consider running the chain longer.")
  }
  
  # Apply burn-in and thinning
  if (burn_in > 0) {
    results <- apply_burn_in(results, burn_in)
  }
  if (thinning_factor > 1) {
    results <- apply_thinning(results, thinning_factor)
  }
  
  # Extract samples from the chains
  combined_chains <- combine_chains(results)
  
  # Initialize variables
  output <- list()
  
  tryCatch(
    {
      if (rejection_rates) {
        # Generate rejection rates
        output$rejection_rates <- printRejectionRates(results)
      }
      
      if (summary_stats) {
        # Generate summary statistics
        output$summary_stats <- generate_summary(combined_chains)
      }
      
      if (density_plots) {
        # Generate density plots
        output$density_plots <- generate_density_plots(combined_chains)
      }
      
      if (penetrance_plot) {
        # Generate penetrance plot
        output$penetrance_plot <- plot_penetrance(combined_chains, prob = probCI, max_age = max_age)
      }
    },
    error = function(e) {
      # Handle errors here
      cat("An error occurred in the output display: ", e$message, "\n")
    }
  )
  
  output$combined_chains <- combined_chains
  output$results <- results
  output$data <- data
  
  return(output)
}