#' Execution of a Single Chain in Metropolis-Hastings for Cancer Risk Estimation
#'
#' Performs a single chain execution in the Metropolis-Hastings algorithm for Bayesian inference,
#' specifically tailored for cancer risk estimation. It estimates parameters related to cancer penetrance
#' based on family data, genetic information, and baseline database estimates.
#'
#' @param seed Integer, the seed for the random number generator to ensure reproducibility.
#' @param n_iter Integer, the number of iterations to perform in the Metropolis-Hastings algorithm.
#' @param burn_in Integer, the number of initial iterations to discard (burn-in period).
#' @param chain_id Integer, the identifier for the chain being executed.
#' @param ncores Integer, the number of cores to use for parallel computation.
#' @param data Data frame, containing family and genetic information used in the analysis.
#' @param max_age Integer, the maximum age considered in the analysis.
#' @param baseline_data Numeric matrix, containing baseline risk estimates for different ages and sexes.
#' @param prior_distributions List, containing prior distributions for the parameters being estimated.
#' @param af Numeric, the allele frequency of the risk allele in the population.
#' @param median_max Numeric, the maximum median age for the Weibull distribution.
#' @param max_penetrance Numeric, the maximum penetrance value allowed.
#' @param BaselineNC Logical, indicates if non-carrier penetrance should be based on SEER data.
#' @param var Numeric, the variance for the proposal distribution in the Metropolis-Hastings algorithm.
#' @param ageImputation Logical, indicates if age imputation should be performed.
#' @param removeProband Logical, indicates if the proband should be removed from the analysis.
#'
#' @return A list containing samples, log likelihoods, acceptance ratio, and rejection rate for each iteration.
#' @export
mhChain <- function(seed, n_iter, burn_in, chain_id, ncores, data, max_age, baseline_data,
                    prior_distributions, af, median_max, max_penetrance, BaselineNC, var,
                    ageImputation, removeProband) {
  # Set seed for the chain
  set.seed(seed)
  
  # Calculate empirical age density for affected individuals
  age_density <- calculateEmpiricalDensity(data, aff_column = "aff", age_column = "age")
  
  # Prepare initial age imputation if enabled
  if (ageImputation) {
    data <- calcPedDegree(data)
    # Initialize ages
    threshold <- prior_distributions$prior_params$threshold$min
    init_result <- imputeAgesInit(data, threshold, max_age)
    data <- init_result$data
    na_indices <- init_result$na_indices
  } else {
    # If age imputation is disabled, set unknown ages to 1 so they are disregarded in likelihood calculation
    data$age[is.na(data$age)] <- 1
  }
  
  # Option to remove the proband after age imputation
  if (removeProband) {
    data <- data[data$isProband != 1, ]
  }
  
  # Process baseline risk data
  baseline_male <- as.numeric(baseline_data[,"Male"])
  baseline_female <- as.numeric(baseline_data[,"Female"])
  
  baseline_male_cum <- cumsum(baseline_male)
  baseline_female_cum <- cumsum(baseline_female)
  
  baseline_male_df <- data.frame(
    age = 1:length(baseline_male),
    cum_prob = baseline_male_cum / max(baseline_male_cum)
  )
  baseline_female_df <- data.frame(
    age = 1:length(baseline_female),
    cum_prob = baseline_female_cum / max(baseline_female_cum)
  )
  
  midpoint_prob_male <- baseline_male_cum[length(baseline_male_cum)] / 2
  midpoint_prob_female <- baseline_female_cum[length(baseline_female_cum)] / 2
  
  midpoint_index_male <- which(baseline_male_cum >= midpoint_prob_male)[1]
  midpoint_index_female <- which(baseline_female_cum >= midpoint_prob_female)[1]
  
  baseline_mid_male <- midpoint_index_male
  baseline_mid_female <- midpoint_index_female
  
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
                                  min(quantile(data_male_affected$age, probs = 0.25, na.rm = TRUE), median_male - 1), NA
    )
    first_quartile_female <- ifelse(length(data_female_affected$age) > 0,
                                    min(quantile(data_female_affected$age, probs = 0.25, na.rm = TRUE), median_female - 1), NA
    )
    
    asymptote_male <- runif(1, max(baseline_male_cum), 1)
    asymptote_female <- runif(1, max(baseline_female_cum), 1)
    
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
  
  # Initialize Parameters using the function draw_initial_params
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
  
  # Function to calculate the (log) prior probabilities using the fixed prior distributions and specified prior parameters. 
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
  
  # Run n_iter iterations of the adaptive Metropolis-Hastings algorithm
  for (i in 1:n_iter) {
    # Calculate Weibull parameters from current parameters
    weibull_params_male <- calculate_weibull_parameters(params_current$median_male, params_current$first_quartile_male, params_current$threshold_male)
    alpha_male <- weibull_params_male$alpha
    beta_male <- weibull_params_male$beta
    delta_male <- params_current$threshold_male
    
    weibull_params_female <- calculate_weibull_parameters(params_current$median_female, params_current$first_quartile_female, params_current$threshold_female)
    alpha_female <- weibull_params_female$alpha
    beta_female <- weibull_params_female$beta
    delta_female <- params_current$threshold_female
    
    # Impute ages at each iteration based on current parameters
    if (ageImputation) {
      data <- imputeAges(
        data, na_indices, baseline_male_df, baseline_female_df, alpha_male, beta_male, delta_male,
        alpha_female, beta_female, delta_female
      )
    }
    
    # Store the current parameter values
    params_vector <- c(
      params_current$asymptote_male, params_current$asymptote_female,
      params_current$threshold_male, params_current$threshold_female,
      params_current$median_male, params_current$median_female,
      params_current$first_quartile_male, params_current$first_quartile_female
    )
    
    # Generate the proposal parameters from a multivariate normal distribution centered around current parameters
    proposal_vector <- mvrnorm(1, mu = params_vector, Sigma = C)
    
    # Ensure the proposals for the asymptote fall within the 0 to 1 range
    proposal_vector[1] <- ifelse(proposal_vector[1] < 0, -proposal_vector[1],
                                 ifelse(proposal_vector[1] > 1, 2 - proposal_vector[1], proposal_vector[1])
    )
    proposal_vector[2] <- ifelse(proposal_vector[2] < 0, -proposal_vector[2],
                                 ifelse(proposal_vector[2] > 1, 2 - proposal_vector[2], proposal_vector[2])
    )
    
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
    
    # Evaluate the current set of parameters
    loglikelihood_current <- mhLogLikelihood_clipp(
      params_current, data, max_age,
      baseline_data, af, BaselineNC, ncores
    )
    logprior_current <- calculate_log_prior(params_current, prior_distributions, max_age)
    # Record the outputs of the evaluation for the current set of parameters
    out$loglikelihood_current[i] <- loglikelihood_current
    out$logprior_current[i] <- logprior_current
    
    out$loglikelihood_proposal[i] <- NA
    out$logprior_proposal[i] <- NA
    out$acceptance_ratio[i] <- NA
    
    # Check that the proposed parameters satisfy the basic requirements
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
        baseline_data, af, BaselineNC, ncores
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
#' @param data A data frame containing the pedigree data in the required format. It should include the following columns:
#'   - `ID`: A numeric value representing the unique identifier for each individual. There should be no duplicated entries.
#'   - `Sex`: A numeric value where `0` indicates female and `1` indicates male. Missing entries are not currently supported.
#'   - `MotherID`: A numeric value representing the unique identifier for an individual's mother.
#'   - `FatherID`: A numeric value representing the unique identifier for an individual's father.
#'   - `isProband`: A numeric value where `1` indicates the individual is a proband and `0` otherwise.
#'   - `CurAge`: A numeric value indicating the age of censoring (current age if the person is alive or age at death if the person is deceased). Allowed ages range from `1` to `94`.
#'   - `isAff`: A numeric value indicating the affection status of cancer, with `1` for diagnosed individuals and `0` otherwise. Missing entries are not supported.
#'   - `Age`: A numeric value indicating the age of cancer diagnosis, encoded as `NA` if the individual was not diagnosed. Allowed ages range from `1` to `94`.
#'   - `isDead`: A numeric value where `1` indicates the individual is deceased and `0` indicates the individual is alive. Missing entries are assumed to be `0`.
#'   - `geno`: A column for germline testing or tumor marker testing results. Positive results should be coded as `1`, negative results as `0`, and unknown results as `NA` or left empty.
#'   - Optional: `race`: A character string with expected values `"All_Races"`, `"AIAN"`, `"Asian"`, `"Black"`, `"White"`, `"Hispanic"`, `"WH"`, and `"WNH"`.
#'   - Optional: `Ancestry`: A character string with expected values `"AJ"`, `"nonAJ"`, and `"Italian"`.
#'   - `Twins`: A numeric value where `0` indicates non-identical/single births, `1` indicates the first set of identical twins/multiple births in the family, `2` indicates the second set, and so on.
#' @param n_chains Number of chains for parallel computation. Default is 1.
#' @param n_iter_per_chain Number of iterations for each chain. Default is 10000.
#' @param ncores Number of cores for parallel computation. Default is 6.
#' @param baseline_data Data for the baseline risk estimates (probability of developing cancer), such as population-level risk from a cancer registry. Default is the allele frequency for MLH1 from the PanelPRO database.
#' @param max_age Maximum age considered for analysis. Default is 94.
#' @param removeProband Logical, indicating whether to remove probands from the analysis. Default is FALSE.
#' @param ageImputation Logical, indicating whether to perform age imputation. Default is FALSE.
#' @param median_max Boolean indicating whether to use the baseline median age or max_age as an upper bound for the median proposal. Default is TRUE.
#' @param BaselineNC Boolean indicating that the non-carrier penetrance is assumed to be the baseline penetrance. Default is TRUE.
#' @param var Vector of variances for the proposal distribution. Default is c(0.1, 0.1, 2, 2, 5, 5, 5, 5).
#' @param burn_in Fraction of results to discard as burn-in (0 to 1). Default is 0 (no burn-in).
#' @param thinning_factor Factor by which to thin the results. Default is 1 (no thinning).
#' @param distribution_data Data for generating prior distributions.
#' @param af Allele frequency for the risk allele. Default is 0.0001.
#' @param max_penetrance Maximum penetrance considered for analysis. Default is 1.
#' @param sample_size Optional sample size for distribution generation.
#' @param ratio Optional ratio parameter for distribution generation.
#' @param prior_params Parameters for prior distributions.
#' @param risk_proportion Proportion of risk for distribution generation.
#' @param summary_stats Boolean indicating whether to include summary statistics in the output. Default is TRUE.
#' @param rejection_rates Boolean indicating whether to include rejection rates in the output. Default is TRUE.
#' @param density_plots Boolean indicating whether to include density plots in the output. Default is TRUE.
#' @param penetrance_plot Boolean indicating whether to include penetrance plots in the output. Default is TRUE.
#' @param probCI Probability level for confidence intervals in penetrance plots. Default is 0.95.
#' @return A list containing combined results from all chains, along with optional statistics and plots.
#' @importFrom stats rbeta runif
#' @importFrom parallel makeCluster stopCluster parLapply
#' @export
PenEstim <- function(data, n_chains = 1,
                     n_iter_per_chain = 10000,
                     ncores = 6,
                     max_age = 94,
                     baseline_data = baseline_data_default,
                     removeProband = FALSE,
                     ageImputation = FALSE,
                     median_max = TRUE,
                     BaselineNC = TRUE,
                     var = c(0.1, 0.1, 2, 2, 5, 5, 5, 5),
                     burn_in = 0,
                     thinning_factor = 1,
                     distribution_data = distribution_data_default,
                     af = 0.0001,
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

  # Apply the transformation to adjust the format for the clipp package
  data <- do.call(rbind, lapply(data, transformDF))

  # Create the prior distributions
  prop <- makePriors(
    data = distribution_data,
    sample_size = sample_size,
    ratio = ratio,
    prior_params = prior_params,
    risk_proportion = risk_proportion, 
    baseline_data = baseline_data
  )

  cores <- parallel::detectCores()

  if (n_chains > cores) {
    stop("Error: 'n_chains exceeds the number of available CPU cores.")
  }
  cl <- parallel::makeCluster(n_chains)

  # Load required packages to the clusters
  parallel::clusterEvalQ(cl, {
    library(clipp)
    library(stats4)
    library(MASS)
    library(parallel)
    library(kinship2)
  })

  parallel::clusterExport(cl, c(
    "mhChain", "mhLogLikelihood_clipp", "calculate_weibull_parameters", "validate_weibull_parameters", "prior_params",
    "transformDF", "lik.fn", "mvrnorm", "var", "calculateEmpiricalDensity", "baseline_data", "calcPedDegree",
    "seeds", "n_iter_per_chain", "burn_in", "imputeAges", "imputeAgesInit", "drawBaseline", "calculateNCPen",
    "data", "prop", "af", "max_age", "BaselineNC", "median_max", "ncores", "removeProband"
  ), envir = environment())

  results <- parallel::parLapply(cl, 1:n_chains, function(i) {
    mhChain(seeds[i],
      n_iter = n_iter_per_chain,
      burn_in = burn_in,
      chain_id = i,
      data = data,
      ncores = ncores,
      prior_distributions = prop,
      max_age = max_age,
      af = af,
      max_penetrance = max_penetrance,
      median_max = median_max,
      baseline_data = baseline_data,
      BaselineNC = BaselineNC,
      var = var,
      ageImputation = ageImputation,
      removeProband = removeProband
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

