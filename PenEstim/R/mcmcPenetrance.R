#' Execution of a Single Chain in Metropolis-Hastings for Cancer Risk Estimation
#'
#' Performs a single chain execution in the Metropolis-Hastings algorithm for Bayesian inference,
#' specifically tailored for cancer risk estimation. It estimates parameters related to cancer penetrance
#' based on family data, genetic information, and SEER database estimates.
#'
#' @param seed Seed value for random number generation.
#' @param n_iter Number of iterations for the chain.
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
#' @param homozygote Boolean indicating whether to exclude the possibility of homozygous carriers of the pathogenic variant in the genetic model for penetrance estimation.
#' @param SeerNC Boolean indicating if non-carrier penetrance is assumed to be the SEER penetrance.
#' @param sex Specifies the sex for which the estimation is performed; can be "NA" (default), "Female", or "Male".
#' @return A list containing samples, log likelihoods, acceptance ratio, and rejection rate for each iteration.
#' @importFrom stats set.seed
#' @importFrom parallel makeCluster stopCluster parLapply
#' @examples
#' result <- mhChain(
#'   seed = 123, n_iter = 1000, chain_id = 1, data = familyData,
#'   max_age = 90, db = database,
#'   prior_distributions = propDist, cancer_type = "breast",
#'   gene_input = "BRCA1", af = 0.0001, median_max = TRUE,
#'   max_penetrance = 1, homozygote = TRUE, SeerNC = TRUE, sex = "NA"
#' )
#' @export
# Main mhChain function
mhChain <- function(
    seed, n_iter, chain_id, data,
    max_age, db,
    prior_distributions, cancer_type, gene_input, af,
    median_max, max_penetrance, homozygote, SeerNC, sex) {
  # Set seed
  set.seed(seed)

  # Calculate SEER baseline and midpoint
  SEER_baseline <- calculate_lifetime_risk(cancer = cancer_type, gene = "SEER", race = "All_Races", sex=sex, type = "Net", db = db)
  midpoint_prob <- SEER_baseline$total_prob / 2
  midpoint_index <- which(SEER_baseline$cumulative_risk >= midpoint_prob)[1]
  baseline_mid <- as.numeric(names(SEER_baseline$cumulative_risk)[midpoint_index])

  draw_initial_params <- function() {
    asymptote_factor <- 2 * max_penetrance - SEER_baseline$total_prob
    if (asymptote_factor > 1) {
      asymptote_factor <- 1 - SEER_baseline$total_prob
    }
    asymptote <- SEER_baseline$total_prob +
      do.call(prior_distributions$asymptote_distribution, list(1)) * asymptote_factor
    asymptote <- max(0, min(1, asymptote))
    threshold <- do.call(prior_distributions$threshold_distribution, list(1))
    median <- if (median_max) {
      do.call(prior_distributions$median_distribution, list(1)) * (baseline_mid - threshold) + threshold
    } else {
      do.call(prior_distributions$median_distribution, list(1)) * (max_age - threshold) + threshold
    }
    first_quartile <- do.call(prior_distributions$first_quartile_distribution, list(1)) * (median - threshold) + threshold
    return(c(first_quartile, median, threshold, asymptote))
  }

  # Draw initial valid parameters
  initial_params <- draw_initial_params()
  first_quartile_current <- initial_params[1]
  median_current <- initial_params[2]
  threshold_current <- initial_params[3]
  asymptote_current <- initial_params[4]

  num_rejections <- 0

  # Set up an object to record the results for one chain
  out <- list(
    median_samples = numeric(n_iter),
    threshold_samples = numeric(n_iter),
    first_quartile_samples = numeric(n_iter),
    asymptote_samples = numeric(n_iter),
    loglikelihood_current = numeric(n_iter),
    loglikelihood_proposal = numeric(n_iter),
    acceptance_ratio = numeric(n_iter),
    rejection_rate = numeric(n_iter)
  )

  cat("Starting Chain", chain_id, "\n")

  for (i in 1:n_iter) {
    # Propose new values using the defined prior distributions for the asymptote, threshold,
    # median and then first quartile
    # Asymptote proposal is scaled to lie between either 2 * max_penetrane (a user-defined asymptote value
    # based on prior information) or 1 as the upper bound and the total cumulative probability of
    # the SEER baseline as the lower bound
    asymptote_factor <- 2 * max_penetrance - SEER_baseline$total_prob
    if (asymptote_factor > 1) {
      asymptote_factor <- 1 - SEER_baseline$total_prob
    }
    asymptote_proposal <- SEER_baseline$total_prob +
      do.call(prior_distributions$asymptote_distribution, list(1)) * asymptote_factor
    asymptote_proposal <- max(0, min(1, asymptote_proposal))
    # threshold proposal
    threshold_proposal <- do.call(prior_distributions$threshold_distribution, list(1))
    # Median proposal scaled to lie between either the median SEER age (basedline_mid), per default,
    # or the max_age as an upper bound and threshold proposal as a lower bound
    median_proposal <- if (median_max) {
      do.call(prior_distributions$median_distribution, list(1)) * (baseline_mid - threshold_proposal) + threshold_proposal
    } else {
      do.call(prior_distributions$median_distribution, list(1)) * (max_age - threshold_proposal) + threshold_proposal
    }
    # First Quartile proposal scaled to lie between the median proposal as an upper bound and threshold
    # proposal as a lower bound
    first_quartile_proposal <- do.call(prior_distributions$first_quartile_distribution, list(1)) *
      (median_proposal - threshold_proposal) + threshold_proposal

      # Compute the likelihood for the current and proposed
    loglikelihood_current <- mhLogLikelihood_clipp(
      paras = c(
        median_current,
        first_quartile_current, threshold_current,
        asymptote_current
      ), families = data,
      max_age = max_age,
      cancer_type = cancer_type,
      db = db,
      af = af,
      homozygote = homozygote,
      SeerNC = SeerNC,
      sex = sex
    )

    loglikelihood_proposal <- mhLogLikelihood_clipp(
      paras =
        c(
          median_proposal, first_quartile_proposal,
          threshold_proposal, asymptote_proposal
        ), families = data,
      max_age = max_age,
      cancer_type = cancer_type,
      db = db,
      af = af,
      homozygote = homozygote,
      SeerNC = SeerNC,
      sex = sex
    )

    # Compute the acceptance ratio (likelihood ratio)
    acceptance_ratio <- exp(loglikelihood_proposal - loglikelihood_current)

    # Accept or reject the proposal
    if (runif(1) < acceptance_ratio) {
      median_current <- median_proposal
      threshold_current <- threshold_proposal
      first_quartile_current <- first_quartile_proposal
      asymptote_current <- asymptote_proposal
    } else {
      num_rejections <- num_rejections + 1 # Increment the rejection counter
    }

    # Update the outputs
    out$median_samples[i] <- median_current
    out$threshold_samples[i] <- threshold_current
    out$first_quartile_samples[i] <- first_quartile_current
    out$asymptote_samples[i] <- asymptote_current
    out$loglikelihood_current[i] <- loglikelihood_current
    out$loglikelihood_proposal[i] <- loglikelihood_proposal
    out$acceptance_ratio[i] <- acceptance_ratio
    out$rejection_rate <- num_rejections / n_iter
  }

  # Return the result as a list
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
#' Unknown or missing values should be explicitly coded as `NA`.
#' * `ID`: A numeric value; ID for each individual. There should not be any
#' duplicated entries.
#' * `Sex`: A numeric value; `0` for female and `1` for male. Missing entries
#' are not currently supported.
#' * `MotherID`: A numeric value; unique ID for someone's mother.
#' * `FatherID`: A numeric value; unique ID for someone's father.
#' * `isProband`: A numeric value; `1` if someone is a proband, `0` otherwise.
#' This will be overridden by the `proband` argument in `PanelPRO`, if it is
#' specified. At least one proband should be specified by either the
#' `isProband` column or `proband`. Multiple probands are supported.
#' * `CurAge`: A numeric value; the age of censoring (current age if the person
#' is alive or age of death if the person is dead). Ages ranging from `1` to
#' `94` are allowed.
#' * `isAffX`: A numeric value; the affection status of cancer `X`, where `X`
#' is a `short` cancer code (see Details). Affection status should be encoded
#' as `1` if the individual was diagnosed, `0 `otherwise. Missing entries are
#' not currently supported.
#' * `AgeX`: A numeric value; the age of diagnosis for cancer `X`, where `X` is
#' a `short` cancer code (see Details). Ages ranging from `1` to `94` are
#' allowed. If the individual was not diagnosed for a given cancer, their
#' affection age should be encoded as `NA` and will be ignored otherwise.
#' * `isDead`: A numeric value; `1` if someone is dead, `0` otherwise. Missing
#' entries are assumed to be `0`.
#' * Columns for germline testing results (e.g. `BRCA1`,
#' `MLH1`) or tumor marker testing results. `ER`, `PR`, `CK14`, `CK5.6` and
#' `HER2` are tumor markers associated with breast cancer that will modify the
#' likelihood of phenotypes associated with `BRCA1` and `BRCA2`. `MSI` is a
#' tumor marker for colorectal cancer that will modify the likelihoods
#' associated with `MLH1`, `MSH2` `MSH6` and `PMS2`. For each of these optional
#' columns, positive results should be coded as `1`, negative results should be
#' coded as `0`, and unknown results should be coded as `NA`.
#' * Optional: `race`: A character string; expected values are `"All_Races"`, `"AIAN"`
#' (American Indian and Alaska Native), `"Asian"`, `"Black"`, `"White"`,
#' `"Hispanic"`, `"WH"` (white Hispanic), and `"WNH"` (non-white Hispanic)
#' (see `PanelPRO:::RACE_TYPES`). Asian-Pacific Islanders should be encoded as
#' `"Asian"`. Race information will be used to select the cancer and death by
#' other causes penetrances used in the model. Missing entries are recoded as
#' the `unknown.race` argument, which defaults to
#' `PanelPRO:::UNKNOWN_RACE`.
#' * Optional: `Ancestry`: A character string; expected values are `"AJ"`(Ashkenazi Jewish),
#' `"nonAJ"`, and `"Italian"` (see `PanelPRO:::ANCESTRY_TYPES`). The ancestry
#' information is used to determine the allele frequencies of BRCA1 and BRCA2,
#' if they are included in the model. Missing entries are recoded as the
#' `unknown.ancestry` argument, which defaults to `PanelPRO:::UNKNOWN_ANCESTRY`.
#' * `Twins`: A numeric value; `0` for non-identical/single births, `1` for
#' the first set of identical twins/multiple births in the family, `2` for the
#' second set, etc. Missing entries are assumed to be `0`.
#' @param cancer_type The type of cancer for which to estimate penetrance.
#' @param gene_input Gene information used for risk estimation.
#' @param n_chains Number of chains for parallel computation.
#' @param n_iter_per_chain Number of iterations for each chain.
#' @param db Database for the the baseline risk estimates. The default uses the the
#' PanelPRODatabase.
#' @param max_age Maximum age considered for analysis, default is 94.
#' @param removeProband Logical, indicating whether to remove probands from the
#' analysis (default is FALSE).
#' @param median_max Boolean indicating whether to use SEER median age or max_age as an
#' upper bound for the median proposal. Defaults to TRUE, i.e. using the SEER median.
#' @param sex Option to run the estimation for each sex seperately. The default is estimation 
#' of one penetrance curve for both sexes (sex = "NA"). Other settings are sex = "Female" and 
#' sex = "Male". 
#' @param homozygote Boolean to exclude the possibility of homozygous carriers of the PV in the 
#' genetic model for the penetrance estimation. Default setting is "True", which sets the likelihood
#' of a homozygous carrier to zero.
#' @param SeerNC Boolean with default value "True" indicating that the non-carrier penetrance
#'  is assumed to be the SEER penetrance. If "False", the non-carrier penetrance is calculated 
#' using the carrier penetrance in calculateNCPen. 
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
#' @param probCI Probability level for confidence intervals in penetrance plots
#' (default is 0.95).
#' @return A list containing combined results from all chains, along with optional
#' statistics and plots.
#' @importFrom stats rbeta runif
#' @importFrom parallel makeCluster stopCluster parLapply
#' @examples
#' # Example usage:
#' result <- PenEstim(
#'   data = familyData, cancer_type = "Breast", gene_input = "BRCA1",
#'   n_chains = 4, n_iter_per_chain = 1000, max_age = 90,
#'   burn_in = 0.1, thinning_factor = 2, summary_stats = TRUE,
#'   rejection_rates = TRUE, density_plots = TRUE, penetrance_plot = TRUE
#' )
#' @export

PenEstim <- function(data, cancer_type, gene_input, n_chains = 4,
                     n_iter_per_chain = 200,
                     db = PPP::PanelPRODatabase,
                     sex ="NA",
                     max_age = 94,
                     removeProband = FALSE,
                     median_max = TRUE,
                     homozygote = TRUE,
                     SeerNC = TRUE,
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

  # create the seeds for the individual chains
  seeds <- sample.int(1000, n_chains)

  # Apply the prepAges function to treat missing data
  data <- prepAges(data, removeProband)

  # Apply the transformation to adjust the format for the clipp package
  data <- do.call(rbind, lapply(data, transformDF,
    cancer_type = cancer_type,
    gene = gene_input
  ))

  #  Create the prior distributions
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
    library(dplyr)
  })

  parallel::clusterExport(cl, c(
    "mhChain", "mhLogLikelihood_clipp", "calculate_lifetime_risk", "calculateNCPen",
    "calculate_weibull_parameters", "validate_weibull_parameters", "calculateBaseline",
    "transformDF", "makePriors", "penet.fn",
    "seeds", "n_iter_per_chain", "sex",
    "data", "prop", "af", "max_age", "homozygote", "SeerNC", "median_max",
    "PanelPRODatabase", "cancer_type", "gene_input", "CANCER_TYPES",
    "GENE_TYPES", "CANCER_NAME_MAP"
  ), envir = environment())

  results <- parallel::parLapply(cl, 1:n_chains, function(i) {
    mhChain(seeds[i],
      n_iter = n_iter_per_chain, chain_id = i,
      data = data,
      db = db,
      prior_distributions = prop,
      max_age = max_age,
      cancer_type = cancer_type,
      gene_input = gene_input,
      af = af,
      max_penetrance = max_penetrance,
      median_max = median_max,
      homozygote = homozygote,
      SeerNC = SeerNC,
      sex = sex
    )
  })

  parallel::stopCluster(cl)

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