#' Execution of a Single Chain in Metropolis-Hastings
#'
#' This function performs a single chain execution in the Metropolis-Hastings algorithm
#' for Bayesian inference, specifically tailored for cancer risk estimation.
#'
#' @param seed Seed value for random number generation.
#' @param n_iter Number of iterations for the chain.
#' @param chain_id Identifier for the chain.
#' @param data List of families data.
#' @param max_age Maximum age to be considered.
#' @param PanelPRODatabase Database containing PanelPRO information.
#' @param proposal_distributions List of the parameters for the distributions of the proposal.
#' @param cancer_type Type of cancer for which risk is being estimated.
#' @param gene_input Gene information for risk estimation.
#' @param median_max Boolean indicating whether to use SEER median or max_age as an upper bound for the median proposal. Defaults to TRUE, i.e. using the SEER median.
#' @return A list containing samples, log likelihoods, acceptance ratio, and rejection rate for each iteration.
#' @importFrom stats set.seed
#' @importFrom parallel makeCluster stopCluster parLapply
#' @examples
#' # Example usage:
#' result <- mhChain(
#'   seed = 123, n_iter = 1000, chain_id = 1, data = familyData,
#'   max_age = 90, PanelPRODatabase = database,
#'   proposal_distributions = propDist, cancer_type = "breast",
#'   gene_input = "BRCA1", median_max = TRUE
#' )
#' @export

# Main mhChain function
mhChain <- function(
    seed, n_iter, chain_id, data,
    max_age, db,
    prior_distributions, cancer_type, gene_input,
    median_max = TRUE, max_penetrance) {
  # Set seed
  set.seed(seed)

  # Calculate SEER baseline and midpoint
  SEER_baseline <- calculate_lifetime_risk(cancer = cancer_type, gene = "SEER", data = db)
  midpoint_prob <- SEER_baseline$total_prob / 2
  midpoint_index <- which(SEER_baseline$cumulative_risk >= midpoint_prob)[1]
  baseline_mid <- as.numeric(names(SEER_baseline$cumulative_risk)[midpoint_index])

  draw_initial_params <- function() {
    repeat {
      asymptote_factor <- 2 * max_penetrance
      if (asymptote_factor > 1) {
        asymptote_factor <- 1 - SEER_baseline$total_prob
      }
      asymptote <- SEER_baseline$total_prob +
        do.call(prior_distributions$asymptote_distribution, list(1)) * asymptote_factor
      asymptote <- max(0, min(1, asymptote))
      shift <- do.call(prior_distributions$shift_distribution, list(1))
      median <- if (median_max) {
        do.call(prior_distributions$median_distribution, list(1)) * (baseline_mid - shift) + shift
      } else {
        do.call(prior_distributions$median_distribution, list(1)) * (max_age - shift) + shift
      }
      first_quartile <- do.call(prior_distributions$first_quartile_distribution, list(1)) * (median - shift) + shift

      if (validate_weibull_parameters(first_quartile, median, shift, asymptote)) {
        return(c(first_quartile, median, shift, asymptote))
      }
    }
  }

  # Draw initial valid parameters
  initial_params <- draw_initial_params()
  first_quartile_current <- initial_params[1]
  median_current <- initial_params[2]
  shift_current <- initial_params[3]
  asymptote_current <- initial_params[4]


  num_rejections <- 0

  # Set up an object to record the results for one chain
  out <- list(
    median_samples = numeric(n_iter),
    shift_samples = numeric(n_iter),
    first_quartile_samples = numeric(n_iter),
    asymptote_samples = numeric(n_iter),
    loglikelihood_current = numeric(n_iter),
    loglikelihood_proposal = numeric(n_iter),
    acceptance_ratio = numeric(n_iter),
    rejection_rate = numeric(n_iter)
  )

  cat("Starting Chain", chain_id, "\n")

  for (i in 1:n_iter) {
    # Propose new values using the prior distributions
    # generate asymptote parameter (gamma)
    repeat {
      asymptote_factor <- 2 * max_penetrance 
      if (asymptote_factor > 1) {
        asymptote_factor <- 1 - SEER_baseline$total_prob
      }
      asymptote_proposal <- SEER_baseline$total_prob +
        do.call(prior_distributions$asymptote_distribution, list(1)) * asymptote_factor
      asymptote_proposal <- max(0, min(1, asymptote_proposal))
      shift_proposal <- do.call(prior_distributions$shift_distribution, list(1))
      median_proposal <- if (median_max) {
        do.call(prior_distributions$median_distribution, list(1)) * (baseline_mid - shift_proposal) + shift_proposal
      } else {
        do.call(prior_distributions$median_distribution, list(1)) * (max_age - shift_proposal) + shift_proposal
      }
      first_quartile_proposal <- do.call(prior_distributions$first_quartile_distribution, list(1)) *
        (median_proposal - shift_proposal) + shift_proposal

      if (validate_weibull_parameters(first_quartile_proposal, median_proposal, shift_proposal, asymptote_proposal)) {
        break
      }
    }

    # Compute the likelihood for the current and proposed
    loglikelihood_current <- mhLogLikelihood_clipp(
      paras = c(
        median_current,
        first_quartile_current, shift_current, 
        asymptote_current
      ), families = data,
      max_age = max_age,
      cancer_type = cancer_type,
      db,
      af = 0.001
    )

    loglikelihood_proposal <- mhLogLikelihood_clipp(
      paras =
        c(
          median_proposal, first_quartile_proposal,
           shift_proposal, asymptote_proposal
        ), families = data,
      max_age = max_age,
      cancer_type = cancer_type,
      db,
      af = 0.001
    )
  

    # Compute the acceptance ratio (likelihood ratio)
    acceptance_ratio <- exp(loglikelihood_proposal - loglikelihood_current)

    # Accept or reject the proposal
    if (runif(1) < acceptance_ratio) {
      median_current <- median_proposal
      shift_current <- shift_proposal
      first_quartile_current <- first_quartile_proposal
      asymptote_current <- asymptote_proposal
    } else {
      num_rejections <- num_rejections + 1 # Increment the rejection counter
    }

    # Update the outputs
    out$median_samples[i] <- median_current
    out$shift_samples[i] <- shift_current
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
#' @param data A list of family data for penetrance estimation.
#' @param cancer_type The type of cancer for which to estimate penetrance.
#' @param gene_input Gene information used for risk estimation.
#' @param n_chains Number of chains for parallel computation.
#' @param n_iter_per_chain Number of iterations for each chain.
#' @param max_age Maximum age considered for analysis, default is 94.
#' @param removeProband Logical, indicating whether to remove probands from the analysis (default is FALSE).
#' @param median_max Boolean indicating whether to use SEER median or max_age as an upper bound for the median proposal. Defaults to TRUE, i.e. using the SEER median.
#' @param burn_in Fraction of results to discard as burn-in (0 to 1). Default is 0 (no burn-in).
#' @param thinning_factor Factor by which to thin the results, default is 1 (no thinning).
#' @param distribution_data Data for generating prior distributions.
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
#' @importFrom parallel makeCluster stopCluster parLapply
#' @examples
#' # Example usage:
#' result <- PenEstim(
#'   data = familyData, cancer_type = "breast", gene_input = "BRCA1",
#'   n_chains = 4, n_iter_per_chain = 1000, max_age = 90,
#'   burn_in = 0.1, thinning_factor = 2, summary_stats = TRUE,
#'   rejection_rates = TRUE, density_plots = TRUE, penetrance_plot = TRUE
#' )
#' @export


PenEstim <- function(data, cancer_type, gene_input, n_chains = 4,
                     n_iter_per_chain = 200,
                     max_age = 94,
                     removeProband = FALSE,
                     median_max = TRUE,
                     burn_in = 0,
                     thinning_factor = 1,
                     distribution_data = distribution_data_default,
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
  if (missing(cancer_type)) {
    stop("Error: 'cancer_type' parameter is missing. Please specify the type of cancer for penetrance estimation.")
  }
  if (missing(gene_input)) {
    stop("Error: 'gene_input' parameter is missing. Please provide the gene.")
  }
  if (missing(n_chains) || !is.numeric(n_chains) || n_chains <= 0) {
    stop("Error: 'n_chains' parameter is missing or invalid.
    Please specify the number of chains to be run (on seperate cores).
    It must be a positive integer.")
  }
  if (missing(n_iter_per_chain) || !is.numeric(n_iter_per_chain) || n_iter_per_chain <= 0) {
    stop("Error: 'n_iter_per_chain' parameter is missing or invalid. It must be a positive integer.")
  }

  # create the seeds for the individual chains
  seeds <- sample.int(1000, n_chains)

  # Apply the prepAges function to treat missing data
  data <- prepAges(data, removeProband)

   # Apply the transformation to adjust the format for the clipp package 
  data <- do.call(rbind, lapply(data, transformDF))

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

  parallel::clusterEvalQ(cl, {
    library(PPP) # Load the "PPP" library
    library(clipp)
    library(stats4)
    library(dplyr)
  })

  parallel::clusterExport(cl, c(
    "mhChain", "mhLogLikelihood_clipp", "calculate_lifetime_risk", 
    "calculate_weibull_parameters", "validate_weibull_parameters", "calculateBaseline",
    "penet.fn", "transformDF", "trans_monogenic2",
    "makePriors",
    "seeds", "n_iter_per_chain",
    "data", "prop", "max_age",
    "PanelPRODatabase", "cancer_type", "gene_input"
  ), envir = environment())

  results <- parallel::parLapply(cl, 1:n_chains, function(i) {
    mhChain(seeds[i],
      n_iter = n_iter_per_chain, chain_id = i,
      data = data,
      db = PanelPRODatabase,
      prior_distributions = prop,
      max_age = max_age,
      cancer_type = cancer_type,
      gene_input = gene_input,
      max_penetrance = max_penetrance
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
        output$penetrance_plot <- plot_penetrance(combined_chains, probCI)
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
