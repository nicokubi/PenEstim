#' Execution of a Single Chain in Metropolis-Hastings
#'
#' @param seed Seed value for random number generation.
#' @param n_iter Number of iterations for the chain.
#' @param chain_id Identifier for the chain.
#' @param data List of families data.
#' @param proposal_params List of the parameters for the distributions of the proposal.
#' @param max_age Maximum age to be considered.
#' @return A list with samples and rejection rate.
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport clusterEvalQ
#' @importFrom PPP PPP

# Main mhChain function
mhChain <- function(
    seed, n_iter, chain_id, data,
    max_age, PanelPRODatabase, proposal_distributions, cancer_type, gene_input,
    median_max = TRUE) {
  # Set the right seed
  set.seed(seed)

  # Calculate SEER baseline and midpoint
  SEER_baseline <- calculate_lifetime_risk(cancer = cancer_type, gene = "SEER")
  midpoint_prob <- SEER_baseline$total_prob / 2
  midpoint_index <- which(SEER_baseline$cumulative_risk >= midpoint_prob)[1]
  baseline_mid <- as.numeric(names(SEER_baseline$cumulative_risk)[midpoint_index])

  # Initialize parameters using random draws from the proposal distributions
  asymptote_start <- SEER_baseline$total_prob +
    do.call(proposal_distributions$asymptote_distribution, list(1)) * (1 - SEER_baseline$total_prob)
  shift_start <- do.call(proposal_distributions$shift_distribution, list(1))
  # Initialize median_start with an if-statement to choose between baseline_mid and max_age
  if (median_max == TRUE) {
    median_start <- do.call(proposal_distributions$median_distribution, list(1)) *
      (baseline_mid - shift_start) + shift_start
  } else {
    median_start <- do.call(proposal_distributions$median_distribution, list(1)) *
      (max_age - shift_start) + shift_start
  }
  first_quartile_start <- do.call(proposal_distributions$first_quartile_distribution, list(1)) *
    (median_start - shift_start) + shift_start

  # Initialize parameters using the provided starting values
  median_current <- median_start
  asymptote_current <- asymptote_start
  shift_current <- shift_start
  first_quartile_current <- first_quartile_start

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
    asymptote_proposal <- SEER_baseline$total_prob +
      do.call(proposal_distributions$asymptote_distribution, list(1)) * (1 - SEER_baseline$total_prob)

    # generate shift parameter (delta)
    shift_proposal <- do.call(proposal_distributions$shift_distribution, list(1))

    # generate median
    if (median_max == TRUE) {
      median_proposal <- do.call(proposal_distributions$median_distribution, list(1)) *
        (baseline_mid - shift_proposal) + shift_proposal
    } else {
      median_proposal <- do.call(proposal_distributions$median_distribution, list(1)) *
        (max_age - shift_proposal) + shift_proposal
    }

    # generate first quartile
    first_quartile_proposal <-
      do.call(proposal_distributions$first_quartile_distribution, list(1)) *
      (median_proposal - shift_proposal) + shift_proposal

    # Compute the likelihood for the current and proposed
    loglikelihood_current <- mhLogLikelihood(
      paras = c(
        median_current,
        first_quartile_current, asymptote_current,
        shift_current
      ), families = data,
      max_age = max_age,
      gene_input = gene_input,
      cancer_type = cancer_type,
      PanelPRODatabase = PanelPRODatabase
    )

    loglikelihood_proposal <- mhLogLikelihood(
      paras =
        c(
          median_proposal, first_quartile_proposal,
          asymptote_proposal, shift_proposal
        ), families = data,
      max_age = max_age,
      gene_input = gene_input,
      cancer_type = cancer_type,
      PanelPRODatabase = PanelPRODatabase
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
#' This function employs a Bayesian approach for penetrance estimation, utilizing the
#' Independent Metropolis-Hastings algorithm. It leverages parallel computing and requires
#' the `stats`, `parallel`, and `PPP` packages.
#'
#' @param data List of families data.
#' @param cancer_type The type of cancer to estimate penetrance for.
#' @param n_chains Number of chains for parallel computation.
#' @param n_iter_per_chain Number of iterations for each chain.
#' @param proposal_params List of the parameters for the distributions of the proposal.
#' @param burn_in The fraction proportion of results to discard as burn-in (0 to 1). The default is no burn-in, burn_in=0.
#' @param thinning_factor The factor by which to thin the results (positive integer). The default thinning factor is 1, which implies no thinning.
#' @param max_age Maximum age to be considered. Default is 94, based on PanelPRO settings.
#' @param summary_stats Includes summary statistics in the function output.
#' @param rejection_rates Includes the rejection rates for each chain in the function output.
#' @param density_plots Includes simple density plots for the posterior samples in the function output.
#' @param trace_plots Includes the trace plots for the individual chains in the function output.
#' @return A list containing results for each chain.
#' @importFrom stats rbeta runif dweibull
#' @importFrom parallel makeCluster stopCluster parLapply
#' @importFrom PPP PPP
#'
#' @export

PenEstim <- function(data, cancer_type, gene_input, n_chains = 4,
                     n_iter_per_chain = 200,
                     max_age = 94,
                     summary_stats = TRUE,
                     rejection_rates = TRUE,
                     density_plots = TRUE,
                     trace_plots = TRUE,
                     burn_in = 0,
                     thinning_factor = 1,
                     distribution_data = distribution_data_default,
                     sample_size = NULL,
                     ratio = NULL,
                     proposal_params = proposal_params_default,
                     risk_proportion = risk_proportion_default
                     )  {
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

 # Create the proposal distributions
 prop <- create_distributions(
   data = distribution_data,
   sample_size = sample_size,
   cancer = cancer_type,
   ratio = ratio,
   proposal_params = proposal_params,
   risk_proportion = risk_proportion
 )

  cl <- parallel::makeCluster(n_chains)

  clusterEvalQ(cl, {
    library(PPP) # Load the "PPP" library
  })

  clusterExport(cl, c(
    "mhChain", "mhLogLikelihood", "calculate_lifetime_risk",
    "generate_proposal",
    "seeds", "n_iter_per_chain",
    "data", "prop", "max_age",
    "PanelPRODatabase", "cancer_type", "gene_input"
  ), envir = environment())

  results <- parallel::parLapply(cl, 1:n_chains, function(i) {
    mhChain(seeds[i],
      n_iter = n_iter_per_chain, chain_id = i,
      data = data,
      PanelPRODatabase = PanelPRODatabase,
      proposal_distributions = prop,
      max_age = max_age,
      cancer_type = cancer_type,
      gene_input = gene_input
    )
  })

  parallel::stopCluster(cl)

  # Check rejection rates and issue a warning if they are all above 90%
  all_high_rejections <- all(sapply(results, function(x) x$rejection_rate > 0.9))
  if (all_high_rejections) {
    warning("Low acceptance rate. Please consider running the chain longer.")
  }

  # Apply burn-in and thinning (assuming you have these functions defined)
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

  if (trace_plots) {
    # Generate trace plots
    output$plot_trace <- plot_trace(results, n_chains)
  }

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

  output$combined_chains <- combined_chains
  output$results <- results

  return(output)
}
