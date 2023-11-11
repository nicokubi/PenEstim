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

mhChain <- function(
  seed, n_iter, chain_id, data,
  max_age, PanelPRODatabase, proposal_distributions, cancer_type, gene_input
) {

  # set seed
  set.seed(seed)

  # Identify the index where cumulative probability crosses the midpoint
  SEER_baseline <- calculate_lifetime_risk(cancer = cancer_type, gene = "SEER")
  midpoint_prob <- SEER_baseline$total_probability/ 2
  midpoint_index <- which(SEER_baseline$cumulative_risk >= midpoint_prob)[1]

  # Identify the age at which the cumulative probability crosses the midpoint
  baseline_mid <- as.numeric(names(lifetime_risk_cum)[midpoint_index])

  # Initialize parameters using random draws from the proposal distributions
  asymptote_start <- total_prob +
    proposal_distributions$asymptote_distribution(1) * (1 - total_prob)
  shift_start <- proposal_distributions$shift_distribution(1)
  median_start <- proposal_distributions$median_distribution(1) *
    (baseline_mid + 5 - shift_start) + shift_start
  first_quartile_start <-
    proposal_distributions$first_quartile_distribution(1) *
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
    # generate aysmptote parameter (gamma)
    asymptote_proposal <- total_prob +
      proposal_distributions$asymptote_distribution(1) * (1 - total_prob)

    # generate shift parameter (delta)
    shift_proposal <- proposal_distributions$shift_distribution(1)

    # generate median
    median_proposal <- proposal_distributions$median_distribution(1) *
      (baseline_mid + 5 - shift_proposal) + shift_proposal

    # generate first quartile
    first_quartile_proposal <-
      proposal_distributions$first_quartile_distribution(1) *
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

PenEstim <- function(data, cancer_type, gene_input, n_chains,
                     n_iter_per_chain,
                     proposal_distributions,
                     max_age = 94,
                     summary_stats = TRUE,
                     rejection_rates = TRUE,
                     density_plots = TRUE,
                     trace_plots = TRUE,
                     burn_in = 0,
                     thinning_factor = 1) {
  seeds <- sample.int(1000, n_chains)

  cl <- parallel::makeCluster(n_chains)

  clusterEvalQ(cl, {
    library(PPP) # Load the "PPP" library
  })

  clusterExport(cl, c(
    "mhChain", "mhLogLikelihood", "seeds", "n_iter_per_chain",
    "data", "proposal_distributions", "max_age",
    "PanelPRODatabase", "cancer_type", "gene_input"
  ), envir = environment())

  results <- parallel::parLapply(cl, 1:n_chains, function(i) {
    mhChain(seeds[i],
      n_iter = n_iter_per_chain, chain_id = i,
      data = data,
      PanelPRODatabase = PanelPRODatabase,
      proposal_distributions = proposal_distributions,
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
