#' Execution of a Single Chain in Metropolis-Hastings
#'
#' @param seed Seed value for random number generation.
#' @param n_iter Number of iterations for the chain.
#' @param chain_id Identifier for the chain.
#' @param data List of families data.
#' @param m1 parameter for the beta distribution for the median.
#' @param m2 parameter for the beta distribution fo the median.
#' @param max_age Maximum age to be considered.
#' @param shift_prior_min Minimum possible value for the shift parameter.
#' @param shift_prior_max Maximum possible value for the shift parameter.
#' @param p0 baseline lifetime risk.
#' @param q1 parameter of the beta distribution for the first quartile.
#' @param q2 parameter of the beta distribution for the first quartile.
#' @param g1 parameter of the beta distribution for the shift.
#' @param g2 parameter of the beta distribution for the shift.
#' @return A list with samples and rejection rate.
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport clusterEvalQ
#' @importFrom PPP PPP


mhChain <- function(seed, n_iter, chain_id, data, save_interval,
                    m1, m2, max_age, shift_prior_min, shift_prior_max,
                    p0, q1, q2, g1, g2,PanelPRODatabase) {

  set.seed(seed)

  # Initialize parameters using random draws from the proposal distributions
  median_start <- 60
  asymptote_start <- 0.8
  shift_start <- 24
  first_quartile_start <- 50


  # Initialize parameters using the provided starting values
  median_current <- median_start
  asymptote_current <- asymptote_start
  shift_current <- shift_start
  first_quartile_current <- first_quartile_start

  # Set up empty vectors
  median_samples <- numeric(n_iter)
  first_quartile_samples <- numeric(n_iter)
  asymptote_samples <- numeric(n_iter)
  shift_samples <- numeric(n_iter)

  num_rejections <- 0

  cat("Starting Chain", chain_id, "\n")

  for (i in 1:n_iter) {
    # Propose new values using the prior distributions
    # generate aysmptote parameter (gamma)
    asymptote_proposal <- rbeta(1,g1,g2)
    asymptote_proposal <- p0 + asymptote_proposal *(1-p0)

    # generate shift parameter (delta)
    shift_proposal <- runif(1, shift_prior_min, shift_prior_max)

    # generate median
    median_proposal <- rbeta(1,m1,m2)
    median_proposal <- (median_proposal)*(max_age-shift_proposal) + shift_proposal

    # generate first quartile
    first_quartile_proposal <- rbeta(1,q1,q2)
    first_quartile_proposal <- (first_quartile_proposal)*(median_proposal-shift_proposal) + shift_proposal

    # Compute the likelihood for the current and proposed
    loglikelihood_current <- mhLogLikelihood(paras = c(median_current,first_quartile_current,asymptote_current,
                                                       shift_current), families  = data,
                                             max_age = max_age,PanelPRODatabase = PanelPRODatabase)
    loglikelihood_proposal <- mhLogLikelihood(paras = c(median_proposal,first_quartile_proposal,
                                                        asymptote_proposal,shift_proposal),
                                              families = data,max_age = max_age,PanelPRODatabase = PanelPRODatabase)

    # Compute the acceptance ratio (likelihood ratio)
    acceptance_ratio <- exp(loglikelihood_proposal - loglikelihood_current)

    # Accept or reject the proposal
    if (runif(1) < acceptance_ratio) {
      median_current <- median_proposal
      shift_current <- shift_proposal
      first_quartile_current <- first_quartile_proposal
      asymptote_current <- asymptote_proposal
    } else {
      num_rejections <- num_rejections + 1  # Increment the rejection counter
    }

    # Store the samples
    median_samples[i] <- median_current
    shift_samples[i] <- shift_current
    first_quartile_samples[i] <- first_quartile_current
    asymptote_samples[i] <- asymptote_current

    if (i %% save_interval == 0) {
      save_file_name <- paste0("MH_chain_state_", chain_id, "_iter_", i, ".RDS")
      save_state <- list(median_samples = median_samples[1:i],
                         shift_samples = shift_samples[1:i],
                         first_quartile_samples = first_quartile_samples[1:i],
                         asymptote_samples = asymptote_samples[1:i])
      saveRDS(save_state, save_file_name)
      cat("Saved state at iteration", i, "\n")
    }

  }
  # Return the result as a list
  list(median_samples = median_samples,
       shift_samples = shift_samples,
       first_quartile_samples = first_quartile_samples,
       asymptote_samples = asymptote_samples,
       rejection_rate = num_rejections / n_iter)
}

#' Bayesian Inference using Independent Metropolis-Hastings for Penetrance Estimation
#'
#' This function employs a Bayesian approach for penetrance estimation, utilizing the
#' Independent Metropolis-Hastings algorithm. It leverages parallel computing and requires
#' the `stats`, `parallel`, and `PPP` packages.
#'
#' @param n_chains Number of chains for parallel computation.
#' @param n_iter_per_chain Number of iterations for each chain.
#' @param save_interval Interval after which states should be saved.
#' @param m1 parameter for the beta distribution for the median.
#' @param m2 parameter for the beta distribution fo the median.
#' @param max_age Maximum age to be considered.
#' @param shift_prior_min Minimum possible value for the shift parameter.
#' @param shift_prior_max Maximum possible value for the shift parameter.
#' @param p0 baseline lifetime risk.
#' @param q1 parameter of the beta distribution for the first quartile.
#' @param q2 parameter of the beta distribution for the first quartile.
#' @param g1 parameter of the beta distribution for the shift.
#' @param g2 parameter of the beta distribution for the shift.
#' @param data List of families data.
#' @return A list containing results for each chain.
#'
#' @importFrom stats rbeta runif dweibull
#' @importFrom parallel makeCluster stopCluster parLapply
#' @importFrom PPP PPP

#' @export

PenEstim <- function(data,n_chains, n_iter_per_chain,save_interval=200,
                      max_age=94, shift_prior_min=0, shift_prior_max=25,
                     p0=0.15, m1=2, m2=2,q1=6, q2=3, g1=9, g2=1) {

  seeds <- sample.int(1000, n_chains)



  cl <- makeCluster(n_chains)

  clusterEvalQ(cl, {
    library(PPP)
    PanelPRODatabase
  })

  clusterExport(cl, c("mhChain", "mhLogLikelihood", "seeds", "n_iter_per_chain",
                      "data", "save_interval",
                      "m1", "m2", "max_age", "shift_prior_min", "shift_prior_max",
                      "p0", "q1", "q2", "g1", "g2"), envir=environment())

  results <- parLapply(cl,1:n_chains, function(i) {
    mhChain(seeds[i], n_iter = n_iter_per_chain, chain_id = i,
            data = data,
            PanelPRODatabase = PanelPRODatabase,
            save_interval = save_interval,
            m1 = m1, m2 = m2, max_age = max_age,
            shift_prior_min = shift_prior_min, shift_prior_max = shift_prior_max,
            p0 = p0, q1 = q1, q2 = q2, g1 = g1, g2 = g2)
  })

  stopCluster(cl)

  # Check rejection rates and issue a warning if they are all above 90%
  all_high_rejections <- all(sapply(results, function(x) x$rejection_rate > 0.9))
  if(all_high_rejections) {
    warning("Low acceptance rate. Please consider running the chain longer.")
  }

  return(results)}
