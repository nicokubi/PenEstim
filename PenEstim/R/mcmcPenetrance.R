#' Execution of a Single Chain in Metropolis-Hastings
#'
#' @param seed Seed value for random number generation.
#' @param n_iter Number of iterations for the chain.
#' @param chain_id Identifier for the chain.
#' @param data List of families data.
#' @param proposal_dist List of the statistical distribution for the proposal distributions for the median, first quartile, shift and asymptote parameters. 
#' @param proposal_params List of the parameters for the distributions of the proposal. 
#' @param max_age Maximum age to be considered. 
#' @return A list with samples and rejection rate.
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport clusterEvalQ
#' @importFrom PPP PPP

mhChain <- function(seed, n_iter, chain_id, data,
                    proposal_dist, proposal_params, max_age, PanelPRODatabase) {

  set.seed(seed)

 # Recover the SEER lifetime risk for the cancer 
  gene <- "SEER"
  cancer <- "Breast"
  race <- "All_Races"
  female <- "Female"
  male <- "Male"
  type <- "Crude"
  
  # Find the indices for the resp. attributes
  dim_names <- attr(PanelPRODatabase$Penetrance, "dimnames")
  gene_index <- which(dim_names$Gene == gene)
  cancer_index <- which(dim_names$Cancer == cancer)
  race_index <- which(dim_names$Race == race)
  sex_index <- which(dim_names$Sex == female)
  type_index <- which(dim_names$PenetType == type)

  # Calculate the cummunlative risk for every age up until max. age 
  lifetime_risk_cum <- cumsum(PanelPRODatabase$Penetrance[cancer_index, gene_index, race_index, sex_index, ,type_index])
  total_prob <- sum(lifetime_risk)
  midpoint_prob <- total_prob / 2

  # Identify the index where cumulative probability crosses the midpoint
  midpoint_index <- which(lifetime_risk_cum >= midpoint_prob)[1]

  # Identify the age at which the cumulative probability crosses the midpoint
  baseline_mid <- as.numeric(names(lifetime_risk_cum)[midpoint_index]

  # Initialize parameters using random draws from the proposal distributions
  shift_start <- proposal_dist$shift(proposal_params$shift)
  median_start <- proposal_dist$median_(proposal_params$median)
  first_quartile_start <- proposal_dist$quartile(proposal_params$quartile)
  asymptote_start <- proposal_dist$asymptote(proposal_params$asymptote)
  

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
    asymptote_proposal <- proposal_dist$asymptote(proposal_params$asymptote)
    shift_proposal <- proposal_dist$shift(proposal_params$shift)
    median_proposal <- proposal_dist$median(proposal_params$median)
    first_quartile_proposal <- proposal_dist$quartile(proposal_params$quartile)
    
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
#' @param data List of families data.
#' @param n_chains Number of chains for parallel computation.
#' @param n_iter_per_chain Number of iterations for each chain.
#' @param proposal_dist List of the statistical distribution for the proposal distributions for the median, first quartile, shift and asymptote parameters. 
#' @param proposal_params List of the parameters for the distributions of the proposal. 
#' @param max_age Maximum age to be considered.
#' @return A list containing results for each chain.
#'
#' @importFrom stats rbeta runif dweibull
#' @importFrom parallel makeCluster stopCluster parLapply
#' @importFrom PPP PPP

#' @export

PenEstim <- function(data, n_chains, n_iter_per_chain,
                     max_age=94, proposal_dist, proposal_params){

  seeds <- sample.int(1000, n_chains)

  cl <- makeCluster(n_chains)

  clusterEvalQ(cl, {
    library(PPP)
    PanelPRODatabase
  })

  clusterExport(cl, c("mhChain", "mhLogLikelihood", "seeds", "n_iter_per_chain",
                      "data", "proposal_dist", "proposal_params"
                      "max_age"), envir=environment())

  results <- parLapply(cl, 1:n_chains, function(i) {
    mhChain(seeds[i], n_iter = n_iter_per_chain, chain_id = i,
            data = data,
            PanelPRODatabase = PanelPRODatabase,
            proposal_dist = proposal_dist,
            proposal_params = proposal_params,
            max_age = max_age)
  })

  stopCluster(cl)

  # Check rejection rates and issue a warning if they are all above 90%
  all_high_rejections <- all(sapply(results, function(x) x$rejection_rate > 0.9))
  if(all_high_rejections) {
    warning("Low acceptance rate. Please consider running the chain longer.")
  }

  return(results)}
