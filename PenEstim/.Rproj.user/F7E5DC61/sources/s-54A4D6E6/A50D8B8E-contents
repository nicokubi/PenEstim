#' Calculate and format carrier probabilities and future risks 
#'
#' Estimates the posterior carrier probabilities using the peeling-paring 
#' algorithm and calculates future risk predictions for the proband(s). If 
#' missing ages were imputed, the average and range of the results from these 
#' imputations are reported. 
#' 
#' @param ped A checked pedigree data frame returned by \code{\link{checkFam}}. 
#' @param lm A list of imputed ages returned by \code{\link{checkFam}}. 
#' @param db A model-specific database returned by \code{\link{buildDatabase}}.
#' @param proband A numeric value or vector of the unique IDs in `ped` 
#' corresponding to the proband(s). 
#' @param ifmod A logical value indicating if cancer penetrance risk modifiers 
#' should be considered. The default is `FALSE`. 
#' @param net A logical value indicating whether to return net or crude future 
#' risk estimates. The default is `FALSE`, in which case crude future risk 
#' estimates will be returned. 
#' @param age.by The fixed age interval for the future risk calculation. The 
#' default is `5`, in which case if the proband's current age is 30, risk 
#' predictions will be calculated at ages 35, 40, 45, etc. 
#' @param max.mut The maximum number of simultaneous mutations allowed. The 
#' default is `2`. 
#' @param debug A logical value indicating whether to debug using browsers and 
#' sequential imputation. The default is `FALSE`. 
#' @param rr.bcrat A data frame with 3 columns: the IDs of the probands (`ID`),
#' the probands' BCRAT relative risks for age < 50 (`rr1`), and the probands' 
#' BCRAT relative risks for age >= 50 (`rr2`). When `net = FALSE` and breast 
#' cancer is included in the model, these values will modify the future risk 
#' calculation. The default is `NULL`, in which case the BCRAT relative risks 
#' will not be incorporated. 
#' @param rr.pop A data frame of race-specific for estimates of 1 /
#' (1 - population attributable risk). The default is `rr.ref`,
#' which was estimated from NHIS 2015.
#' 
#' @return A list with two components: 
#' * `posterior.prob`: A list where each element corresponds to a proband in 
#' `proband` and consists of a data frame of posterior carrier probabilities. 
#' There are three columns: `estimate` (the estimate if no ages were imputed, 
#' or the estimate averaged over all age imputations), `lower` (`NA` if no ages 
#' were imputed, or the minimum result for all age imputations), and `upper` 
#' (`NA` if no ages were imputed, or the maximum result for all age 
#' imputations). The rows correspond to the different genotypes. If no genes 
#' were in the model specification, each proband's data frame of carrier 
#' probabilities is replaced by the message 
#' `"No carrier probabilities were requested by the model specification."` 
#' If ages were imputed, `ImputeRange` messages giving the lower and upper
#' bounds for the probability of carrying any pathogenic variant will be 
#' printed for each proband. 
#' A `ZeroProb` warning is raised if one or more of the estimated 
#' probabilities is `0` (or `1`), germline testing results were incorporated, 
#' and default sensitivities/specifities of `1` were used. 
#' * `future.risk`: A list where each element corresponds to a proband in 
#' `proband` and consists of a list of data frames of future risk estimates at 
#' fixed intervals for each cancer in the model. Each data frame has three 
#' columns: `estimate` (the estimate if no ages were imputed, or the estimate 
#' averaged over all age imputations), `lower` (`NA` if no ages were imputed, 
#' or the minimum result for all age imputations), and `upper` (`NA` if no 
#' ages were imputed, or the maximum result for all age imputations). Each row 
#' corresponds to the risk estimate for at a different age, e.g. 35, 40, 45, 
#' etc. when `age.by = 5` for a proband whose current age is 30. If no cancers 
#' were in the model specification, each proband's list of future risk 
#' estimates is replaced by the message 
#' `"No future risk estimates were requested by the model specification."`
#' 
#' @import foreach
#' @import Rcpp
PanelPROCalc <- function(ped, lm, db, proband, ifmod = FALSE, net = FALSE, 
                         age.by = 5, max.mut = 2, 
                         debug = FALSE, rr.bcrat = NULL, rr.pop = rr.ref) {
  
  # Get model specification from database
  MS_cancers <- db$MS$CANCERS
  MS_homozygous <- db$MS$HOMOZYGOUS
  MS_all_gene_variants <- db$MS$ALL_GENE_VARIANTS
  
  # Indices for genes with heterozygous and homozygous mutations
  if (length(MS_homozygous) == 0) {
    # Empty list if there are no genes with heterozygous and homozygous mutations
    homo_idx <- list()
  } else {
    # Otherwise, identify the index pairs for heterozygous and homozygous mutations
    homo_idx <- lapply(MS_homozygous, function(g) {
      which(MS_all_gene_variants %in% g)
    })
  }
  
  # Generate individual-specific allele frequency matrix
  af <- aperm(db$af, c("Ancestry", "Gene"))[ped$Ancestry, , drop = FALSE]
  
  # Define possible genotypes
  if (length(MS_all_gene_variants) > 0) {
    PGs <- .getPossibleGenotype(MS_all_gene_variants, max_mut = max.mut, 
                                homo_genes = db$MS$HOMOZYGOUS)
  } else {
    PGs = list(df = as.data.frame(0), list = "SEER")
  }
  
  if (length(MS_cancers) > 0) { # If at least one cancer is specified
    # Extract genotypes with no more than 1 mutation
    direct_fill_PGs <- unname(PGs$list[rowSums(PGs$df) < 2])
    # Extract genotypes with 2 or more mutations
    # Character strings used for naming
    multi_PGs <- unname(PGs$list[rowSums(PGs$df) >= 2]) 
    # In list format
    multi_muts <- strsplit(PGs$list, split = "\\.")[rowSums(PGs$df) >= 2] 
    
    # Subset net penetrances for likelihood calculation
    sub_dens_lik = subsetCancerPenetrance(fam = ped, lm.ages = lm, db = db, 
                                          proband = NULL, 
                                          consider.modification = ifmod,
                                          net = TRUE, doc = FALSE)
    
    # Cancer penetrance for probands (used for calculating future risk)
    if (net == TRUE) { 
      # If net future risk estimates are requested, subset sub_dens_lik 
      # when calculating the penetrances
      cplist.fr <- calcCancerPenetrance(
        ped, proband, 
        sub_dens = sub_dens_lik[as.character(proband),,,,drop = FALSE], 
        PGs = PGs, 
        direct_fill_PGs = direct_fill_PGs, 
        multi_PGs = multi_PGs, 
        multi_muts = multi_muts, 
        net = TRUE,
        consider.modification = ifmod, 
        doc = FALSE,
        lm.ages = lm
      )
    } else { 
      # Crude penetrances otherwise
      cplist.fr <- calcCancerPenetrance(ped, proband, db, PGs = PGs, 
                                        direct_fill_PGs = direct_fill_PGs, 
                                        multi_PGs = multi_PGs, 
                                        multi_muts = multi_muts, 
                                        net = FALSE,
                                        consider.modification = ifmod, 
                                        doc = FALSE,
                                        lm.ages = lm
      )
    }
    # Death by other causes for each proband
    doclist <- calcCancerPenetrance(ped, proband, db, PGs = PGs, 
                                    direct_fill_PGs = direct_fill_PGs, 
                                    multi_PGs = multi_PGs, 
                                    multi_muts = multi_muts, 
                                    net = FALSE,
                                    consider.modification = FALSE, 
                                    doc = TRUE,
                                    lm.ages = lm
    )
    
  } else { # If no cancers are specified
    sub_dens_lik = NULL
    cplist.fr = NULL
    doclist = NULL
  }
  
  if (length(MS_all_gene_variants) > 0) { # If at least one gene is specified
    # Germline testing used to modify likelihood
    if (ifmod && any(MS_all_gene_variants %in% colnames(ped))) {
      germ <- germlineContrib(ped, db = db, PGs = PGs)
    } else {
      germ <- 1
    }
    
    # Marker testing used to modify likelihood
    if (ifmod) {
      marker <- markerContrib(ped, db = db, PGs = PGs)
    } else {
      marker <- 1
    }
  } else {
    germ = NULL
    marker = NULL
  }
  
  # If the pedigree is complete, no age imputations or averaging are necessary
  if (is.null(lm) || sum(sapply(lm, length)) == 0) {
    
    if (length(MS_all_gene_variants) > 0) { # If at least one gene is specified
      if (length(MS_cancers) > 0) { # If at least one cancer is specified
        
        # Calculate likelihood (using net penetrances)
        lik <- calcLik(ped, sub_dens = sub_dens_lik, PGs = PGs, 
                       direct_fill_PGs = direct_fill_PGs, 
                       multi_PGs = multi_PGs, multi_muts = multi_muts)
      } else { # If no cancers are specified
        
        # Create a dummy likelihood of 1s
        lik = matrix(1, nrow = nrow(ped), ncol = length(PGs$list), 
                     dimnames = list(ped$ID, PGs$list))
      }
      
      # Modify likelihood with germline and marker testing
      lik <- lik * germ * marker
      
      # If pedigree contains twins
      if (any(ped$Twins != 0)) {
        
        # Modify the likelihood for twins
        temp_res <- .twinsLikMod(lik = lik, ped = ped, proband = proband)
        
        # Extract collapsed likelihood, pedigree, probands, and data frame 
        # keeping track of twins
        lik <- temp_res$lik
        collapsed_ped <- temp_res$ped
        collapsed_proband <- temp_res$proband
        twin_labels_df <- temp_res$twin_labels_df
        
        # Run peeling-paring on collapsed likelihood and family
        pmat <- pp.peelingParing(collapsed_ped,
                                 proband = collapsed_proband,
                                 af = af, LIK = lik, 
                                 K = length(MS_all_gene_variants),
                                 M = max.mut, 
                                 homo_idx = homo_idx
        )
        
        # Identify twins who are probands but were dropped
        dropped_twin_ids <- proband[!(proband %in% collapsed_proband)]
        if (length(dropped_twin_ids) > 0) {
          for (dropped_twin in dropped_twin_ids) {
            # Identify the corresponding twin who was kept
            kept_twin <- twin_labels_df$ID[twin_labels_df$isKept == 1 &
                                             twin_labels_df$Label ==
                                             twin_labels_df$Label[twin_labels_df$ID == dropped_twin]]
            
            # Duplicate probability of kept twin
            pmat <- rbind(pmat, pmat[as.character(kept_twin), ])
            rownames(pmat)[nrow(pmat)] <- dropped_twin
            pmat <- pmat[match(proband, rownames(pmat)), ]
          }
        }
      } else {
        # Run peeling-paring
        pmat <- pp.peelingParing(ped,
                                 proband = proband,
                                 af = af, LIK = lik, 
                                 K = length(MS_all_gene_variants),
                                 M = max.mut, 
                                 homo_idx = homo_idx
        )
      }
    } else { # If no genes are specified
      # Create a dummy probability matrix of 1s
      pmat = array(1, dim = c(length(proband), 1), dimnames = list(proband, "SEER"))
    }
    
    # Calculate future risk
    if (length(MS_cancers) > 0) { # If at least one cancer is specified
      frisk <- calcFutureRisk(
        ped = ped, proband = proband,
        MS_cancers = MS_cancers,
        cplist = cplist.fr,
        doclist = doclist,
        posterior_probs = pmat,
        net = net, age.by = age.by,
        rr.bcrat = rr.bcrat, rr.pop = rr.pop,
        cbc_penets = db$penet$penet_cbc, 
        lm.ages = lm
      )
    } else { # If no cancers are specified, return message
      frisk <- replicate(length(proband),
                         "No future risk estimates were requested by the model specification.",
                         simplify=FALSE)
      names(frisk) = proband
    }
    
    # Format posterior probabilities
    if (length(MS_all_gene_variants) > 0) { # If at least one gene is specified
      pmat <- apply(pmat, 1, function(probi) {
        output <- data.frame(genes = colnames(pmat))
        output$estimate <- probi
        output$lower <- NA
        output$upper <- NA
        
        return(output)
      })
      
      # Raise a warning if a zero probability was estimated and germline 
      # testing results were incorporated using sensitivities/specificities 
      # of 1. 
      if (length(germ) > 1 & 
          any(sapply(pmat, function(x) {any(x == 0, na.rm = TRUE)})) & 
          1 %in% db$germline) {
        rlang::warn("Some of the estimated probabilities are 0 or 1, most likely due to the incorporation of reported germline testing results (by default, we assume that the sensitivity and specificity for all germline tests is 1). However, this does not guarantee that the person is a noncarrier or mutation carrier; a more accurate estimate may be achieved with user-specified sensitivities and specificities.", 
                    level = "ZeroProb")
      }
      
    } else { # If no genes are specified, return message
      pmat_dn = dimnames(pmat)[[1]]
      pmat <- replicate(nrow(pmat), 
                        "No carrier probabilities were requested by the model specification.", 
                        simplify=FALSE)
      names(pmat) = pmat_dn
    }
    
    # Get the number of cancers from the lists
    nCancers <- length(dimnames(doclist$Dens)$cancers)
    
    # Get the names of the cancers
    cancers <- dimnames(doclist$Dens)$cancers
    
    # Format future risk
    frisk <- lapply(frisk, function(friski) {
      output <- list()
      for (k in seq_len(nCancers)) {
        if (is.character(friski) == TRUE) {
          output[[cancers[k]]] <- friski
        } else {
          output[[cancers[k]]] <- data.frame(
            ByAge = friski$ByAge,
            estimate = friski[[cancers[k]]],
            lower = NA,
            upper = NA
          )
        }
      }
      return(output)
    })
    
    return(list(posterior.prob = pmat, future.risk = frisk))
  }
  
  # Check that the ages are consistent across imputations
  stopifnot(length(unique(sapply(lm, nrow))) == 1)
  # Number of times missing ages were imputed
  impute.times <- nrow(lm$CurAge)
  
  # Get results for multiple imputations depending on debug or parallel option
  if (!debug) { 
    # This it the default option which takes advantage of parallel imputation
    par_res <- foreach::foreach(
      tt = iterators::icount(impute.times),
      .inorder = FALSE, .export = ".parallelCalc"
    ) %dopar% {
      
      # Return the list of lik, pmat, frisk into par_res
      .parallelCalc(
        tt, lm, ped, db, sub_dens_lik, cplist.fr, 
        PGs, direct_fill_PGs, multi_PGs, multi_muts, 
        germ, marker, proband,
        max.mut, af, net, rr.bcrat, rr.pop,
        doclist, age.by, homo_idx
      )
    }
  } else {
    # Sequential version, especially used for debugging
    par_res <- list()
    for (tt in 1:impute.times) {
      
      # Fill in par_res list
      par_res[[tt]] <- .parallelCalc(
        tt, lm, ped, db, sub_dens_lik, cplist.fr, 
        PGs, direct_fill_PGs, multi_PGs, multi_muts, 
        germ, marker, proband, 
        max.mut, af, net, rr.bcrat, rr.pop,
        doclist, age.by, homo_idx
      )
    }
  }
  
  # Extract results
  pmats <- lapply(par_res, `[[`, 1)
  frisks <- lapply(par_res, `[[`, 2)
  pmats_dn <- dimnames(pmats[[1]])
  
  # Get the number of cancers from the lists
  nCancers <- length(dimnames(doclist$Dens)$cancers)
  # Get the names of the cancers
  cancers <- dimnames(doclist$Dens)$cancers
  
  # Posterior probabilities
  if (length(MS_all_gene_variants) > 0) { # If at least one gene is specified
    # Re-format posterior probability data frames
    if (length(MS_all_gene_variants) > 0) {
      genes <- colnames(pmats[[1]])
    }
    row_num <- nrow(pmats[[1]])
    
    # Format average estimates
    pmats_avg <- matrix(apply(do.call(rbind, lapply(pmats, function(m) {
      vec <- as.vector(m)
    })), 2, mean), nrow = row_num, byrow = FALSE, dimnames = pmats_dn)
    # Format min/lower estimates
    pmats_min <- matrix(apply(do.call(rbind, lapply(pmats, function(m) {
      vec <- as.vector(m)
    })), 2, min), nrow = row_num, byrow = FALSE, dimnames = pmats_dn)
    # Format max/upper estimates
    pmats_max <- matrix(apply(do.call(rbind, lapply(pmats, function(m) {
      vec <- as.vector(m)
    })), 2, max), nrow = row_num, byrow = FALSE, dimnames = pmats_dn)
    
    # Transform pmats results
    pmat_res <- lapply(seq_len(nrow(pmats_avg)), function(i) {
      output <- data.frame(genes = genes)
      output$estimate <- pmats_avg[i, ]
      output$lower <- pmats_min[i, ]
      output$upper <- pmats_max[i, ]
      return(output)
    })
    
    # Print messages with the percent change in the upper/lower bounds for 
    # each proband's carrier probability
    sapply(seq_len(length(pmat_res)), function(i) {
      lower = 1 - pmat_res[[i]][which(pmat_res[[i]][,"genes"]=="noncarrier"),"upper"]
      upper = 1 - pmat_res[[i]][which(pmat_res[[i]][,"genes"]=="noncarrier"),"lower"]
      rlang::inform(sprintf(
        "ID %s's lower and upper bounds for the probability of carrying any pathogenic variant are (%s, %s).", 
        pmats_dn[[1]][i], round(lower, 5), round(upper, 5)), 
        level = "ImputeRange")
    })
    # Print general message explaining that range is due to age imputations
    rlang::inform(
      "The value ranges in the results are driven by the different imputed missing ages at each imputation, so please consider including more age information when possible to reduce the range widths.", 
      level = "ImputeRange")
    
    # Raise a warning if a zero probability was estimated and germline 
    # testing results were incorporated using sensitivities/specificities 
    # of 1. 
    if (length(germ) > 1 & 
        any(sapply(pmat_res, function(x) {any(x == 0, na.rm = TRUE)})) & 
        1 %in% db$germline) {
      rlang::warn("Some of the estimated probabilities are 0 or 1, most likely due to the incorporation of reported germline testing results (by default, we assume that the sensitivity and specificity for all germline tests is 1). However, this does not guarantee that the person is a noncarrier or mutation carrier; a more accurate estimate may be achieved with user-specified sensitivities and specificities.", 
                  level = "ZeroProb")
    }
  } else { # If no genes are specified, return message
    pmat_res <- replicate(length(pmats_dn[[1]]), 
                          "No carrier probabilities were requested by the model specification.", 
                          simplify=FALSE)
  }
  
  # Name the list components after proband IDs
  names(pmat_res) <- pmats_dn[[1]]
  
  # Format future risk
  frisk_aggregate <- lapply(1:length(frisks[[1]]), function(i) {
    lfrisk <- lapply(frisks, `[[`, i)
    if (is.character(lfrisk[[1]])) {
      return(lfrisk[[1]])
    } else {
      frisk <- do.call(rbind, lfrisk)
      frisk_min <- stats::aggregate(frisk, list(frisk$ByAge), min)
      frisk_min <- frisk_min[order(frisk_min$ByAge), -1]
      frisk_avg <- stats::aggregate(frisk, list(frisk$ByAge), mean)
      frisk_avg <- frisk_avg[order(frisk_avg$ByAge), -1]
      frisk_max <- stats::aggregate(frisk, list(frisk$ByAge), max)
      frisk_max <- frisk_max[order(frisk_max$ByAge), -1]
      return(list(
        estimate = frisk_avg,
        lower = frisk_min,
        upper = frisk_max
      ))
    }
  })
  
  if (nCancers > 0) { # If at least one cancer is specified
    frisk_res <- lapply(seq_len(length(frisks[[1]])), function(i) {
      output <- list()
      for (k in seq_len(nCancers)) {
        if (is.character(frisk_aggregate[[i]]) == TRUE) {
          output[[cancers[k]]] <- frisk_aggregate[[i]]
        } else {
          output[[cancers[k]]] <- data.frame(
            ByAge = frisk_aggregate[[i]]$estimate$ByAge,
            estimate = frisk_aggregate[[i]]$estimate[[cancers[k]]],
            lower = frisk_aggregate[[i]]$lower[[cancers[k]]],
            upper = frisk_aggregate[[i]]$upper[[cancers[k]]]
          )
        }
      }
      return(output)
    })
  } else { # If no cancers are specified
    frisk_res <- frisk_aggregate
  }
  
  # Name the list components after proband IDs
  names(frisk_res) <- pmats_dn[[1]]
  
  return(list(
    posterior.prob = pmat_res,
    future.risk = frisk_res
  ))
}


#' Calculate carrier probabilities and future risks for a single age imputation 
#'
#' This function is called within \code{\link{PanelPROCalc}} to calculate 
#' values when multiply imputing missing ages. It wraps the routine that can be 
#' executed in parallel, allowing for flexibility when toggling the `parallel` 
#' argument. It is also useful for debugging (`debug = TRUE`), in which the 
#' non-parallel version should be executed. 
#' 
#' @param tt A numeric value that indexes which set of imputed ages in `lm` is 
#' being used in the calculation. 
#' @param lm A list of imputed ages returned by \code{\link{checkFam}}. 
#' @param ped A checked pedigree data frame returned by \code{\link{checkFam}}. 
#' @param db A model-specific database returned by \code{\link{buildDatabase}}.
#' @param sub_dens_lik A numeric array of subsetted cancer penetrances returned 
#' by \code{\link{subsetCancerPenetrance}}, to be used for likelihood 
#' calculation. 
#' @param cplist.fr A list of density and survival arrays for cancer 
#' penetrances for the probands in `proband`, returned by 
#' \code{\link{calcCancerPenetrance}}.  (If `net == TRUE`, they should be net 
#' penetrances, and if `net == FALSE`, they should be crude.) They will be used 
#' to calculate future risk. 
#' @param PGs Possible genotypes in both list and data frame format, returned 
#' by \code{\link{.getPossibleGenotype}}. 
#' @param direct_fill_PGs A vector of genotype names with no more than 1 
#' mutation (including `"noncarrier"`), corresponding to `PGs`. 
#' @param multi_PGs A vector of genotype names with 2 or more simultaneous 
#' mutations, corresponding to `PGs`.  
#' @param multi_muts A list where each component corresponds to a genotype in 
#' `multi_PGs`, represented as a vector of gene names to indicate the mutations. 
#' mutations, corresponding to `PGs`. 
#' @param germ A matrix of numeric values used to modify the likelihood using 
#' germline testing, returned by \code{\link{germlineContrib}}, or `1`. 
#' @param marker A matrix of numeric values used to modify the likelihood using 
#' marker testing, returned by \code{\link{markerContrib}}, or `1`. 
#' @param proband A numeric value or vector of the unique IDs in `ped` 
#' corresponding to the proband(s). 
#' @param max.mut The maximum number of simultaneous mutations allowed. 
#' @param af A matrix of allele frequencies, where each row corresponds to the 
#' allele frequencies of a relative in `ped`. 
#' @param net A logical value indicating whether to return net or crude future 
#' risk estimates. 
#' @param rr.bcrat A numeric vector of length 2 containing the proband's BCRAT 
#' relative risks for ages < 50 and ages >= 50. 
#' @param rr.pop A data frame of race-specific for estimates of 1 /
#' (1 - population attributable risk). 
#' @param doclist A list of density and survival arrays for death by other 
#' causes for the probands in `proband`, returned by 
#' \code{\link{calcCancerPenetrance}}. They will be used to calculate future 
#' risk. 
#' @param age.by The fixed age interval for the future risk calculation. 
#' @param homo_idx List of indices for genes with heterozygous and homozygous 
#' mutations. Each element consists of a vector with the index pair 
#' corresponding to the heterozygous and homozygous mutations for a given gene. 
#' 
#' @return A list with two components: 
#' * A matrix of posterior carrier probabilities returned by 
#' \code{\link{pp.peelingParing}}. The rows represent the relatives with IDs 
#' in `proband` and the columns represent the genotypes. If no genes were in 
#' the model specification, the values in the matrix are set to `1` for all 
#' probands and genotypes. 
#' * `future.risk`: A list where each element corresponds to a proband in 
#' `proband` and consists of a list of future risk estimates, as returned by 
#' \code{\link{calcFutureRisk}}. If no cancers were in the model specification, 
#' each proband's list of future risk estimates is replaced by the message 
#' `"No future risk estimates were requested by the model specification."`
.parallelCalc <- function(tt, lm, ped, db, sub_dens_lik, cplist.fr, 
                          PGs, direct_fill_PGs, multi_PGs, multi_muts, 
                          germ, marker, proband, 
                          max.mut, af, net, rr.bcrat, rr.pop,
                          doclist, age.by, homo_idx) {
  
  # Get model specification from database
  MS_cancers <- db$MS$CANCERS
  MS_homozygous <- db$MS$HOMOZYGOUS
  MS_all_gene_variants <- db$MS$ALL_GENE_VARIANTS
  
  # Fill in imputed ages
  for (ll in seq_along(lm)) {
    col_to_fill <- names(lm)[ll]
    ids_to_fill <- colnames(lm[[ll]])
    ped[[col_to_fill]][match(ids_to_fill, ped$ID)] <- lm[[ll]][tt, ]
  }
  
  # Check that there are no missing current ages
  stopifnot(all(ped$CurAge > 0))
  
  if (length(MS_all_gene_variants) > 0) { # If at least one gene is specified
    if (length(MS_cancers) > 0) { # If at least one cancer is specified
      # Calculate likelihood (using net penetrances)
      lik <- calcLik(ped, sub_dens = sub_dens_lik, PGs = PGs, 
                     direct_fill_PGs = direct_fill_PGs, 
                     multi_PGs = multi_PGs, multi_muts = multi_muts)
    } else {
      lik = matrix(1, nrow = nrow(ped), ncol = length(PGs$list), 
                   dimnames = list(ped$ID, PGs$list))
    }
    # Modify likelihood with germline and marker testing
    lik <- lik * germ * marker
    
    # If pedigree contains twins
    if (any(ped$Twins != 0)) {
      # Modify the likelihood for twins
      temp_res <- .twinsLikMod(lik = lik, ped = ped, proband = proband)
      
      # Extract collapsed likelihood, pedigree, probands, and data frame 
      # keeping track of twins
      lik <- temp_res$lik
      collapsed_ped <- temp_res$ped
      collapsed_proband <- temp_res$proband
      twin_labels_df <- temp_res$twin_labels_df
      
      # Run peeling-paring on collapsed likelihood and family
      pmat <- pp.peelingParing(collapsed_ped,
                               proband = collapsed_proband,
                               af = af, LIK = lik, 
                               K = length(MS_all_gene_variants),
                               M = max.mut, 
                               homo_idx = homo_idx
      )
      
      # Identify twins who are probands but were dropped
      dropped_twin_ids <- proband[!(proband %in% collapsed_proband)]
      if (length(dropped_twin_ids) > 0) {
        for (dropped_twin in dropped_twin_ids) {
          # Identify the corresponding twin who was kept
          kept_twin <- twin_labels_df$ID[twin_labels_df$isKept == 1 &
                                           twin_labels_df$Label ==
                                           twin_labels_df$Label[twin_labels_df$ID == dropped_twin]]
          
          # Duplicate probability of kept twin
          pmat <- rbind(pmat, pmat[as.character(kept_twin), ])
          rownames(pmat)[nrow(pmat)] <- dropped_twin
          pmat <- pmat[match(proband, rownames(pmat)), ]
        }
      }
    } else {
      # Run peeling-paring
      pmat <- pp.peelingParing(ped,
                               proband = proband,
                               af = af, LIK = lik, 
                               K = length(MS_all_gene_variants),
                               M = max.mut, 
                               homo_idx = homo_idx
      )
    }
  } else {# If no genes are specified
    # Create a dummy probability matrix of 1s
    pmat = array(1, dim = c(length(proband), 1), 
                 dimnames = list(proband, "SEER"))
  }
  
  # Calculate future risk
  if (length(MS_cancers) > 0) { # If at least one cancer is specified
    frisk <- calcFutureRisk(
      ped = ped, proband = proband,
      MS_cancers = MS_cancers,
      cplist = cplist.fr,
      doclist = doclist,
      posterior_probs = pmat,
      net = net, age.by = age.by,
      rr.bcrat = rr.bcrat, rr.pop = rr.pop,
      cbc_penets = db$penets$penet_cbc, 
      lm.ages = lm
    )
  } else { # If no cancers are specified, return message
    frisk <- replicate(length(proband),
                       "No future risk estimates were requested by the model specification.",
                       simplify=FALSE)
    names(frisk) = proband
  }
  
  return(list(pmat, frisk))
}


#' PPP master function
#'
#' This is the main function of the `PPP` package. It combines the 
#' pre-processing steps (checking the user-supplied pedigree, imputing 
#' missing ages, and setting up the database of model parameters) with the 
#' calculation of the genotypic probability distributions and future risks of 
#' cancer for the specified probands. 
#' 
#' @param pedigree A data frame of family history information with the 
#' following columns. Unknown or missing values should be explicitly coded 
#' as `NA`. 
#' * `ID`: A numeric value; ID for each individual. There should not be any 
#' duplicated entries. 
#' * `Sex`: A numeric value; `0` for female and `1` for male. Missing entries 
#' are not currently supported. 
#' * `MotherID`: A numeric value; unique ID for someone's mother. 
#' * `FatherID`: A numeric value; unique ID for someone's father. 
#' * `isProband`: A numeric value; `1` if someone is a proband, `0` otherwise. 
#' This will be overridden by the `proband` argument in `PPP`, if it is 
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
#' * `race`: A character string; expected values are `"All_Races"`, `"AIAN"` 
#' (American Indian and Alaska Native), `"Asian"`, `"Black"`, `"White"`, 
#' `"Hispanic"`, `"WH"` (white Hispanic), and `"WNH"` (non-white Hispanic) 
#' (see `PPP:::RACE_TYPES`). Asian-Pacific Islanders should be encoded as 
#' `"Asian"`. Race information will be used to select the cancer and death by 
#' other causes penetrances used in the model. Missing entries are recoded as 
#' the `unknown.race` argument, which defaults to 
#' `PPP:::UNKNOWN_RACE`. 
#' * `Ancestry`: A character string; expected values are `"AJ"`(Ashkenazi Jewish), 
#' `"nonAJ"`, and `"Italian"` (see `PPP:::ANCESTRY_TYPES`). The ancestry 
#' information is used to determine the allele frequencies of BRCA1 and BRCA2, 
#' if they are included in the model. Missing entries are recoded as the 
#' `unknown.ancestry` argument, which defaults to `PPP:::UNKNOWN_ANCESTRY`. 
#' * `Twins`: A numeric value; `0` for non-identical/single births, `1` for 
#' the first set of identical twins/multiple births in the family, `2` for the 
#' second set, etc. Missing entries are assumed to be `0`. 
#' * `riskmod`: A character list; expected values are `"Mastectomy"`, 
#' `"Hysterectomy"`, and `"Oophorectomy"` 
#' (see `data.frame(PPP:::RISKMODS)`). These preventative interventions 
#' will be used to modify the cancer penetrances for breast, endometrial, and 
#' ovarian cancer, respectively. Currently, `"Mastectomy"` refers to a 
#' bilateral prophylactic mastectomy (occurring prior to first breast cancer). 
#' * `interAge`: A numeric list; the age of intervention for each risk modifier 
#' in `riskmod`. For example, `riskmod = list("Mastectomy", "Hysterectomy")` 
#' and `interAge = {45, 60}` indicates that the individual had a mastectomy at 
#' age 45 and a hysterectomy at age 60.
#' * There can be optional columns for germline testing results (e.g. `BRCA1`, 
#' `MLH1`) or tumor marker testing results. `ER`, `PR`, `CK14`, `CK5.6` and 
#' `HER2` are tumor markers associated with breast cancer that will modify the 
#' likelihood of phenotypes associated with `BRCA1` and `BRCA2`. `MSI` is a 
#' tumor marker for colorectal cancer that will modify the likelihoods 
#' associated with `MLH1`, `MSH2` `MSH6` and `PMS2`. For each of these optional 
#' columns, positive results should be coded as `1`, negative results should be 
#' coded as `0`, and unknown results should be coded as `NA`. 
#' @param model_spec One of the following character strings indicating a 
#' pre-specified model: `"BRCAPRO"`, `"BRCAPRO5"`, `"BRCAPRO6"`, `"MMRPRO"`, 
#' `"PanPRO11"`, and `"PanPRO22"`. See `PPP:::MODELPARAMS` for more 
#' information. Either `model_spec` or `genes` and `cancers` should be 
#' specified by the user. When both input types are supplied, the input from 
#' `genes` and `cancers` will take precedence. The default is `"PanPRO22"`, 
#' the largest pre-specified model currently supported by the package. 
#' @param proband A numeric value or vector of the unique IDs in `pedigree` 
#' corresponding to the proband(s) for whom carrier probabilities and future 
#' risk predictions should be estimated. The default is 
#' `pedigree$ID[pedigree$isProband == 1]`, in which case probands will be 
#' inferred from the `isProband` column in `pedigree`, if it exists and is 
#' valid. Otherwise, the IDs in `proband` will override the contents of 
#' `pedigree$isProband`. The user can also specify `proband = "All"`, in which 
#' case all IDs in the `ped` will be treated as probands. 
#' @param genes A character vector of the genes of interest. The default is 
#' `NULL`. Available options are listed in `PPP:::GENE_TYPES`. Either 
#' `model_spec` or `genes` and `cancers` should be specified by the user. When 
#' both input types are supplied, the input from `genes` and `cancers` will 
#' take precedence. If at least one cancer is specified in `cancers`, but no 
#' genes are specified in `genes`, this will be treated as a valid model that 
#' calculates future risk predictions based on the model parameters in 
#' `database`, without incorporating the family history in `pedigree`, and does 
#' not return any carrier probabilities. 
#' @param cancers A character vector of the cancers of interest. Available 
#' options are listed in `PPP:::CANCER_TYPES`. The default is `NULL`. 
#' Either `model_spec` or `genes` and `cancers` should be specified by the 
#' user. When both input types are supplied, the input from `genes` and 
#' `cancers` will take precedence. If at least one gene is specified in 
#' `genes`, but no cancers are specified in `cancers`, this will be treated as 
#' a valid model that calculates carrier probability estimates based on the 
#' model parameters in `database`, without incorporating the family history in 
#' `pedigree`, and does not return any future risks. 
#' @param database A structured list providing the parameters of the model.
#' The default is \code{\link{PanelPRODatabase}}. Users who wish to provide 
#' their own model parameters should supply a list with the same structure as 
#' `PanelPRODatabase`. 
#' @param unknown.race A character string indicating the default race to use 
#' when `pedigree$race` is missing or unsupported. The default is 
#' `PPP:::UNKNOWN_RACE`. 
#' @param unknown.ancestry A character string indicating the default ancestry 
#' to use when `pedigree$Ancestry` is missing or unsupported. The default is 
#' `PPP:::UNKNOWN_ANCESTRY`. 
#' @param impute.missing.ages A logical value indicating whether missing ages 
#' should be multiply imputed. The default is `TRUE`. If set to `FALSE`, an 
#' error will be raised if any individuals have missing current or cancer 
#' affection ages. 
#' @param allow.intervention A logical value indicating if cancer penetrance 
#' risk modifiers should be considered. The default is `TRUE`. Currently 
#' supported interventions are preventative mastectomies, hysterectomies, and 
#' oophorectomies (see `data.frame(PPP:::RISKMODS)`). 
#' @param ignore.proband.germ A logical value indicating if the proband(s)'s 
#' germline testing results (if provided in `ped`) should be considered. The 
#' default is `FALSE`, in which case the proband(s)'s germline testing results 
#' will not be ignored. 
#' @param max.mut The maximum number of simultaneous mutations allowed. The 
#' default is `NULL`, in which case a maximum of `2` simultaneous mutations 
#' will be used when `pedigree` has 30 or fewer members, and `1` otherwise. 
#' @param remove.miss.cancers A logical value indicating whether cancers in 
#' `cancers` or implied by `model_spec` should be removed from the model 
#' specification when they are missing in `pedigree`. The default is `TRUE`. 
#' @param iterations A numeric value indicating the number of iterations that 
#' should be run when multiply imputing missing ages. The default is `20`.
#' @param max.iter.tries A numeric value indicating the maximum number of 
#' iterations that should be attempted when multiply imputing missing ages. 
#' (Invalid iterations typically occur when an individual is missing both their 
#' current and cancer ages. Current ages for living individuals are imputed 
#' first; if a valid cancer affection cannot be found, subject to the upper 
#' bound set by the imputed current age, this imputation run will be 
#' discarded.) The default is `NULL`, in which case the maximum number of 
#' tries will be set to five times `iterations`. 
#' @param parallel A logical value indicating whether the peeling-paring 
#' algorithm should be run in parallel when missing ages are imputed multiple 
#' times. The default is `TRUE`. 
#' @param debug A logical value indicating whether to debug using browsers and 
#' sequential imputation. The default is `FALSE`. Setting `debug = TRUE` is 
#' only useful for package developers. 
#' @param net A logical value indicating whether to return net or crude future 
#' risk estimates. The default is `FALSE`, in which case crude future risk 
#' estimates will be returned. Net risk is the probability of developing a 
#' given cancer by the specified age, in a world where there is no death from 
#' other causes. Crude risk is the probability of developing a given cancer by 
#' the specified age and not having died from other causes beforehand. 
#' @param age.by The fixed age interval for the future risk calculation. The 
#' default is `5`, in which case if the proband's current age is 30, risk 
#' predictions will be calculated at ages 35, 40, 45, etc. 
#' @param random.seed The random seed (a numeric value) to set for imputing 
#' missing ages. The default is `42L`. 
#' @param plusBCRAT A logical value indicating whether to use the BRCAPRO+BCRAT 
#' model. The default is `FALSE`. If `plusBCRAT = TRUE`, the user should 
#' specify one of `bcrat.vars` or `rr.bcrat`. The BCRAPRO+BCRAT model is only 
#' run if the user also sets `net = FALSE` and breast cancer is included in the 
#' model specification. 
#' @param bcrat.vars A data frame of BCRAT covariates, where each row 
#' corresponds to a proband in `pedigree`. The expected columns are: 
#' * `ID`: Numeric ID for each proband. Should not contain duplicated 
#' entries. 
#' * `Age`: Initial age (should be >=20 and <90). 
#' * `NBiops`: Number of benign breast biopsies (`99` if unknown). 
#' * `Hyp`: A logical value indicating whether the individual has had a breast 
#' biopsy with atypical hyperplasia (`0` if no, `1` if yes, `99` if unknown). 
#' * `AgeMen`: Age at menarche (`99` if unknown)
#' * `AgeFLB`: Age at first live birth (`98` if nulliparous, `99` if unknown). 
#' * `NumRel`: The number of female first-degree relatives with breast cancer 
#' (`99` if unknown). 
#' * `Race`: A character string indicating race (`"White"`, `"Black"`, 
#' `"Hispanic"`, `"FHispanic"`, `"Asian"`, or `"Other"`; `"Hispanic"` 
#' corresponds to U.S.-born Hispanic while `"FHispanic"` corresponds to 
#' foreign-born Hispanic). 
#' 
#' The default is `NULL`. `bcrat.vars` is only used if `plusBCRAT = TRUE`, 
#' `net = FALSE`, breast cancer is included in the model specification, and 
#' the user does not input `rr.bcrat`. The user only needs to input one of 
#' `rr.bcrat` and `bcrat.vars` to run BCRAPRO+BCRAT; if both are specified, 
#' `rr.bcrat` will be used for the model. 
#' @param rr.bcrat A data frame with 3 columns: the IDs of the probands (`ID`),
#' the probands' BCRAT relative risks for age < 50 (`rr1`), and the probands' 
#' BCRAT relative risks for age >= 50 (`rr2`). When `net = FALSE` and breast 
#' cancer is included in the model, these values will modify the future risk 
#' calculation. The default is `NULL`. `rr.bcrat` is only used if 
#' `plusBCRAT = TRUE`, `net = FALSE`, and breast cancer is included in the 
#' model specification. The user only needs to input one of `rr.bcrat` and 
#' `bcrat.vars` to run BCRAPRO+BCRAT; if both are specified, `rr.bcrat` will be 
#' used for the model. 
#' @param rr.pop A data frame of race-specific for estimates of 1 /
#' (1 - population attributable risk). The default is `rr.ref`,
#' which was estimated from NHIS 2015.
#' 
#' @return A list with two components: 
#' * `posterior.prob`: A list where each element corresponds to a proband in 
#' `proband` and consists of a data frame of posterior carrier probabilities. 
#' There are three columns: `estimate` (the estimate if no ages were imputed, 
#' or the estimate averaged over all age imputations), `lower` (`NA` if no ages 
#' were imputed, or the minimum result for all age imputations), and `upper` 
#' (`NA` if no ages were imputed, or the maximum result for all age 
#' imputations). The rows correspond to the different genotypes. If no genes 
#' were in the model specification, each proband's data frame of carrier 
#' probabilities is replaced by the message 
#' `"No carrier probabilities were requested by the model specification."`
#' If ages were imputed, `ImputeRange` messages giving the lower and upper
#' bounds for the probability of carrying any pathogenic variant will be 
#' printed for each proband. 
#' A `ZeroProb` warning is raised if one or more of the estimated 
#' probabilities is `0` (or `1`), germline testing results were incorporated, 
#' and default sensitivities/specifities of `1` were used. 
#' * `future.risk`: A list where each element corresponds to a proband in 
#' `proband` and consists of a list of data frames of future risk estimates at 
#' fixed intervals for each cancer in the model. Each data frame has three 
#' columns: `estimate` (the estimate if no ages were imputed, or the estimate 
#' averaged over all age imputations), `lower` (`NA` if no ages were imputed, 
#' or the minimum result for all age imputations), and `upper` (`NA` if no 
#' ages were imputed, or the maximum result for all age imputations). Each row 
#' corresponds to the risk estimate for at a different age, e.g. 35, 40, 45, 
#' etc. when `age.by = 5` for a proband whose current age is 30. If no cancers 
#' were in the model specification, each proband's list of future risk 
#' estimates is replaced by the message 
#' `"No future risk estimates were requested by the model specification."`
#' 
#' @details 
#' The `PPP` master function combines pre-processing of model parameters 
#' and the pedigree with calculation of the posterior carrier probabilities and 
#' future risk estimates. First, a model-specific database is built by 
#' subsetting `database` for the genes and cancers implied by `model_spec` or 
#' `genes` and `cancers` (see \code{\link{buildDatabase}}). Then, `pedigree` is 
#' checked for having a valid, consistent structure, and missing ages are 
#' imputed (see \code{\link{checkFam}}). If `pedigree` fails a check, 
#' corrections will be made when possible, and an informative error message 
#' will be issued if not. Finally, the carrier probabilities and future risk 
#' predictions for the user-specified probands are estimated and returned 
#' (see \code{\link{PanelPROCalc}}). 
#' 
#' For convenience, the `PPP` package exports wrapper functions for all 
#' of the pre-specified models: \code{\link{BRCAPRO}}, \code{\link{BRCAPRO5}}, 
#' \code{\link{BRCAPRO6}}, \code{\link{PanelPRO11}}, and 
#' \code{\link{PanelPRO22}}. 
#' 
#' To represent a cancer in `pedigree`, we need to use `short` cancer codes. 
#' See `data.frame(PPP:::CANCER_NAME_MAP)` for available options. For 
#' example, `BC` is the `short` name that maps to `Breast`, so breast cancer 
#' affection status and age of diagnosis should be recorded in columns 
#' named `isAffBC` and `AgeBC`. 
#' 
#' 
#' When the corresponding cancer penetrance information is available in 
#' `PanelPRODatabase`, `PPP` calculates and outputs carrier probabilities 
#' for both homozygous and heterozygous mutations of the genes included in model. 
#' MUTYH is the only gene for which we currently have both homozygous and 
#' heterozygous penetrances; homozygous mutations for all other genes are 
#' considered non-viable (see `PPP:::ALL_GENE_VARIANT_TYPES`).
#'
#' Sets of identical twins/multiple births only contribute to the likelihood 
#' once; see \code{\link{.twinsLikMod}} for more information. Carrier 
#' probabilities returned for individuals in the same set will be identical. 
#' 
#' The R package `cbcrisk`, which is on GitHub but is not on CRAN, will be 
#' automatically installed when the function runs, if is is not installed already. 
#' 
#' @seealso \code{\link{buildDatabase}}, \code{\link{checkFam}}, 
#' \code{\link{PanelPROCalc}}, \code{\link{BRCAPRO}}, \code{\link{BRCAPRO5}}, 
#' \code{\link{BRCAPRO6}}, \code{\link{PanelPRO11}}, \code{\link{PanelPRO22}}
#' @examples
#' PPP(test_fam_1, model_spec = "BRCAPRO", proband = 5, parallel = FALSE)
#' PPP(test_fam_1, genes = c("BRCA1", "BRCA2"), 
#'          cancers = c("Breast", "Ovarian"), proband = 5, parallel = FALSE)
#' 
#' output <- PanelPRO22(test_fam_2, parallel = FALSE, max.mut = 1)
#' # Visualise the risk profile of the proband(s)
#' visRisk(output)
#' @md
#' @export
PPP <- function(pedigree, model_spec = "PanPRO22",
                     proband = pedigree$ID[pedigree$isProband == 1],
                     genes = NULL, cancers = NULL, database = PanelPRODatabase,
                     unknown.race = UNKNOWN_RACE, 
                     unknown.ancestry = UNKNOWN_ANCESTRY,
                     impute.missing.ages = TRUE, 
                     allow.intervention = TRUE, ignore.proband.germ  = FALSE, 
                     max.mut = NULL,  
                     remove.miss.cancers = TRUE, 
                     iterations = 20, max.iter.tries = NULL, 
                     parallel = TRUE, debug = FALSE, net = FALSE, age.by = 5,
                     random.seed = 42L, plusBCRAT = FALSE, bcrat.vars = NULL,
                     rr.bcrat = NULL, rr.pop = rr.ref) {
  
  # check for cbcrisk package, and if it is not present then install it
  if(system.file(package = "cbcrisk") == ""){
    message("The R package cbcrisk v2.0 is required and is being installed.")
    devtools::install_github("ihsajal/CBCRisk@44f7b5ee801a0b4e09977e3f5645e17a38cc7598")
    if(system.file(package = "cbcrisk") == ""){
      stop("Error: cbcrisk package v2.0 installation unsuccessful. Please install cbcrisk v2.0 manually from https://github.com/ihsajal/CBCRisk using commit 44f7b5ee801a0b4e09977e3f5645e17a38cc7598, to use PPP.")
    }
  }
  
  # Extract cancers in pedigree
  ped_cancers = .mapCancerNames(short = .getCancersFromFam(pedigree))
  
  # Attempt to map the model specification to genes and cancers
  if (!is.null(model_spec) && is.null(genes) && is.null(cancers)) {
    model_spec <- match.arg(model_spec, names(MODELPARAMS))
    genes <- MODELPARAMS[[model_spec]]$GENES
    cancers <- MODELPARAMS[[model_spec]]$CANCERS
  } 
  
  # Remove cancers from model spec if they are not in the pedigree
  if (remove.miss.cancers == TRUE) {
    is_missing_cancer = !(cancers %in% ped_cancers)
    if (sum(is_missing_cancer) > 0) {
      rlang::inform(sprintf(
        "%s cancer(s) not in pedigree, so they will be removed from the model specification.", 
        paste(cancers[is_missing_cancer], collapse=", ")), 
        level = "RemoveCancers")
      cancers = cancers[!is_missing_cancer]
    }
  }
  
  # Build model-specific database
  db <- buildDatabase(
    genes = genes, cancers = cancers,
    ppd = database
  )
  
  # Number of relatives in pedigree
  n_rels = nrow(pedigree)
  
  # Check that max.mut makes sense
  # If the user supplied max.mut and it exceeds the number of genes, print a warning
  if (!is.null(max.mut) && max.mut > length(db$MS$GENES)) {
    rlang::inform(sprintf(
      "Maximum allowed mutations cannot exceed total number of genes specified in model; switching from %s to %s", 
      max.mut, length(db$MS$GENES)), level = "CapMaxMut")
    max.mut <- length(db$MS$GENES)
  }
  # If no max.mut was specified, set max.mut = 2 if there are <=30 relatives, 1 otherwise
  if (is.null(max.mut)) {
    if (n_rels <= 30) {
      max.mut = 2
    } else {
      max.mut = 1
    }
    if (max.mut > length(db$MS$GENES)) {
      max.mut <- length(db$MS$GENES)
    } 
    rlang::inform(sprintf("Setting max.mut to %s.", max.mut), level = "SetMaxMut")
  } 
  
  # Check family structure and impute missing ages
  checkl <- checkFam(
    ped = pedigree, db = db,
    proband = proband, 
    unknown.race = unknown.race,
    unknown.ancestry = unknown.ancestry,
    ignore.proband.germ = ignore.proband.germ, 
    impute.missing.ages = impute.missing.ages,
    impute.times = iterations, 
    max.iter.tries = max.iter.tries, 
    random.seed = random.seed
  )
  
  # If breast cancer is in the model, the user specifies BRCAPRO+BCRAT 
  # arguments, and if the user asks for crude risk, set up BRCAPRO+BCRAT
  if (plusBCRAT == TRUE) {
    if (!is.null(rr.bcrat) | !is.null(bcrat.vars)) {
      if ("Breast" %in% db$MS$CANCERS & net == FALSE) {
        if (is.null(rr.bcrat) & !is.null(bcrat.vars)) {
          rr.bcrat <- calc.rr.bcrat(bcrat.vars)
        }
        if (any(is.na(rr.bcrat)) | ncol(rr.bcrat) != 3){
          rlang::abort("Neither bcrat.vars nor rr.bcrat is in a valid format.")
        }
      } else if (!("Breast" %in% db$MS$CANCERS)) {
        rlang::inform("BRCAPRO+BCRAT arguments not used because breast cancer is not in the model.")
      } else {
        rlang::inform("BRCAPRO+BCRAT arguments not used because net risk is being reported.")
      }
    }
  }
  
  # Register parallel cores if ages should be imputed in parallel
  if (parallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 2)
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
  doRNG::registerDoRNG()
  
  # Calculate posterior genotype distribution and future risk
  res <- Map(function(ped, lm) {
    PanelPROCalc(ped, lm,
                 proband = checkl$proband, db = db,
                 ifmod = allow.intervention, net = net,
                 age.by = age.by, max.mut = max.mut,
                 debug = debug,
                 rr.bcrat = rr.bcrat, rr.pop = rr.pop
    )
  }, ped = checkl$ped_list, lm = checkl$lms
  )
  
  if (parallel) on.exit(parallel::stopCluster(cl))
  
  # Collect results among potentially separated families
  if (length(res) == 1) {
    return(res[[1]])
  } else {
    pprob <- do.call(rbind, lapply(1:length(res), 
                                   function(i) res[[i]]$posterior.prob))
    frisk <- unlist(lapply(1:length(res), 
                           function(i) res[[i]]$future.risk), recursive = FALSE)
    return(list(posterior.prob = pprob, future.risk = frisk))
  }
}


#' BRCAPRO
#'
#' Runs PPP with `genes = c("BRCA1", "BRCA2")` and 
#' `cancers = c("Breast", "Ovarian", "Contralateral")`, as in 
#' `PPP:::MODELPARAMS$BRCAPRO`. 
#' @param ... Additional arguments to be passed into `PPP`.
#' @examples 
#' data(test_fam_1)
#' BRCAPRO(test_fam_1, parallel = FALSE)
#' @family predefined
#' @seealso \code{\link{PPP}}
#' @export
BRCAPRO <- function(...) {
  PPP(model_spec = "BRCAPRO", genes = NULL, cancers = NULL, ...)
}


#' MMRPRO
#'
#' Runs PPP with `genes = c("MLH1", "MSH2", "MSH6")` and
#' `cancers = c("Colorectal", "Endometrial")`, as in
#' `PPP:::MODELPARAMS$MMRPRO`. 
#' @param ... Additional arguments to be passed into `PPP`.
#' @examples 
#' data(test_fam_3)
#' MMRPRO(test_fam_3, parallel = FALSE)
#' @family predefined
#' @seealso \code{\link{PPP}}
#' @export
MMRPRO <- function(...) {
  PPP(model_spec = "MMRPRO", genes = NULL, cancers = NULL, ...)
}


#' BRCAPRO-5
#'
#' Runs PPP with `genes = c("BRCA1", "BRCA2", "ATM", "CHEK2", "PALB2")` and
#' `cancers = c("Breast", "Ovarian", "Contralateral")`, as in
#' `PPP:::MODELPARAMS$BRCAPRO5`.
#' @param ... Additional arguments to be passed into `PPP`.
#' @examples 
#' data(test_fam_1)
#' BRCAPRO5(test_fam_1, parallel = FALSE)
#' @family predefined
#' @seealso \code{\link{PPP}}
#' @export
BRCAPRO5 <- function(...) {
  PPP(model_spec = "BRCAPRO5", genes = NULL, cancers = NULL, ...)
}


#' BRCAPRO-6
#'
#' Runs PPP with 
#' `genes = c("BRCA1", "BRCA2", "MLH1", "MSH2", "MSH6", "CDKN2A")`
#' and `cancers = c("Breast", "Ovarian", "Contralateral")`, as in 
#' `PPP:::MODELPARAMS$BRCAPRO6`.
#' @param ... Additional arguments to be passed into `PPP`.
#' @examples
#' data(test_fam_1)
#' BRCAPRO6(test_fam_1, parallel = FALSE)
#' @family predefined
#' @seealso \code{\link{PPP}}
#' @export
BRCAPRO6 <- function(...) {
  PPP(model_spec = "BRCAPRO6", genes = NULL, cancers = NULL, ...)
}


#' PPP-11
#'
#' Runs PPP with 11 genes 
#' (`genes = c("BRCA1", "BRCA2", "ATM", "PALB2", "CHEK2", "EPCAM", "PMS2", "MLH1", "MSH2", "MSH6", "CDKN2A")`)
#' and 11 cancers 
#' (`cancers = c("Brain", "Breast", "Colorectal", "Contralateral", "Endometrial", "Gastric", "Kidney", "Melanoma", "Ovarian", "Pancreas", "Prostate", "Small intestine")`),
#' as in `PPP:::MODELPARAMS$PanPRO11`.
#' @param ... Additional arguments to be passed into `PPP`.
#' @examples 
#' data(test_fam_3)
#' PanelPRO11(test_fam_3, parallel = FALSE)
#' @family predefined
#' @seealso \code{\link{PPP}}
#' @export
PanelPRO11 <- function(...) {
  PPP(model_spec = "PanPRO11", genes = NULL, cancers = NULL, ...)
}


#' PPP-22
#'
#' Runs PPP with 22 genes 
#' (`genes = c("ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CDK4", "CDKN2A", "CHEK2", "EPCAM", "MLH1", "MSH2", "MSH6", "MUTYH", "NBN", "PALB2", "PMS2", "PTEN", "RAD51C", "RAD51D", "STK11", "TP53")`)
#' and 17 cancers 
#' (`cancers = c("Brain", "Breast", "Colorectal", "Contralateral", "Endometrial", "Gastric", "Kidney", "Leukemia", "Melanoma", "Ovarian", "Osteosarcoma", "Pancreas", "Prostate", "Small Intestine", "Soft Tissue Sarcoma", "Thyroid", "Urinary Bladder", "Hepatobiliary")`),
#' as in `PPP:::MODELPARAMS$PanPRO22`.
#' @param ... Additional arguments to be passed into `PPP`.
#' @examples 
#' data(test_fam_2)
#' PanelPRO22(test_fam_2, parallel = FALSE, max.mut = 1)
#' @family predefined
#' @seealso \code{\link{PPP}}
#' @export
PanelPRO22 <- function(...) {
  PPP(model_spec = "PanPRO22", genes = NULL, cancers = NULL, ...)
}
