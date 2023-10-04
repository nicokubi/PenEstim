#' Modify individual cancer penetrances using a hazard ratio 
#' 
#' @param hr Hazard ratio to be used to modify cancer penetrances `rawp`. 
#' @param rawp A numeric vector of age-specific cancer penetrances for a given 
#' cancer, genotype, race, and sex. 
#' @param intage A numeric value indicating the age at which the individual 
#' received the intervention that modifies risk. 
#' @return A numeric vector of modified cancer penetrances for an individual. 
#' @family modifypen
postIntervHR <- function(hr, rawp, intage) {
  asur <- c(1, 1 - cumsum(rawp))
  hra <- c(rep(1, times = intage - 1), rep(hr, times = MAXAGE - intage + 1))
  haza <- hra * rawp / asur[1:MAXAGE]
  asurn <- c(1, cumprod(1 - haza))
  return(asurn[1:MAXAGE] - asurn[-1])
}


#' Modify individual cancer penetrances using a relative risk
#'
#' @param rr Relative risk to be used to modify cancer penetrances `rawp`. 
#' @param rawp A numeric vector of age-specific cancer penetrances for a given 
#' cancer, genotype, race, and sex. 
#' @param intage A numeric value indicating the age at which the individual 
#' received the intervention that modifies risk. 
#' @return A numeric vector of modified cancer penetrances for an individual. 
#' @family modifypen
postIntervRR <- function(rr, rawp, intage) {
  rr <- c(rep(1, times = intage - 1), rep(rr, times = MAXAGE - intage + 1))
  return(rawp * rr)
}


#' Modify cancer penetrances based on hazard ratios or relative risks for interventions
#'
#' @param raw_famp A numeric array of cancer penetrances for individuals in 
#' `fam_distinct`, created inside of the `calcCancerPenetrance` function. 
#' @param fam_distinct A checked pedigree data frame returned by 
#' \code{\link{checkFam}}. 
#' @param riskmod Risk modifier information contained in the `riskmod` 
#' component in the database returned by \code{\link{buildDatabase}}. 
#' @details Risk modifiers are assumed to be independent. 
#' @return A numeric array of modified cancer penetrances with the same 
#' dimensions as `raw_famp`. 
#' @family modifypen
riskmod_perfam <- function(raw_famp, fam_distinct, riskmod) {
  # List of interventions for each member of pedigree
  interv_list <- fam_distinct$riskmod 
  # List of intervention ages for each member of pedigree
  age_list <- fam_distinct$interAge 
  
  # Cancers extracted from pedigree
  all_cancers <- .getCancersFromFam(fam_distinct)
  # Get cancer affection ages for each member of pedigree
  aff_cancers <- lapply(1:nrow(fam_distinct), function(i) {
    vec <- unlist(fam_distinct[i, paste0("Age", all_cancers)])
    names(vec) <- all_cancers
    vec[vec > 1 & vec < 94]
  })

  # Helper function for individual risk modification
  riskmod_ind <- function(rlist, alist, aff_cancer_vec, 
                          rawp, riskmod, gender) {
    # If individual had no interventions, do not modify penetrance
    if ("None" %in% rlist) {
      return(rawp)
    } else { # Otherwise
      # Ensure that there is the same number of interventions and intervention 
      # ages
      if (length(rlist) != length(alist)) {
        rlang::abort("The number of interventions does not match the number of 
                     intervention ages.")
      }
      # Convert lists of interventions/intervention ages into vectors
      alist <- as.numeric(alist)[order(as.numeric(alist))]
      rlist <- rlist[order(alist)] # order interventions by age

      # Subset risk modifier information for each individual
      riskmod_list <- Map(function(interv, g, age) {
        subset_array(riskmod, 
                     c("Intervention", "Sex", "IntervAge", "DataType"), 
                     list(interv, g, age, c("HR", "RR")), drop = FALSE)
      }, interv = rlist, g = gender, age = alist)

      riskmod_list <- lapply(riskmod_list, function(arr) {
        abind::adrop(arr, drop = c(3, 4, 5))
      })

      # Check that riskmod_list has the same structure for cancer and genes
      id <- which(lengths(riskmod_list) > 0)[1] 
      if (!all(dimnames(rawp)$Cancer == dimnames(riskmod_list[[id]])$Cancer)) {
        rlang::abort("Risk modifications structure does not match cancer structure.")
      }
      if (!all(dimnames(rawp)$Gene_type == dimnames(riskmod_list[[id]])$Gene_type)) {
        rlang::abort("Risk modifications structure does not match cancer structure.")
      }

      # Iterate through individual's interventions
      # Assume that none of the risk modifiers have both HR and RR information
      for (i in seq_along(rlist)) {
        # Intervention age
        mod_age <- alist[i]
        # Check that intervention age is not missing
        if (is.na(mod_age) || mod_age == -999) {
          rlang::abort("Risk modification age is missing.")
        }
        
        # Risk will only modify cancers that were diagnosed after the intervention
        aff_cancer_vec <- abs(aff_cancer_vec) # turn -999 to 999
        cancers_affected <- .mapCancerNames(
          short = names(aff_cancer_vec)[aff_cancer_vec > mod_age]
        )

        # Iterate through affected cancers
        for (cancer in cancers_affected) {
          genetypes <- dimnames(riskmod_list[[id]])$Gene
          for (genetype in genetypes) { # Iterate through genes
            # Extract HR and RR information from database input
            hr <- riskmod_list[[i]][cancer, genetype, "HR"]
            rr <- riskmod_list[[i]][cancer, genetype, "RR"]
            if ((hr < 1) & (rr == 1)) { 
              # HR is stored in database, RR is not
              rawp[i, cancer, genetype, ] <- 
                postIntervHR(hr, rawp = rawp[i, cancer, genetype, ], 
                             intage = alist[i])
            } else if ((hr == 1) & (rr < 1)) { 
              # RR is stored in database, HR is not
              rawp[i, cancer, genetype, ] <- 
                postIntervRR(rr, rawp = rawp[i, cancer, genetype, ], 
                             intage = alist[i])
            }
          }
        }
      }
      return(rawp)
    }
  }
  
  # Initialize container for storing modified penetrances
  mod_famp <- raw_famp 

  # Modify the risk for each individual 
  for (i in seq(nrow(fam_distinct))) {
    mod_famp[i, , , ] <- riskmod_ind(
      rlist = interv_list[[i]],
      alist = age_list[[i]],
      aff_cancer_vec = aff_cancers[[i]],
      gender = fam_distinct$Sex[i],
      rawp = raw_famp[i, , , , drop = FALSE],
      riskmod = riskmod
    )
  }
  return(mod_famp)
}
