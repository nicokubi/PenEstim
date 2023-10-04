#' Get current age-bounding information for an individual
#'
#' @param ped A pedigree data frame returned by 
#' \code{\link{.pedStructureCheck}}. 
#' @param ind An individual ID. 
#' @param rel_l The list of relationships in `ped` returned by 
#' \code{\link{.pedStructureCheck}}. 
#' @return A list with current age-bounding information for the individual with 
#' ID `ind`, based on their relatives' ages in `ped. The 9 list components are: 
#' * `ID`: The individual's ID (the same as `ind`). 
#' * `curage`: The individual's current age, if non-missing; `NA` otherwise. 
#' * `max_aff`: The individual's maximum age of cancer affection. 
#' * `sex`: The individual's sex. 
#' * `interv`: The individual's maximum age at which they received an 
#' intervention. 
#' * `mage`: The individual's mother's current age, if she is alive and her 
#' current age is not missing; `NA` otherwise. 
#' * `fage`: The individual's father's current age, if he is alive and his 
#' current age is not missing; `NA` otherwise. 
#' * `cage`: The current age of the individual's oldest living child with a 
#' known age; `NA` if no living children have non-missing ages. 
#' * `sage`: The current ages of the individual's living siblings; `NA` if no 
#' living siblings have non-missing ages. 
#' @family impute
.getAgeFromPed <- function(ped, ind, rel_l) {
  if (!ind %in% ped$ID) {
    rlang::abort(sprintf("ID %s must be presented in pedigree column ID", ind))
  }
  # Extract row index corresponding to individual ID
  i <- which(ped$ID == ind) 
  
  # Initialize output list for storing relatives' ages
  rel <- list()
  
  # Individual's ID
  rel$ID <- ind
  # Individual's current age, if non-missing
  rel$curage <- ifelse(ped[i, "CurAge"] == -999, NA, ped[i, "CurAge"])

  # Individual's maximum age of affection
  cancers <- .getCancersFromFam(ped)
  if (length(cancers) > 0) {
    rel$max_aff <- max(
      ped[i, paste0("Age", cancers)],
      ped[i, paste0("Age", cancers, "_upper")]
    )
  } else { 
    rel$max_aff = NA
  }

  # Individual's sex 
  rel$sex <- ped[i, "Sex"]
  
  # Individual's maximum age of intervention
  rel$interv <- .safeGet(max, ped$interAge[[i]][ped$riskmod[[i]] != "None"])
  
  # Mother's current age, if alive and current age is not missing
  mother_id <- rel_l[[i]]$mid
  
  if (length(mother_id) == 0) {
    rel$mage <- NA
  } else {
    if (mother_id == -999) {
      rel$mage <- NA
    }
    else {
      mage <- ped[ped$ID == mother_id, "CurAge"]
      if (ped[ped$ID == mother_id, "isDead"] == 1) {
        rel$mage <- NA
      }
      else {
        rel$mage <- ifelse(mage == -999, NA, mage)
      }
    }
  }
  
  # Father's current age, if alive and current age is not missing
  father_id <- rel_l[[i]]$fid
  if (length(father_id) == 0) {
    rel$fage <- NA
  } else {
    if (father_id == -999) {
      rel$fage <- NA
    }
    else {
      fage <- ped[ped$ID == father_id, "CurAge"]
      if (ped[ped$ID == father_id, "isDead"] == 1) {
        rel$fage <- NA
      }
      else {
        rel$fage <- ifelse(fage == -999, NA, fage)
      }
    }
  }

  # Current age of oldest living child with non-missing age
  children_ages <- ped[ped$ID %in% rel_l[[i]]$cid, "CurAge"]
  children_dead <- ped[ped$ID %in% rel_l[[i]]$cid, "isDead"]
  children_ages <- children_ages[children_ages != -999 & children_dead != 1]
  rel$cage <- .safeGet(max, children_ages)

  # Current ages of living siblings with non-missing ages
  sib_ages <- ped[ped$ID %in% rel_l[[i]]$fsibid, "CurAge"]
  sib_dead <- ped[ped$ID %in% rel_l[[i]]$fsibid, "isDead"]
  rel$sage <- sib_ages[sib_ages != -999 & sib_dead != 1]
  if (length(rel$sage) == 0) rel$sage <- NA
  
  return(rel)
}


#' Get lower and upper current age bounds for an individual
#'
#' @details 
#' When the current age for the individual is non-missing, both the lower and 
#' upper bounds are set to the known current age. 
#' 
#' Otherwise, the lower current age bound is the maximum of: 
#' * their maximum age of cancer affection
#' * their maximum intervention age
#' * the current age of their oldest child + `min.fertile.age.female` (if this 
#' individual is female) or `min.fertile.age.male` (if this individual is male)
#' * the current age of their youngest sibling - `max.sib.diff`
#' 
#' If this process does not yield a non-missing lower bound greater than `1`, 
#' it will be set to `1`. 
#' 
#' The upper current age bound is the minimum of
#' * their mother's current age - `min.fertile.age.female`
#' * their father's current age - `min.fertile.age.male`
#' * the current age of their oldest sibling + `max.sib.diff`
#' 
#' `min.fertile.age.female` (minimum female fertility age), 
#' `min.fertile.age.male`(minimum male fertility age),  and `max.sib.diff` 
#' (maximum  sibling age difference) are calculated based on the distributions 
#' in `AgeDistribution`. 
#' 
#' If this process does not yield a non-missing upper bound less than 
#' `PPP:::MAXAGE`, it will be set to `PPP:::MAXAGE`. 
#'
#' @param rel_ages A list of current age-bounding information for an individual 
#' returned by \code{\link{.getAgeFromPed}}. 
#' @return A vector of length two with the lower and upper current age bounds 
#' for the individual with ID `rel_ages$ID`. 
#' @family impute
#' @md
.boundAgeFromRelatives <- function(rel_ages) {
  
  # Calculate minimum birth ages and maximum difference in sibling ages based 
  # on data in AgeDistribution
  min.fertile.age.female = min(AgeDistribution$female_birth)
  min.fertile.age.male = min(AgeDistribution$male_birth)
  max.sib.diff = max(AgeDistribution$sib_diff)
  
  # The individual's ID
  id <- rel_ages$ID

  # If the individual's current age is known, return it as both the lower and 
  # upper bounds for current age
  if (!is.na(rel_ages$curage)) {
    return(rep(rel_ages$curage, 2))
  }

  # When the individual's current age is not known, try to find bounds
  # Possible lower bounds for current age
  possible_lower_bounds <- c(
    # Individual is at least their maximum age of cancer affection
    rel_ages$max_aff,
    # Individual is at least their maximum intervention age
    rel_ages$interv,
    # Individual is at least min.fertile.age.female (if female) or 
    # min.fertile.age.male (if male) years older than their oldest child
    rel_ages$cage + ifelse(rel_ages$sex == "Male", 
                           min.fertile.age.male, min.fertile.age.female),
    # Individual is at least max.sib.diff years younger than their youngest sibling
    min(rel_ages$sage) - max.sib.diff
  )
  # Possible upper bounds for current age
  possible_upper_bounds <- c(
    # Individual is at most min.fertile.age.female years younger than their mother
    rel_ages$mage - min.fertile.age.female,
    # Individual is at most min.fertile.age.male years younger than their father
    rel_ages$fage - min.fertile.age.male,
    # Individual is at most max.sib.diff years older than their oldest sibling
    max(rel_ages$sage) + max.sib.diff
  )
  # Convert negative (missing) values to NA
  possible_lower_bounds[possible_lower_bounds <= 0] <- NA
  possible_upper_bounds[possible_upper_bounds <= 0] <- NA

  # Choose the smallest of possible_upper_bounds as the upper bound
  ub <- .safeGet(min, possible_upper_bounds)
  # Choose the largest of possible_lower_bounds as the lower bound
  lb <- .safeGet(max, possible_lower_bounds)

  # If one of the bounds is still unknown or invalid, replace with extreme values
  if (is.na(ub) | ub < 0 | ub > MAXAGE | is.infinite(ub)) {
    ub <- MAXAGE
  }
  if (is.na(lb) | lb < 1 | is.infinite(lb)) {
    lb <- 1
  }

  # Ensure that the choice of bounds is valid
  if (ub < lb) {
    rlang::inform(sprintf(
      "ID %s cannot get valid capping estimate, with upper bound %s and lower bound %s",
      id, ub, lb
    ))
    
    # If upper bound is lower than lower bound, swap them
    temp = lb
    lb = ub
    ub - lb
    
    # If one of the bounds is still invalid, replace with extreme values
    if (ub < 0 | ub > MAXAGE) {
      ub <- MAXAGE
    }
    if (lb < 1) {
      lb <- 1
    }
  }
  
  return(c(lb, ub))
}

#' Impute cancer age for an individual
#'
#' Samples age of cancer affection using cancer penetrances as sample 
#' weights, subject to lower and upper age bounds for the individual. 
#'
#' @param max_age_allowed A numeric value; the upper bound to be used for 
#' imputing cancer age.  
#' @param penet A numeric vector of cancer penetrances for the individual with 
#' length equal to `PPP:::MAXAGE`. 
#' @param least_possible_age A numeric value; the lower bound to be used for 
#' imputing cancer age. The default is `1`. 
#' @return An integer value representing the imputed cancer affection age. If 
#' no cancer age can be sampled (i.e. all probabilities are 0), NA will be 
#' returned. 
#' @family impute
.cancerAgeImpute <- function(max_age_allowed, penet, least_possible_age = 1) {
  # If max_age_allowed is NA, use PPP:::MAXAGE 
  if (is.na(max_age_allowed)) max_age_allowed <- length(penet)
  
  # Subset the penetrance vector to only include ages in the bounded range
  capped_vec <- penet[least_possible_age:max_age_allowed]
  if (all(capped_vec == 0)) {
    # If the subsetted vector is all zeros, return NA
    return(NA)
  }
  
  # Renormalize penetrances
  capped_vec <- capped_vec / sum(capped_vec) 
  
  # Check that subsetted penetrance vector is valid
  if (any(capped_vec < 0)) {
    rlang::abort("Sampling portion should not smaller than 0!")
  }
  
  # Impute cancer age, use penetrances as sampling weights
  sample(least_possible_age:max_age_allowed, size = 1, 
         replace = TRUE, prob = capped_vec)
}


#' Impute current age for an individual
#'
#' @details 
#' Samples current age based on the empirical age distributions in `dist_age`, 
#' subject to the age bounds in `ab` and age information in `rel_ages`. 
#' 
#' First, attempt to impute the following ages when the age information is 
#' known: 
#' 
#' 1. If the individual's mother's current age is known, calculate the 
#' difference between the mother's age and an age randomly sampled from the 
#' distribution of female birth ages
#' 2. If the individual's father's current age is known, calculate the 
#' difference between the father's age and an age randomly sampled from the 
#' distribution of male birth ages
#' 3. If at least one of the individual's children's current ages is 
#' known, calculate the sum of the oldest child's current age and an age 
#' randomly sampled from the distribution of female ages at first birth, if the 
#' individual is female; calculate the sum of the oldest child's current age 
#' and an age randomly sampled from the distribution of male ages at first 
#' birth, if the individual is male 
#' 4. If at least one of the individual's siblings' current ages is known, 
#' calculate the sum of the mean sibling age and +/- an age difference randomly 
#' sampled from the distribution of sibling age differences, where the +/- sign 
#' is determined by a fair coin toss 
#' 
#' The final imputed age returned by the function is the mean of the ages 
#' successfully imputed by the above approaches, bounded by the lower and upper 
#' age bounds in `ab`. If no ages were successfully imputed by any of the 
#' approaches, sample from a uniform distribution bounded by the lower and 
#' upper age bounds in `ab`. 
#'
#' @param rel_ages A list of the individual's current age-bounding information 
#' returned by \code{\link{.getAgeFromPed}}. 
#' @param ab A numeric vector of length two with the individual's current age 
#' bounds. 
#' @param dist_age A list containing empirical distributions of female and male 
#' ages at first birth; female and male birth ages; and sibling age 
#' differences. The default is \code{\link{AgeDistribution}}. 
#' @return An integer value representing the imputed current age. 
#' @family impute
#' @md
.currentAgeImpute <- function(rel_ages, ab, 
                              dist_age = AgeDistribution) {
  
  # Extract individual's lower and upper current age bounds
  bmin <- min(ab)
  bmax <- max(ab)
  
  # Initialize vector for storing imputed ages from different approaches
  imputed_ages <- rep(NA, 4)
  names(imputed_ages) = c("moth", "fath", "child", "sib")
  
  # If mother's current age is known, use her age - age randomly sampled from 
  # the distribution of female birth ages
  if (!is.na(rel_ages$mage)) { 
    imputed_ages["moth"] <- rel_ages$mage - 
      sample(dist_age$female_birth, size = 1)
  } 
  
  # If father's current age is known, use his age - age randomly sampled from 
  # the distribution of male birth ages
  if (!is.na(rel_ages$fage)) { 
    imputed_ages["fath"] <- rel_ages$fage - 
      sample(dist_age$male_birth, size = 1)
  } 
  
  # If the oldest child's current age is known, use their age + age randomly 
  # sampled from the distribution of (male or female) ages at first birth
  if (!is.na(rel_ages$cage)) { 
    if (rel_ages$sex == "Male") {
      first_birth_dist <- dist_age$male_first_birth
    } else {
      first_birth_dist <- dist_age$female_first_birth
    }
    imputed_ages["child"] <- rel_ages$cage + sample(first_birth_dist, size = 1)
  }
  
  # If at least one of the siblings' current ages is known, use the mean of 
  # their ages +/- age difference randomly sampled from distribution of 
  # sibling age differences
  if (!any(is.na(rel_ages$sage))) {
    imputed_ages["sib"] <- mean(rel_ages$sage) + sample(c(-1,1), 1, 0.5) *
      sample(dist_age$sib_diff, size = 1)
  } 
  
  # Set imputed ages above upper bound to upper bound
  imputed_ages <- ifelse(imputed_ages > bmax, bmax, imputed_ages)
  # Set imputed ages lower bound to lower bound
  imputed_ages <- ifelse(imputed_ages < bmin, bmin, imputed_ages)
  
  # Calculate average of available imputed ages
  final_age <- ceiling(mean(imputed_ages, na.rm = TRUE))
  
  # If no ages could be imputed using any of the approaches, sample from a 
  # uniform distribution bounded by the lower and upper age bounds. 
  if (is.na(final_age)){
    final_age <- sample(bmin:bmax, size = 1)
  }
  
  return(final_age)
}


#' Find the indices of the first-degree relatives of a given individual in a pedigree
#'
#' @param pedigree A pedigree data frame that includes at least the `ID`, 
#' `MotherID`, and `FatherID` columns. 
#' @param ii The index (not ID) of the individual of interest in the pedigree. 
#' @return A data frame with two columns: `index` of first-degree relatives and 
#' their `lineage`. 
#' @details Adapted from the function of the same name in `BayesMendel` 
#' package, with little modification. 
#' @family impute
.firstDegreeRelative <- function(pedigree, ii) {
  numSample <- dim(pedigree)[1]
  rlt <- NULL
  FatherID <- pedigree$FatherID[ii]
  MotherID <- pedigree$MotherID[ii]
  ID <- pedigree$ID[ii]
  for (i in 1:numSample) {
    if (i == ii) {
      next
    }
    if ((pedigree$FatherID[i] == FatherID && 
         pedigree$MotherID[i] == MotherID) &&
      (FatherID != -999 && MotherID != -999)) {
      rltTmp <- data.frame(matrix(data = NA, nrow = 1, ncol = 2))
      names(rltTmp) <- c("index", "lineage")
      rltTmp[1] <- i
      rltTmp[2] <- "S"
      rlt <- rbind(rlt, rltTmp)
      next
    }
    if (pedigree$ID[i] == FatherID && FatherID != -999) {
      rltTmp <- data.frame(matrix(data = NA, nrow = 1, ncol = 2))
      names(rltTmp) <- c("index", "lineage")
      rltTmp[1] <- i
      rltTmp[2] <- "F"
      rlt <- rbind(rlt, rltTmp)
      next
    }
    if (pedigree$ID[i] == MotherID && MotherID != -999) {
      rltTmp <- data.frame(matrix(data = NA, nrow = 1, ncol = 2))
      names(rltTmp) <- c("index", "lineage")
      rltTmp[1] <- i
      rltTmp[2] <- "M"
      rlt <- rbind(rlt, rltTmp)
      next
    }
    if ((pedigree$FatherID[i] == ID || pedigree$MotherID[i] == ID) &&
      (ID != -999)) {
      rltTmp <- data.frame(matrix(data = NA, nrow = 1, ncol = 2))
      names(rltTmp) <- c("index", "lineage")
      rltTmp[1] <- i
      rltTmp[2] <- "C"
      rlt <- rbind(rlt, rltTmp)
      next
    }
  }
  rlt
}


#' Find the indices of the second-degree relatives of a given individual in a pedigree
#'
#' @param pedigree A pedigree data frame that includes at least the `ID`, 
#' `MotherID`, and `FatherID` columns. 
#' @param ii The index (not ID) of the individual of interest in the pedigree. 
#' @return A data frame with two columns: `index` of second-degree relatives 
#' and their `lineage`. 
#' @details Adapted from the function of the same name in `BayesMendel` 
#' package, with little modification. Calls `.firstDegreeRelative`, so it is 
#' inefficient if both \code{\link{.firstDegreeRelative}} and 
#' \code{\link{.secondDegreeRelative}} are called for the same pedigree. 
#' @family impute
.secondDegreeRelative <- function(pedigree, ii) {
  numSample <- dim(pedigree)[1]
  fdRelative <- .firstDegreeRelative(pedigree, ii)
  rlt <- NULL
  if (is.null(fdRelative)) {
    return(rlt)
  }
  numFDR <- dim(fdRelative)[1]
  for (i in 1:numFDR) {
    if (fdRelative[i, 2] == "F") {
      FatherID <- pedigree$FatherID[fdRelative[i, 1]]
      MotherID <- pedigree$MotherID[fdRelative[i, 1]]
      ID <- pedigree$ID[fdRelative[i, 1]]
      for (j in 1:numSample) {
        if (j == fdRelative[i, 1] || j == ii) {
          next
        }
        if ((pedigree$ID[j] == FatherID || pedigree$ID[j] == MotherID) &&
          (pedigree$ID[j] != -999)) {
          rltTmp <- data.frame(matrix(
            data = NA, nrow = 1,
            ncol = 2
          ))
          names(rltTmp) <- c("index", "lineage")
          rltTmp[1] <- j
          rltTmp[2] <- "F"
          rlt <- rbind(rlt, rltTmp)
          next
        }
        if ((pedigree$FatherID[j] == FatherID && 
             pedigree$MotherID[j] == MotherID) &&
          (FatherID != -999 && MotherID != -999)) {
          rltTmp <- data.frame(matrix(
            data = NA, nrow = 1,
            ncol = 2
          ))
          names(rltTmp) <- c("index", "lineage")
          rltTmp[1] <- j
          rltTmp[2] <- "F"
          rlt <- rbind(rlt, rltTmp)
          next
        }
        if (pedigree$FatherID[j] == ID && ID != -999) {
          if (!(j %in% fdRelative[, 1])) {
            rltTmp <- data.frame(matrix(
              data = NA, nrow = 1,
              ncol = 2
            ))
            names(rltTmp) <- c("index", "lineage")
            rltTmp[1] <- j
            rltTmp[2] <- "F"
            rlt <- rbind(rlt, rltTmp)
            next
          }
        }
      }
      next
    }
    if (fdRelative[i, 2] == "M") {
      FatherID <- pedigree$FatherID[fdRelative[i, 1]]
      MotherID <- pedigree$MotherID[fdRelative[i, 1]]
      ID <- pedigree$ID[fdRelative[i, 1]]
      for (j in 1:numSample) {
        if (j == fdRelative[i, 1] || j == ii) {
          next
        }
        if ((pedigree$ID[j] == FatherID || pedigree$ID[j] == MotherID) &&
          (pedigree$ID[j] != -999)) {
          rltTmp <- data.frame(matrix(
            data = NA, nrow = 1,
            ncol = 2
          ))
          names(rltTmp) <- c("index", "lineage")
          rltTmp[1] <- j
          rltTmp[2] <- "M"
          rlt <- rbind(rlt, rltTmp)
          next
        }
        if ((pedigree$FatherID[j] == FatherID && 
             pedigree$MotherID[j] == MotherID) &&
          (FatherID != -999 && MotherID != -999)) {
          rltTmp <- data.frame(matrix(
            data = NA, nrow = 1,
            ncol = 2
          ))
          names(rltTmp) <- c("index", "lineage")
          rltTmp[1] <- j
          rltTmp[2] <- "M"
          rlt <- rbind(rlt, rltTmp)
          next
        }
        if (pedigree$MotherID[j] == ID && ID != -999) {
          if (!(j %in% fdRelative[, 1])) {
            rltTmp <- data.frame(matrix(
              data = NA, nrow = 1,
              ncol = 2
            ))
            names(rltTmp) <- c("index", "lineage")
            rltTmp[1] <- j
            rltTmp[2] <- "M"
            rlt <- rbind(rlt, rltTmp)
            next
          }
        }
      }
      next
    }
    if (fdRelative[i, 2] == "S") {
      ID <- pedigree$ID[fdRelative[i, 1]]
      for (j in 1:numSample) {
        if (j == fdRelative[i, 1] || j == ii) {
          next
        }
        if ((pedigree$FatherID[j] == ID || pedigree$MotherID[j] == ID) &&
          (ID != -999)) {
          rltTmp <- data.frame(matrix(
            data = NA, nrow = 1,
            ncol = 2
          ))
          names(rltTmp) <- c("index", "lineage")
          rltTmp[1] <- j
          rltTmp[2] <- "S"
          rlt <- rbind(rlt, rltTmp)
          next
        }
      }
      next
    }
    if (fdRelative[i, 2] == "C") {
      ID <- pedigree$ID[fdRelative[i, 1]]
      for (j in 1:numSample) {
        if (j == fdRelative[i, 1] || j == ii) {
          next
        }
        if ((pedigree$FatherID[j] == ID || pedigree$MotherID[j] == ID) &&
          (ID != -999)) {
          rltTmp <- data.frame(matrix(
            data = NA, nrow = 1,
            ncol = 2
          ))
          names(rltTmp) <- c("index", "lineage")
          rltTmp[1] <- j
          rltTmp[2] <- "C"
          rlt <- rbind(rlt, rltTmp)
          next
        }
      }
      next
    }
  }
  rlt
}


#' Impute missing ages in pedigree
#'
#' Imputes missing current, cancer affection, and death ages for relatives in 
#' `ped`, subject to age-bounding information. 
#'
#' @details 
#' First, we impute the missing censoring ages (stored as `ped$CurAge`) of 
#' living first-degree relatives, followed by second-degree relatives and then 
#' more distant relatives. If more than one proband is supplied, the first one 
#' will be used to determine first- and second-degree relatives. 
#' 
#' Within each group of relatives, we estimate current age bounds for each 
#' individual (see \code{\link{.getAgeFromPed}} documentation). The relative 
#' with the narrowest age bound interval in the relative group has their 
#' current age imputed first, from a normal distribution truncated by their age 
#' bounds (see \code{\link{.currentAgeImpute}}  documentation). The pedigree is 
#' updated with this relative's imputed age and the age bounds are then 
#' re-calculated for the remaining relatives in the relative group. This 
#' process repeats until no living individuals with missing censoring ages are 
#' left in the relative group.
#'
#' After all missing censoring ages for living relatives have been imputed, we 
#' impute missing cancer affection ages. For each missing age, the upper bound 
#' is the minimum of the individual's censoring age, their upper bound for this 
#' cancer's affection age, and `PPP:::MAXAGE`; the lower bound is the 
#' maximum of their lower bound for this cancer's affection age, their ages of 
#' preventative interventions (if relevant), and `1`. The affection age is 
#' then sampled based on the crude SEER cancer penetrance estimates (subset by 
#' the individual's race and sex) within the age bounds (see 
#' \code{\link{.cancerAgeImpute}} documentation). If a missing cancer age 
#' cannot be sampled such that it satisfies its age bounds, this imputation run 
#' will be discarded and re-attempted, with a warning message issued. If the 
#' desired number of runs (`iterations`, default `20`) cannot be imputed within 
#' `max.iter.tries` (default `5*iterations`) attempts, an error will be raised. 
#'
#' Finally, missing death ages (stored in `CurAge` in the pedigree) of 
#' individuals who are affected for at least one cancer and/or reported at 
#' least one intervention are imputed, using the maximum of their cancer and 
#' intervention ages as the lower bound (see \code{\link{.currentAgeImpute}} 
#' documentation). If a lower imputation bound cannot be found using this 
#' approach and the individual has children in the pedigree, we assign their 
#' `CurAge` to be `15` plus the maximum difference in their living children's 
#' ages in the current imputation. (If the individual has one or zero living 
#' children, `CurAge` is simply set to `15`.) In cases where the individual 
#' has no children in the pedigree, their `CurAge` is set to `1`. 
#' 
#' @param ped A pedigree data frame returned by 
#' \code{\link{.pedStructureCheck}}. 
#' @param rel_l The list of relationships in `ped` returned by 
#' \code{\link{.pedStructureCheck}}. 
#' @param proband A numeric value or vector of the unique IDs in `ped` 
#' corresponding to the proband(s). 
#' @param impute_times The number of iterations to be imputed for missing ages. 
#' @param max_iter_tries A numeric value indicating the maximum number of 
#' iterations that should be attempted when multiply imputing missing ages. 
#' @param penet_db The cancer penetrances for carriers returned by 
#' \code{\link{buildDatabase}} in the `penet$penet_c` list component. 
#' @param dist_age A list containing empirical distributions of female and male 
#' ages at first birth; female and male birth ages; and sibling age 
#' differences. The default is \code{\link{AgeDistribution}}. 
#' @param penet_cbc_db A numeric array of contralateral cancer penetrance for 
#' carriers returned by \code{\link{buildDatabase}} in the `penet$penet_cbc` 
#' list component. 
#' @return A list where each component is a matrix of imputed ages 
#' corresponding to `CurAge` (if one or more individuals in `ped` is missing 
#' their current age) or `AgeX`, where `X` is a cancer for which one or more 
#' individuals in `ped` is missing the affection age. Each matrix has 
#' `impute_times` rows; the column names are the IDs of the individuals for 
#' whom ages were imputed. 
#' @family impute
ageImpute <- function(ped, rel_l, proband, impute_times, max_iter_tries, 
                      penet_db, dist_age = AgeDistribution, penet_cbc_db = NULL) {
  
  # Ensure that the cancer penetrance array is non-empty
  stopifnot(!is.null(penet_db))

  # Short cancer names in pedigree
  cancer_short_names <- intersect(
    .getCancersFromFam(ped),
    .mapCancerNames(long = dimnames(penet_db)$Cancer)
  )
  
  # add CBC if it is in the pedigree (it has its own penet_cbc object)
  if("CBC" %in% .getCancersFromFam(ped)){
    cancer_short_names <- c(cancer_short_names,"CBC")
  }
  
  # For each cancer, extract IDs of those with missing affection age
  m_aff <- list()
  for (cancer in cancer_short_names) {
    m_aff[[cancer]] <- ped$ID[(ped[[paste0("isAff", cancer)]] == 1) & 
                                (ped[[paste0("Age", cancer)]] == -999)]
  }
  
  # IDs of those with missing current ages who are alive
  m_cur <- ped$ID[ped$isDead == 0 & ped$CurAge == -999]
  # IDs of those with missing current ages who are dead
  m_cur_dead <- ped$ID[ped$isDead == 1 & ped$CurAge == -999]

  # Initialize list container for storing matrices of imputed ages
  lm <- list()
  # Initialize matrix for storing imputed current ages
  lm$CurAge <- matrix(
    nrow = impute_times, ncol = length(c(m_cur, m_cur_dead)),
    dimnames = list(seq(impute_times), sort(c(m_cur, m_cur_dead)))
  )
  # Initialize matrices for storing imputed cancer ages
  for (i in seq_along(m_aff)) {
    lname <- paste0("Age", names(m_aff)[i])
    lm[[lname]] <- matrix(
      nrow = impute_times,
      ncol = length(m_aff[[i]]),
      dimnames = list(seq(impute_times), m_aff[[i]])
    )
  }

  # Get first and second-degree relatives with respect to the proband
  # If there is more than one proband, use the one that appears first (smallest ID)
  proband_idx <- which(ped$ID == sort(proband)[1])
  first_degree <- ped$ID[c(proband_idx, .firstDegreeRelative(ped, proband_idx)$index)]
  second_degree <- setdiff(
    ped$ID[.secondDegreeRelative(ped, proband_idx)$index],
    first_degree
  )
  
  # Initialize vector for counting the number of failed imputations 
  # associated with each ID
  cancer_age_fails = rep(0, length(ped$ID))
  names(cancer_age_fails) = ped$ID
  
  # Initialize the failure flag, number of attempts, and index of imputations
  iter_failed = FALSE
  iter_tries = 0
  tt = 1
  while (tt <= impute_times) {
    # Create a temporary copy of the pedigree, to be updated with imputed ages
    ped_temp <- ped

    # First-degree relatives with missing current ages
    m_cur_first_degree_rels <- intersect(first_degree, m_cur)
    # Second-degree relatives with missing current ages
    m_cur_second_degree_rels <- intersect(second_degree, m_cur)
    # More distant relatives with missing current ages
    m_cur_distant_rels <- setdiff(
      setdiff(m_cur, m_cur_first_degree_rels),
      m_cur_second_degree_rels
    )
    # List of all relatives with missing current ages
    m_cur_list <- list(
      m_cur_first_degree_rels,
      m_cur_second_degree_rels,
      m_cur_distant_rels
    )
    
    # Impute missing current ages for people who are alive
    # Impute ages of first-degree relatives first, then second-degree 
    # relatives, then more distant relatives
    for (mcur_temp in m_cur_list) { # Iterate through relative groups
      if (length(mcur_temp) > 1) {
        mcur_temp = sample(mcur_temp)
      }
      # Keep imputing current ages until their are no relatives in the set who 
      # still have missing ages
      while (length(mcur_temp) > 0) {
        # Get age-bounding information for all relatives in this set with 
        # missing current ages
        rel_ages_list <- lapply(mcur_temp, function(id) {
          .getAgeFromPed(ped_temp, id, rel_l)
        })
        # Get all of their age bounds
        ab_list <- lapply(rel_ages_list, function(rel_ages) {
          .boundAgeFromRelatives(
            rel_ages = rel_ages
          )
        })
        # Calculate the interval width of each relative's age bounds
        ab_widths <- sapply(ab_list, diff)
        
        # Identify the relative with the narrowest interval 
        # (ties are broken randomly)
        narrowest_ab <- ifelse(length(ab_widths) == 1, 1, 
                               which.min(ab_widths))
        id <- mcur_temp[narrowest_ab]

        # Impute the current age of the relative with the narrowest interval
        if (ab_list[[narrowest_ab]][1] == ab_list[[narrowest_ab]][2]) {
          # If the upper and lower age bounds are equal, use this value as the 
          # current age instead of imputing it. 
          rlang::inform(sprintf(
            "Upper and lower current age bounds for ID %s are both %s, so skipping this age imputation.",
            id, ab_list[[narrowest_ab]][1]
          ))
          cur_age <- ab_list[[narrowest_ab]][1]
        } else { 
          # Otherwise, impute the current age based on the age bounds and age 
          # bounding information
          cur_age <- .currentAgeImpute(
            rel_ages = rel_ages_list[[narrowest_ab]],
            ab = ab_list[[narrowest_ab]],
            dist_age = dist_age
          )
        }

        # Update the temporary pedigree and matrix container for current ages
        ped_temp[["CurAge"]][ped_temp$ID == id] <- cur_age
        lm$CurAge[tt, as.character(id)] <- cur_age

        # Remove this ID from the vector of relative with missing current age
        mcur_temp <- mcur_temp[mcur_temp != id]
        if (length(mcur_temp) > 1) {
          mcur_temp = sample(mcur_temp)
        }
      }
    }
    
    # Impute missing cancer ages
    # Iterate through relatives with one or more missing cancer ages
    for (id in unique(unlist(m_aff))) {
      # Identify cancers with missing ages for this relative
      cancers_missed <- sort(names(m_aff)[sapply(m_aff, function(vec) id %in% vec)])
      # Information for extracting their cancer penetrances
      fnd_cols <- c("Sex", "race", "Ancestry", "CurAge")
      if("AgeBC" %in% colnames(ped_temp)){ fnd_cols <- c(fnd_cols, "AgeBC")}
      if("ER" %in% colnames(ped_temp)){ fnd_cols <- c(fnd_cols, "ER")}
      penet_info <- ped_temp[match(id, ped_temp$ID), fnd_cols]
      
      # Iterate through cancers with missing ages
      for (cancer in cancers_missed) {
        cancer_long <- .mapCancerNames(short = cancer)
        
        # Extract cancer penetrances for CBC
        if (cancer == "CBC" && ("BC" %in% cancer_short_names)) {
          cbc.params <- getCBCParams(ped_temp, id, penet_cbc_db, NULL)
          tmp.gene <- cbc.params$Gene
          if(tmp.gene == "noncarrier"){ # if negative for CBC genes or untested
            penet_vec <- CBCnoncarrierCrude(diag.age = cbc.params$FirstBCAge, 
                                            bc1.type = cbc.params$FirstBCType, 
                                            fdr.status = cbc.params$FDRBC, 
                                            ER.status = cbc.params$ER, 
                                            birth1.age = cbc.params$FirstBirth, 
                                            antiEst = cbc.params$AntiEstrogen, 
                                            HighRiskP = cbc.params$HRPreneoplasia, 
                                            BD = cbc.params$BreastDensity,
                                            bc1.t.size = cbc.params$BC1TSize,
                                            sex = cbc.params$Sex,
                                            race = cbc.params$race)
            
          } else { # if positive for a CBC associated gene
            penet_vec <- penet_cbc_db$cbc_carrier[cbc.params$Gene, cbc.params$Sex, 
                                                  cbc.params$FirstBCAge, "Crude",]
          }
          
          # Extract cancer penetrances for non-CBC cancers
        } else {
          gene.cols <- intersect(colnames(ped_temp), GENE_TYPES)
          if(length(gene.cols)>0){
            tmp.genes <- ped_temp[match(id, ped_temp$ID), gene.cols]
            if(sum(tmp.genes, na.rm = T) > 1){
              # if multiple positive genes, select the one with the highest penetrance at age 85
              tmp.genes <- tmp.genes[which(!is.na(tmp.genes) & tmp.genes==1)]
              tmp.genes <- unlist(DEFAULT_VARIANTS)[which(names(unlist(DEFAULT_VARIANTS)) == tmp.genes)]
              tmp.gene.85.penets <- penet_db[cancer_long, tmp.genes,
                                             penet_info$race, penet_info$Sex,
                                             85, "Crude"]
              tmp.gene <- names(which(tmp.gene.85.penets == max(tmp.gene.85.penets, na.rm = T)))
            } else if(sum(tmp.genes, na.rm = T) == 1){
              tmp.gene <- names(tmp.genes[which(!is.na(tmp.genes) & tmp.genes == 1)])
              tmp.gene <- unlist(DEFAULT_VARIANTS)[which(names(unlist(DEFAULT_VARIANTS)) == tmp.gene)]
            }
          } else {
            tmp.gene <- "SEER"
          }
          penet_vec <- penet_db[cancer_long, tmp.gene, penet_info$race,
                                penet_info$Sex, , "Crude"]
        }
        
        # Identify intervention age for preventative surgery, if relevant for 
        # the cancer
        riskMods <- unlist(ped_temp[ped_temp$ID == id, "riskmod"])
        intvAges <- unlist(ped_temp[ped_temp$ID == id, "interAge"])
        if ((cancer == "BC" | cancer == "CBC") && "Mastectomy" %in% riskMods) {
          interv_age <- intvAges[riskMods == "Mastectomy"]
        } else if (cancer == "OC" && "Oophorectomy" %in% riskMods) {
          interv_age <- intvAges[riskMods == "Oophorectomy"]
        } else if (cancer == "ENDO" && "Hysterectomy" %in% riskMods) {
          interv_age <- intvAges[riskMods == "Hysterectomy"]
        } else {
          interv_age <- -999
        }

        # The upper bound is the minimum of the individual's current age, the 
        # upper bound of this cancer's affection age, PPP:::MAXAGE 
        max_age_allowed <- .safeGet(min, c(
          ped_temp[["CurAge"]][ped_temp$ID == id],
          ifelse(
            ped_temp[ped_temp$ID == id, paste0("Age", cancer, "_upper")] == -999,
            MAXAGE,
            ped_temp[ped_temp$ID == id, paste0("Age", cancer, "_upper")]
          )
        ))
        
        # The lower bound is the maximum of the lower bound of this cancer's 
        # affection age, the age of preventative intervention (if relevant), 
        # and 1
        least_possible_age <- .safeGet(max, c(
          ped_temp[ped_temp$ID == id, paste0("Age", cancer, "_lower")],
          interv_age, 1
        ))
        
        # adjust min age for CBC
        if(cancer == "CBC"){
          if(least_possible_age < (penet_info$AgeBC+1)){
            least_possible_age <- penet_info$AgeBC+1
          }
        }
        
        # If the upper bound is less than the lower bound, swap them and 
        # print a message
        if (max_age_allowed < least_possible_age) { 
          rlang::inform(sprintf(
            "ID %s cannot get valid capping estimate for cancer %s, with upper bound %s and lower bound %s.",
            id, cancer, least_possible_age, max_age_allowed
          ))
          temp <- max_age_allowed
          max_age_allowed <- least_possible_age
          least_possible_age <- temp
        }
        
        # Impute the affection age for this cancer based on the age bounds
        cancer_age <- .cancerAgeImpute(
          max_age_allowed = max_age_allowed,
          penet = penet_vec,
          least_possible_age = least_possible_age
        )
        
        if (is.na(cancer_age) == TRUE) { # If a cancer age could not be imputed
          # Set failure flag to TRUE
          iter_failed = TRUE
          # Keep track of whose fault it was
          cancer_age_fails[as.character(id)] = cancer_age_fails[as.character(id)] + 1
          # Exit loop
          break
        } else { # Otherwise
          # Update the temporary pedigree and matrix container for this cancer's 
          # affection ages
          ped_temp[[paste0("Age", cancer)]][ped_temp$ID == id] <- cancer_age
          lm[[paste0("Age", cancer)]][tt, as.character(id)] <- cancer_age
        }
      }
      
      # Exit loop if the imputation failed
      if (iter_failed == TRUE) {
        break
      }
    }
    
    # Increment number of tries
    iter_tries = iter_tries + 1
    # Check if the number of tries exceeds the maximum allowed
    if (iter_tries == max_iter_tries) {
      # Warning message about which individuals could not have ages imputed
      msg_fails = cancer_age_fails[cancer_age_fails!=0]
      msg = paste("Imputing ages failed during", 
                  paste(msg_fails, "iterations for", 
                        "ID", names(msg_fails), ";", 
                        collapse = " "))
      # Error message that the maximum number of tries has been exceeded
      rlang::abort(sprintf(
        "Could not generate %s age imputations after the maximum of %s tries. \n %s \n Consider including more age information in the pedigree.", 
        impute_times, max_iter_tries, msg))
    }
    # Skip to the top of the loop and reset the failure flag 
    if (iter_failed == TRUE) {
      iter_failed = FALSE
      next
    }

    # Impute current age (age of death) for dead people if at least one
    # cancer affection or intervention age is known
    # Iterate through dead relatives with missing current ages
    for (id in m_cur_dead) {
      # Get age-bounding information
      rel_ages <- .getAgeFromPed(ped_temp, id, rel_l)
      
      # IDs of children
      cids = rel_l[[which(ped_temp$ID == id)]]$cid

      # The lower bound is the maximum of their cancer affection ages and any 
      # intervention ages
      ab_l <- .safeGet(max, c(rel_ages$max_aff, rel_ages$interv))

      # If a non-missing lower age bound can be found, use it to impute age of 
      # death
      if (!is.na(ab_l)) {
        # Set all relative ages to NA in the age-bounding information, since 
        # this is not informative for setting bounds for age of death
        rel_ages$mage <- NA
        rel_ages$fage <- NA
        rel_ages$cage <- NA
        rel_ages$sage <- NA

        # Impute current age (age of death) based on lower bound
        cur_age <- .currentAgeImpute(
          rel_ages = rel_ages,
          ab = c(ab_l, MAXAGE),
          dist_age = dist_age
        )
      } else if (length(cids) > 0) {
        # Otherwise, check if the person has children
        
        # If the person has children, set their death age to 15
        cur_age = 15
        
        # If the person has more than one living child, set the death age 
        # to 15 plus the maximum difference between the children's ages, 
        # but no more than MAXAGE
        children_ages = ped_temp[ped_temp$ID %in% cids, "CurAge"]
        children_dead = ped_temp[ped_temp$ID %in% cids, "isDead"]
        children_ages = children_ages[children_dead != 1]
        if (length(children_ages) >= 2) {
          cur_age = min(cur_age + (max(children_ages) - min(children_ages)), 
                        MAXAGE)
          if (tt == 1) {
            rlang::inform(sprintf(
              "The CurAge (death age) for dead individual ID %s is missing and has been set to 15 plus the maximum difference between their living children's ages",
              id
            ))
          }
        } else {
          if (tt == 1) {
            rlang::inform(sprintf(
              "The CurAge (death age) for dead individual ID %s is missing and has been set to 15 (assumed that they lived long enough to have the children included in the pedigree)",
              id
            ))
          }
        }
      } else {
        # Otherwise, set the person's death age to 1
        cur_age = 1
        if (tt == 1) {
          rlang::inform(sprintf(
            "The CurAge (death age) for dead and childless individual ID %s is missing and has been set to 1",
            id
          ))
        }
      }
      
      # Update the temporary pedigree and matrix container for current ages
      ped_temp[["CurAge"]][ped_temp$ID == id] <- cur_age
      lm$CurAge[tt, as.character(id)] <- cur_age
    }
    
    # Increment the number of (successful) imputations
    tt = tt + 1
  }
  
  if (sum(cancer_age_fails)) {
    # Warning message about which individuals could not have ages imputed
    msg_fails = cancer_age_fails[cancer_age_fails!=0]
    msg = paste("Imputing ages failed during", 
                paste(msg_fails, "iterations for", 
                      "ID", names(msg_fails), ";", 
                      collapse = " "))
    # Error message that the maximum number of tries has been exceeded
    rlang::warn(sprintf(
      "%s \n Consider including more age information in the pedigree.", 
      msg))
  }
  
  return(lm)
}
