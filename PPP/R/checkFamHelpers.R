#' Check pedigree structure and connectivity
#' 
#' Checks that all required columns in pedigree are present, valid, and 
#' consistent, printing informative error/warning messages and making 
#' corrections when possible. It also removes family members that are not 
#' connected to the proband(s) and splits the pedigree into a list of connected 
#' sub-families, if necessary. A list of relationships for each connected 
#' pedigree is also returned. 
#'
#' @details 
#' This function performs the following structural and connectivity checks, 
#' issuing informative error/warning messages and making corrections when 
#' possible: 
#' * Check expected columns
#'   * Check that all required columns (`ID`, `Sex`, `MotherID`, `FatherID`, 
#'   `CurAge`, `Twins`, `isDead`, `race`, `riskmod`, `interAge`, `Ancestry`) 
#'   are present in `ped`, and raise a `ColumnNameNotMatch` error if not. 
#' * Check probands
#'   * If a non-`NULL` proband argument was specified, the `proband` argument 
#'   will override the contents of the `ped$isProband` column. If 
#'   `proband == "All"`, all IDs in the pedigree will be used. Otherwise, check 
#'   if all IDs in `proband` are present in the pedigree, and raise a 
#'   `ProbandNotValid` error if not. If the `proband` argument does not match 
#'   with the probands implied by `ped$isProband`, print a 
#'   `ProbandOverwrite` message informing the user that the probands in 
#'   `ped$isProband` have been overridden by `proband`. 
#'   * If the `proband` argument is `NULL` or unspecified by the user, search 
#'   `ped` for the `isProband` column, and raise a `ProbandNotValid` error if 
#'   it is not present in the data frame. Then, attempt to infer probands from 
#'   the `isProband` column, treating `NA` values as `0` (not probands). If no 
#'   probands (encoded as `1`) are specified by `isProband`, raise a 
#'   `ProbandNotValid` error. If `isProband` contains unexpected values, raise a 
#'   `TypeCheck` error. 
#'   * Raise a `ProbandMissingAge` error if any proband's `CurAge` is missing. 
#'   * Raise a `ProbandAgeOutOfBounds` warning if any proband's `CurAge` is 
#'   equal to PPP:::MAXAGE, in which case their future risk will not be 
#'   estimated. 
#' * Check cancers and ages
#'   * Columns in `ped` with the naming convention `isAffX` are used to encode 
#'   cancer affection status and columns with the naming convention `AgeX` are
#'    used to encode affection ages. `X` should be one of the cancer short 
#'    names in `PPP::CANCER_NAME_MAP$short`. 
#'   * Raise a `CancerNotFound` error if no `isAffX` columns can be found in 
#'   `ped`. 
#'   * Raise a `BCNotExist` error if CBC specified in the model but `isAffBC` or
#'   `AgeBC` columns cannot be found.
#'   * Raise a `CBC_BCAge` error if at least one subject has a BC age greater 
#'   than their CBC Age.
#'   * Raise a `CBC_Risk_Indicator` error if codings of `FirstBCType`, `AntiEstrogen`, 
#'   `HRPreneoplasia`, and/or `BreastDensity` columns have unrecognized values.
#'   * Raise a `ModelCancersNotExist` error if `ped` does not contain `isAffX` 
#'   columns for all of the cancers in in the model specification `cancers`. 
#'   * Raise an `AffectionAgeNotFound` error if any cancer does not have an 
#'   `AgeX` column. 
#'   * If the `keep.only.model.cancers` argument is set to `TRUE` (the default 
#'   behavior), remove the `isAffX` and `AgeX` columns from `ped` 
#'   corresponding to cancers that are not in `cancers`. 
#'   * Values other than `0` and `1` in the `isAffX` columns will be recoded 
#'   as `-999`, accompanied by an `InvalidNA` message.  
#'   * If any individual is unaffected for a cancer but reports an affection 
#'   age for the cancer, the affection age will be set to `NA` and a
#'   `CancerAgeMisSpecified` message will be raised. 
#'   * If any individual develops cancer after their current age, a 
#'   `CancerAfterCurrentAge` error will be raised along with a suggestion to 
#'   substitute the maximum cancer age for the current age, if the user is 
#'   unable to resolve the age inputs otherwise. 
#'   * If the `allow.age.zero` argument is set to `FALSE` (the default 
#'   behavior), set values in the `AgeX` and `CurAge` columns that are `0` to 
#'   `NA` and print a `FalseAgeZero` message. These values will be treated as 
#'   missing instead of `0`. 
#'   Values of `AgeX` and `CurAge` between `0` and `1` will be rounded up to 
#'   `1`; other decimal ages will be rounded to the nearest integer. In both 
#'   cases, an `AgeRounding` message will be printed. 
#' * Check IDs and personal information
#'   * IDs should be unique identifiers in the pedigree. If any IDs in the 
#'   `ped$ID` column are duplicated, raise a `PedIDInvalid` error. 
#'   * Check that all IDs in the `MotherID` and `FatherID` columns can be 
#'   found in `ped$ID`. Set parent IDs that cannot be found in the `ID` to be 
#'   `NA`, so that they will be treated as missing from the pedigree, and raise 
#'   an `IDNotExist` warning. 
#'   * Set missing or unsupported values of `ped$race` and `ped$Ancestry` to 
#'   the default race and ancestry specified by `unknown.race` and 
#'   `unknown.ancestry`, and raise an `InvalidNA` message.  
#'   * Set unexpected values of `ped$Sex` to `NA` and raise an 
#'   `InvalidInputForcingNA` warning. 
#'   * If `ped$isDead` is unknown (`NA`), assume that the individual is alive, 
#'   set the value to `0`, and raise an `InvalidNA` message. 
#'   * Check that sexes of parents are valid, i.e. that all IDs in `MotherID` 
#'   are female and all IDs in `Father ID` are male. Raise an `InvalidParentSex` 
#'   error if this is not the case. 
#' * Check identical twins
#'   * If `ped$Twins` is unknown (`NA`), assume that the individual is not an 
#'   identical twin, set the value to `0`, and raise an `InvalidNA` message.   
#'   * The following steps describe consistency checks among homozygotic twins 
#'   and other identical multiple births. For more details, 
#'   see the \code{\link{.checkTwins}} documentation. 
#'   * If only one twin is present in a twin set encoded by `ped$Twins`, remove 
#'   the individual's twin status by setting it to `0` and print a 
#'   `MissingTwinPair` message 
#'   * Assume that twins must have the same `FatherID`, `MotherID`, and `Sex`. 
#'   If this information is missing for one (or more) twins but present in the 
#'   other(s), infer the missing information from the individuals where it is 
#'   known and print a `MissingTwinInfoModified` message. If the twins have 
#'   unresolvable mismatched information, remove the twin label and raise a 
#'   `TwinInfoConflict` warning. 
#'   * Assume that living twins must have the same `CurAge`. If there is a 
#'   disagreement in `CurAge` among living twins, the proband's `CurAge` 
#'   overrides the others if one of the twins is a proband, and a 
#'   `TwinInfoConflict` warning will be raised. Otherwise, the average `CurAge` 
#'   will be assigned to all twins, and a `TwinInfoConflict` warning will be 
#'   raised. 
#'   * If the `CurAge` of a living twin is less than the `CurAge` of a dead 
#'   twin, the average `CurAge` will be assigned to all twins, and a 
#'   `TwinInfoConflict` warning will be raised. 
#'   * Assume that twins must have the same `race` and `Ancestry`. If there is 
#'   a disagreement, assign the default race and ancestry specified  by 
#'   `unknown.race` and `unknown.ancestry` to all twins, and raise a 
#'   `TwinInfoConflict` warning. 
#' * Check risk modifiers/interventions
#'   * Individuals who report `NA` for `riskmod` and `interAge` are assumed to 
#'   have no risk modifiers, and an `InvalidInputForcingNone` warning will be 
#'   raised. 
#'   * If an individual's `riskmod` and `interAge` are of different lengths, 
#'   the longer will be truncated to match the length of the shorter. A 
#'   `RiskmodInconsistency` message will be raised. 
#'   Intervention ages between `0` and `1` will be rounded up to `1`; other 
#'   decimal ages will be rounded to the nearest integer. In both cases, an 
#'   `AgeRounding` message will be printed. 
#'   * Unsupported risk modifiers and invalid intervention ages will be removed 
#'   and an `InvalidInputForcingNone` message will be printed. 
#'   * Raise a `DuplicatedRiskmod` error if an individual has duplicated risk 
#'   modifiers stored in `riskmod`. 
#'   * Raise an `InvalidRiskmod` error if an individual's intervention is not 
#'   preventative for the relevant cancer. 
#'   * For individuals who do not report any risk modifiers in `ped`, set their 
#'   value of `riskmod` to `"None"` and `interAge` to `1`. 
#' * Check tumor biomarker testing
#'   * If an individual reports marker testing results in `ped`, but they are 
#'   not affected for the related cancer, remove their marker testing 
#'   results and raise an `InvalidTests` warning. 
#'  * If `ignore.proband.germ = TRUE`, remove proband(s)'s germline testing 
#'  information from `ped`. 
#' * Check germline testing
#'   * If germline testing results are reported in `ped`, print a 
#'   `GermlineTestAssumption` message to inform the user that the default 
#'   variant (see `PPP:::DEFAULT_VARIANTS`) will be assumed for each 
#'   mutation. 
#' * Typechecks
#'   * If any value of `AgeX` or `CurAge` exceeds `PPP:::MAXAGE`, set the 
#'   values to `PPP:::MAXAGE` and raise an `ExceedsMAXAGE` warning. 
#'   * Raise a `TypeCheck` error if `Sex`, `Twins`, or `isAffX` is missing 
#'   (`NA`) for any individual. 
#'   * Convert `NA` values to `-999` and raise an `InvalidNA` message.  
#' * Check connectivity 
#'   * Attempt to correct inconsistent hereditary features (`race` and 
#'   `Ancestry`) among parents and children, and print a `HeredityFix` message. 
#'   For more details, see the \code{\link{.checkHeredity}} documentation. 
#'   * Remove family members who are disconnected from the proband(s), raising 
#'   a `DisconectedFamilyRemoval` warning if so. For more details, see the 
#'   \code{\link{.removeDisconnected}} documentation. 
#'   * Resolve inconsistent sex information. If an individual's sex is 
#'   incompatible with their parental identity, set the corresponding parent ID 
#'   for their children to be missing. If an individual's sex is incompatible 
#'   with their cancer affection status, set the corresponding `isAffX` 
#'   value to `0` and raise a `GenderCancerMismatch` warning. For more details, 
#'   see the \code{\link{.checkGender}} documentation. Finally, encode the 
#'   `Sex` column as `"Female"`/`"Male"` instead of `0`/`1`. 
#'   * Validate the `AgeX_lower` and `AgeX_upper` columns for bounding cancer 
#'   affection age, or create them if necessary. If a lower bound is greater 
#'   than the corresponding upper bound, raise an `AgeBounds` error. For more 
#'   details, see the \code{\link{.boundAffectionAges}} documentation. 
#'   * Check that `ped` does not have any invalid mating patterns, i.e. 
#'   intermarriages that create "loops", and raise a `BadMating` error if not. 
#'   For more details, see the \code{\link{.checkMating}} documentation. 
#'   
#' @param ped A pedigree data frame. 
#' @param cancers A vector of short cancer names corresponding to the cancers 
#' specified by the model that should be checked in `ped`. See 
#' `PPP:::CANCER_NAME_MAP$short` for the list of options. The default is 
#' `.getCancersFromFam(ped)`, which is a vector of all the short cancer names 
#' included in `ped`.  
#' @param proband A numeric value or vector of the unique IDs in `ped` 
#' corresponding to the proband(s). The default is `NULL`, in which case 
#' probands will be inferred from the `isProband` column in `ped`, if it exists 
#' and is valid. Otherwise, the IDs in `proband` will override the contents of 
#' `ped$isProband`. The user can also specify `proband = "All"`, in which case 
#' all IDs in the `ped` will be treated as probands. 
#' @param unknown.race A character string indicating the default race to use 
#' when `ped$race` is missing or unsupported. The default is 
#' `PPP:::UNKNOWN_RACE`. 
#' @param unknown.ancestry A character string indicating the default ancestry 
#' to use when `ped$Ancestry` is missing or unsupported. The default is 
#' `PPP:::UNKNOWN_ANCESTRY`. 
#' @param ignore.proband.germ A logical value indicating if the proband(s)'s 
#' germline testing results (if provided in `ped`) should be considered. The 
#' default is `FALSE`, in which case the proband(s)'s germline testing results 
#' will not be ignored. 
#' @param allow.age.zero A logical value indicating whether zero ages can be 
#' reported in `ped`. The default is `FALSE`, in which case `0` values in age 
#' columns will be treated as missing. `allow.age.zero` should only be `TRUE` 
#' if individuals are confirmed to born with a cancer. 
#' @param keep.only.model.cancers A logical value indicating if only the 
#' cancers in `cancers` should be kept in `ped`. The default is `TRUE`. 
#' @family check
#' @seealso \code{\link{checkFam}}, \code{\link{.checkTwins}}, 
#' \code{\link{.checkHeredity}}, \code{\link{.removeDisconnected}}, 
#' \code{\link{.checkGender}}, \code{\link{.boundAffectionAges}}, 
#' \code{\link{.checkMating}}
#' @return A list with three components: 
#' * `ped_list`: A list where each element is a connected pedigree that is a 
#' subset of the relatives in `ped`. Typically, `ped_list` has length `1`, 
#' unless there are multiple probands who are related to multiple disconnected 
#' pedigrees. Relatives who are not connected to the proband(s) will be removed. 
#' * `proband`: A numeric value or vector of the unique IDs in the pedigree
#' corresponding to the proband(s). If the `proband` argument was specified by 
#' the user, the returned list component will be identical. Otherwise, the IDs 
#' will correspond to individuals in the `ped` argument where 
#' `ped$isProband == 1`. 
#' * `rel_ls`: A list of the same length as `ped_list`, where each component is 
#' a list of relationships for each pedigree in `ped_list`. 
.pedStructureCheck <- function(ped, cancers = .getCancersFromFam(ped),
                               proband = NULL,
                               unknown.race = UNKNOWN_RACE,
                               unknown.ancestry = UNKNOWN_ANCESTRY,
                               ignore.proband.germ = FALSE, 
                               allow.age.zero = FALSE,
                               keep.only.model.cancers = TRUE) {
  
  # Extract column names in pedigree
  famcols <- colnames(ped)
  
  # ----- CHECK EXPECTED COLUMNS -----
  err_code <- "ColumnNameNotMatch"
  
  # Required columns (two riskmod columns)
  expected_colnames1 <- c("ID", "Sex", "MotherID", "FatherID", "CurAge", "Twins",
                         "isDead", "race", "riskmod", "interAge", "Ancestry")
  ifcontains_colnames1 <- expected_colnames1 %in% famcols
  
  # Required columns (six riskmod columns)
  expected_colnames2 <- c("ID", "Sex", "MotherID", "FatherID", "CurAge", "Twins",
                         "isDead", "race", "riskmodHyst", "riskmodOoph", 
                         "riskmodMast", "interAgeHyst", "interAgeOoph",
                         "interAgeMast", "Ancestry")
  ifcontains_colnames2 <- expected_colnames2 %in% famcols
  
  # Error if any of these columns are not present in the pedigree
  if (any(!ifcontains_colnames1) & any(!ifcontains_colnames2)) {
    # One of the two riskmod columns are missing  
    if(sum(ifcontains_colnames1) > sum(ifcontains_colnames2)) {
      which_not_present <- expected_colnames1[!ifcontains_colnames1]
      msg <- sprintf(" %s are expected column names but cannot be found.", 
                           paste(which_not_present, collapse = ","))
      rlang::abort(msg, level = err_code)
    }
    # One or more of the six riskmod columns are missing 
    if(sum(ifcontains_colnames1) < sum(ifcontains_colnames2)) {
      which_not_present <- expected_colnames2[!ifcontains_colnames2]
      msg <- sprintf(" %s are expected column names but cannot be found.", 
                                 paste(which_not_present, collapse = ","))
      rlang::abort(msg, level = err_code)
    }
    # No riskmod columns (default: ask for two) or missing other columns  
    if(sum(ifcontains_colnames1) == sum(ifcontains_colnames2)) {
      which_not_present <- expected_colnames1[!ifcontains_colnames1]
      msg <- sprintf(" %s are expected column names but cannot be found.", 
                                 paste(which_not_present, collapse = ","))
      rlang::abort(msg, level = err_code)
    }     
  }
  
  # ----- CHECK PROBANDS -----
  err_code <- "ProbandNotValid"
  
  if (is.null(proband)) { # If no probands were specified by the user
    if (!"isProband" %in% famcols) {
      # Error if the pedigree does not have an isProband column
      msg <- "\"isProband\" does not exist as column of the pedigree dataframe and `proband` is not specified."
      rlang::abort(msg, level = err_code)
    } else {
      # Assign NA values in the isProband column to 0
      ped$isProband[is.na(ped$isProband)] <- 0
      # Error if no probands were specified by the isProband column
      if (sum(ped$isProband) == 0) {
        msg <- "\"isProband\" column exists but it is not informative for proband identification and `proband` is not specified."
        rlang::abort(msg, level = err_code)
      }
      # Error if the isProband column does not only contain zeros and ones
      if (!all(ped$isProband %in% c(0, 1))) {
        rlang::abort("isProband column should only contain zeros and ones", 
                     level = "TypeCheck")
      }
      # Assign proband to IDs corresponding to ped$isProband==1
      proband <- ped$ID[ped$isProband == 1]
    }
  } else { # Otherwise, use the information in the proband argument
    if (length(proband) == 1 && proband == "All") {
      # If proband=="All", use all IDs in the pedigree
      proband = ped$ID
      ped$isProband = 1
      msg <- sprintf("`proband` argument is set to 'All', so overwriting with all IDs in pedigree")
      rlang::inform(msg, level = "ProbandOverwrite")
    } else if (!all(proband %in% ped$ID)) {
      # Error if any of probands are not present
      not_in_ids <- proband[!proband %in% ped$ID]
      msg <- sprintf("IDs {%s} in `proband` argument do not exist in the pedigree", 
                     paste0(not_in_ids, collapse = ","))
      rlang::abort(msg, level = err_code)
    }
    
    # Assign NA values in the isProband column to 0
    ped$isProband[is.na(ped$isProband)] <- 0
    # Identify IDs corresponding to ped$isProband==1
    pids <- ped$ID[ped$isProband == 1]
    # Warning if the proband argument does not match with the isProband column
    if ("isProband" %in% famcols && !setequal(pids, proband)) {
      msg <- sprintf("`proband` argument does not match with \"isProband\" column in the pedigree, using `proband`")
      rlang::inform(msg, level = "ProbandOverwrite")
      ped$isProband <- as.numeric(ped$ID %in% proband)
    }
  }
  
  # Error if the proband(s)' current age(s) are not known
  if (any(is.na(ped[ped$ID %in% proband, "CurAge"]))) {
    rlang::abort("CurAge of proband(s) must be provided", 
                 level = "ProbandMissingAge")
  }
  
  # Warning that future risk will not be calculated if proband(s)' current 
  # age is the max age
  proband_maxage = ped$ID[which(ped$isProband == 1 & ped$CurAge == MAXAGE)]
  if (length(proband_maxage) > 0) {
    msg <- sprintf("Proband ID(s) %s have current age equal to the maximum age of %s, so their future risks will not be estimated.", 
                   paste(proband_maxage, collapse = ","), MAXAGE)
    rlang::warn(msg, level = "ProbandAgeOutOfBounds")
  }
  
  # ----- CHECK CANCERS AND AGES -----
  err_code <- "CancerNotMatch"
  
  # Extract cancers from the pedigree
  affect_cancers <- .getCancersFromFam(ped)
  
  if ("CBC" %in% cancers & 
      (!"BC" %in% affect_cancers | !"AgeBC" %in% colnames(ped))){
    rlang::abort(sprintf(
      "CBC is in the model but pedigree is missing isAffBC and/or AgeBC"), 
      level = "BCNotExist")
  }
  
  # check for contralateral breast cancer
  if (all(c("BC", "CBC") %in% affect_cancers) &&
      all("AgeCBC" %in% colnames(ped)) &&
      (sum(ped$isAffCBC == 1, na.rm = TRUE) > 0)) {
    
    # check for subjects affected with CBC but not BC
    if(length(which(ped$isAffBC == 0 & ped$isAffCBC == 1)) > 0){
      rlang::abort(sprintf(
        "At least one subject has CBC but not BC"), 
        level = "BCNotExist")
    }
    
    # check if CBC ages are all greater than BC ages
    if(sum(!(ped$AgeCBC > ped$AgeBC), na.rm = T) > 0){
      rlang::abort(sprintf(
        "At least one subject has a BC age greater than their CBC Age"), 
        level = "CBC_BCAge")
    }
    
    # change invalid CBC entries to 0 and NA
    CBC_valid <- which((ped$isAffBC * ped$isAffCBC) == 1 &
                       !.forceFalse(ped$AgeCBC >= ped$AgeBC) == 1)
    ped$isAffCBC[-CBC_valid] <- 0
    ped$AgeCBC[-CBC_valid] <- NA
    
    # check for CBC specific columns relevant for non-carriers
    cbc_risk_indicators <- c("FirstBCType", "AntiEstrogen", "HRPreneoplasia", 
                             "BreastDensity", "FirstBCTumorSize")
    cbc_risk_indicators <- cbc_risk_indicators[cbc_risk_indicators %in% famcols]
    if(length(cbc_risk_indicators) > 0){
      for(cri in cbc_risk_indicators){
        if(cri %in% c("AntiEstrogen", "HRPreneoplasia")){
          if(any(!ped[, cri] %in% c(NA, 0, 1)) == TRUE){
            rlang::abort(sprintf(
              paste0(cri, " column must only contain values NA, 0, or 1 which
                     correspond to unknown, no, or yes respectively.")), 
              level = "CBC_Risk_Indicator")
          }
        } else if(cri == "BreastDensity"){
          if(any(!ped$BreastDensity %in% c(NA, "a", "b", "c", "d")) == TRUE){
            rlang::abort(sprintf(
              "Breast Density column must only contain values NA or the letters 
              a to d which correspond with BI-RADS breast density descriptions."), 
              level = "CBC_Risk_Indicator")
          }
        } else if(cri == "FirstBCType"){
          if(any(!ped$FirstBCType %in% c(NA, "Invasive", "Invasive_DCIS")) == TRUE){
            rlang::abort(sprintf(
              "FirstBCType column must only contain values NA, Invasive, or Invasive_DCIS."), 
              level = "CBC_Risk_Indicator")
          }
        } else if(cri == "FirstBCTumorSize"){
          if(any(!ped$FirstBCTumorSize %in% c(NA, "T0/T1/T2", "T3/T4", "Tis")) == TRUE){
            rlang::abort(sprintf(
              "FirstBCTumorSize column must only contain values NA, 'T0/T1/T2', 'T3/T4', or 'Tis'."), 
              level = "CBC_Risk_Indicator")
          }
        }
      }
    }
  }
  
  # Error if no isAff* columns can be found in the pedigree
  if (length(affect_cancers) == 0 && length(cancers) > 0) {
    rlang::abort(sprintf(
      "Please use \"isAff*\" to denote cancer affection status. e.g. isAffBC for breast cancer"), 
      level = "CancerNotFound")
  }
  
  # Error if pedigree does not contain all cancers in the model specification
  if_contains_MS_cancers <- cancers %in% affect_cancers
  if (any(!if_contains_MS_cancers)) {
    rlang::abort(sprintf(
      "%s cancer is required in model but cannot be found in the pedigree", 
      paste0(cancers[!if_contains_MS_cancers], collapse = ",")), 
      level = "ModelCancersNotExist")
  }
  
  # Extract cancer affection and age columns
  if (length(cancers) > 0) {
    expected_cancer_age_cols <- paste0("Age", cancers)
    expected_cancer_cols <- paste0("isAff", cancers)
    if_have_cancerage_cols <- expected_cancer_age_cols %in% famcols
  } else {
    expected_cancer_age_cols <- c()
    expected_cancer_cols <- c()
    if_have_cancerage_cols <- c()
  }
  
  # Error if any of the cancers do not have affection age columns in pedigree
  if (length(cancers) > 0 && any(!if_have_cancerage_cols)) {
    rlang::abort(sprintf(
      "%s cannot be found in pedigree", 
      paste0(expected_cancer_age_cols[!if_have_cancerage_cols], collapse = ","),
      level = "AffectionAgeNotFound"
    ))
  }
  
  # If only keeping model cancers in the pedigree, drop other cancer 
  # affection and age columns
  if (keep.only.model.cancers) {
    cancers_to_remove <- affect_cancers[!affect_cancers %in% cancers]
    ped[, paste0("isAff", cancers_to_remove)] <- NULL
    ped[, paste0("Age", cancers_to_remove)] <- NULL
    affect_cancers <- .getCancersFromFam(ped)
  }
  
  # Values other than c(0,1) in cancer affection columns will be recoded as 
  # -999
  for (col in expected_cancer_cols) {
    ped[[col]] <- forcingNA_ifnot_contains(ped, col, c(1, 0), -999)
  }
  
  ncancers <- length(expected_cancer_cols)
  if (ncancers > 0) {
    for (i in seq(ncancers)) {
      # Force affection age to NA for those who are unaffected for the cancer
      aff <- ped[[expected_cancer_cols[i]]]
      age <- ped[[expected_cancer_age_cols[i]]]
      if (any(.forceFalse(aff == 0) & !is.na(age))) {
        who <- ped$ID[.forceFalse(aff == 0) & !is.na(age)]
        msg <- sprintf(
          "ID %s is unaffected for %s, forcing affection age to NA.",
          paste0(who, collapse = ","), cancers[i]
        )
        rlang::inform(msg, level = "CancerAgeMisSpecified")
        age[ped$ID %in% who] <- NA
        ped[expected_cancer_age_cols[i]] <- age
      }
      
      isCancerAgeAboveCurAge = ped[expected_cancer_age_cols[i]] > ped$CurAge
      if (any(isCancerAgeAboveCurAge, na.rm = TRUE)) {
        who <- ped$ID[isCancerAgeAboveCurAge == TRUE & !is.na(isCancerAgeAboveCurAge)]
        msg <- sprintf(
          "The %s age for ID %s occurs after their current age. If you cannot resolve the age inputs, consider using the individual's maximum cancer affection age as their current age.",
          cancers[i], paste0(who, collapse = ",")
        )
        rlang::abort(msg, level = "CancerAfterCurrentAge")
      }
    }
    
    # All age columns
    age_cols <- c(expected_cancer_age_cols, "CurAge")
    
    # If zero ages are not allowed, set them to NA
    if (!allow.age.zero) {
      if (any(ped[age_cols] == 0, na.rm = TRUE)) {
        msg <- "Forcing age = 0 to be NA. "
        rlang::inform(msg, level = "FalseAgeZero")
        # Note that fix_with_rep also replaces negative values
        ped[age_cols] <- lapply(ped[age_cols], fix_with_rep(0, NA))
      }
    }
    
    # Round non-integer ages to integers and issue warning (all values < 1 
    # get rounded up to 1)
    for (age_col in age_cols) {
      
      baby_idx = which(ped[age_col] > 0 & ped[age_col] < 1)
      if (length(baby_idx) > 0) {
        rlang::inform(sprintf(
          "ID %s have %s between 0 and 1; rounding up to 1.",
          paste(ped$ID[baby_idx], collapse = ","), age_col
        ), level = "AgeRounding")
        ped[baby_idx,age_col] = 1
      }
      
      dec_idx = which(ped[age_col] %% 1 != 0)
      if (length(dec_idx) > 0) {
        rlang::inform(sprintf(
          "ID %s have decimal-value %s; rounding to nearest integer.",
          paste(ped$ID[dec_idx], collapse = ","), age_col
        ), level = "AgeRounding")
        ped[dec_idx,age_col] = round(ped[dec_idx,age_col])
      }
    }
  }
  
  # ----- CHECK IDS AND PERSONAL INFORMATION  -----
  # Error if IDs are duplicated
  if (any(duplicated(ped$ID))) {
    rlang::abort(sprintf(
      "ID %s have multiple entries",
      paste(ped$ID[duplicated(ped$ID)], collapse = ",")
    ), level = "PedIDInvalid")
  }
  
  # Warning if any IDs in MotherID column cannot be found in the ID column
  if_mothers_contained <- ped$MotherID %in% ped$ID
  if (any(!if_mothers_contained)) {
    mout <- ped$MotherID[!if_mothers_contained]
    mout <- mout[!is.na(mout)]
    if (length(mout) != 0) {
      msg <- sprintf(
        "Mother ID %s cannot be found in pedigree, forcing it to NA.",
        paste0(unique(mout), collapse = ",")
      )
      rlang::warn(msg, "IDNotExist")
      ped$MotherID[ped$MotherID %in% mout] <- NA
    }
  }
  # Warning if any IDs in FatherID column cannot be found in the ID column
  if_fathers_contained <- ped$FatherID %in% ped$ID
  if (any(!if_fathers_contained)) {
    fout <- ped$FatherID[!if_fathers_contained]
    fout <- fout[!is.na(fout)]
    if (length(fout) != 0) {
      msg <- sprintf(
        "Father ID %s cannot be found in pedigree, forcing it to NA.",
        paste0(unique(fout), collapse = ",")
      )
      rlang::warn(msg, "IDNotExist")
      ped$FatherID[ped$FatherID %in% fout] <- NA
    }
  }
  
  # Set missing or unexpected race types to default race
  ped$race <- forcingNA_ifnot_contains(ped, "race", RACE_TYPES, unknown.race)
  # Set missing or unexpected ancestry types to default ancestry
  ped$Ancestry <- forcingNA_ifnot_contains(ped, "Ancestry", ANCESTRY_TYPES, 
                                           unknown.ancestry)
  # Set missing or unexpected sex types to NA
  ped$Sex <- forcingNA_ifnot_contains(ped, "Sex", c(1, 0))
  
  # When death status is not known, we assume they are alive
  ped$isDead <- forcingNA_ifnot_contains(ped, "isDead", c(1, 0), 
                                         default_NA_rep = 0)
  
  # Error if any of the IDs in MotherID are not female (0)
  moth_sexes <- ped$Sex[unlist(sapply(unique(ped$MotherID[!is.na(ped$MotherID)]), 
                                      function(id) {which(id==ped$ID)}))]
  if (any(moth_sexes != 0)) {
    mout <- unique(ped$MotherID[!is.na(ped$MotherID)])[moth_sexes != 0]
    mout <- mout[!is.na(mout)]
    if (length(mout) != 0) {
      msg <- sprintf(
        "Mother ID %s has invalid sex.",
        paste0(unique(mout), collapse = ",")
      )
      rlang::abort(msg, "InvalidParentSex")
    }
  }
    
  # Error if any of the IDs in FatherID are not male (1)
  fath_sexes <- ped$Sex[unlist(sapply(unique(ped$FatherID[!is.na(ped$FatherID)]), 
                                      function(id) {which(id==ped$ID)}))]
  if (any(fath_sexes != 1)) {
    fout <- unique(ped$FatherID[!is.na(ped$FatherID)])[fath_sexes != 1]
    if (length(fout) != 0) {
      msg <- sprintf(
        "Father ID %s has invalid sex.",
        paste0(unique(fout), collapse = ",")
      )
      rlang::abort(msg, "InvalidParentSex")
    }
  }
  
  # ----- CHECK TWINS -----
  # Twins that are NA get forced to 0
  if (sum(is.na(ped$Twins)) > 0) {
    msg <- sprintf("Column %s does not support NA. Forcing it to %s.", 
                   "Twins", 0)
    rlang::inform(msg, level = "InvalidNA")
    ped$Twins[is.na(ped$Twins)] <- 0
  }
  
  # Check for consistency among identical twins/multiple births
  if (any(ped$Twins != 0)) {
    ped <- .checkTwins(ped, unknown.race = unknown.race, 
                       unknown.ancestry = unknown.ancestry)
  }
  
  # ----- CHECK RISK MODIFIERS/INTERVENTIONS -----
  # Whether or not message about risk modifiers that are NA has been printed
  isPrintedRiskModNA = FALSE 
  
  # If both 'riskmod' and 'interAge' columns don't already exist
  if( sum(c("riskmod", "interAge") %in% colnames(ped)) == 0 ) {
          
          # If missing any of the 6 intervention columns
          if(sum(c("riskmodHyst", "interAgeHyst", "riskmodOoph", "interAgeOoph",
                   "riskmodMast", "interAgeMast") %in% colnames(ped)) != 6) {
                  rlang::inform("At least one intervention column is missing. Adding dummy columns.")
                  
                  if (!("riskmodHyst" %in% colnames(ped))) {
                          # Add riskmodHyst if not present
                          ped <- cbind(ped, data.frame(riskmodHyst = 0))
                  }
                  if (!("interAgeHyst" %in% colnames(ped))) {
                          # Add interAgeHyst if not present
                          ped <- cbind(ped, data.frame(interAgeHyst = NA))
                  }
                  if (!("riskmodOoph" %in% colnames(ped))) {
                          # Add riskmodOoph if not present
                          ped <- cbind(ped, data.frame(riskmodOoph = 0))
                  }
                  if (!("interAgeOoph" %in% colnames(ped))) {
                          # Add interAgeOoph if not present
                          ped <- cbind(ped, data.frame(interAgeOoph = NA))
                  }
                  if (!("riskmodMast" %in% colnames(ped))) {
                          # Add riskmodMast if not present
                          ped <- cbind(ped, data.frame(riskmodMast = 0))
                  }
                  if (!("interAgeMast" %in% colnames(ped))) {
                          # Add interAgeMast if not present
                          ped <- cbind(ped, data.frame(interAgeMast = NA))
                  }
          }
          
          # Non-binary riskmod (catches NAs too)
          if(sum( !(ped$riskmodHyst %in% c(0,1)) ) != 0 |
             sum( !(ped$riskmodOoph %in% c(0,1)) ) != 0 |
             sum( !(ped$riskmodMast %in% c(0,1)) ) != 0) {
                  
                  if(sum( !(ped$riskmodHyst %in% c(0,1)) ) != 0) {
                          rlang::abort("Entries in 'riskmodHyst' need to be binary (0 or 1).")
                  } 
                  if(sum( !(ped$riskmodOoph %in% c(0,1)) ) != 0) {
                          rlang::abort("Entries in 'riskmodOoph' need to be binary (0 or 1).") 
                  }
                  if(sum( !(ped$riskmodMast %in% c(0,1)) ) != 0) {
                          rlang::abort("Entries in 'riskmodMast' need to be binary (0 or 1).")
                  }
          }
          
          # 0's instead of NAs - interAge
          if(sum((ped$riskmodHyst == 0) & (ped$interAgeHyst == 0), na.rm=T) != 0 |  # if NA here, it's b/c age was NA
             sum((ped$riskmodOoph == 0) & (ped$interAgeOoph == 0), na.rm=T) != 0 |
             sum((ped$riskmodMast == 0) & (ped$interAgeMast == 0), na.rm=T) != 0) {
                  # Handled in condense() function
                  rlang::inform("Converting 0s in intervention age to NAs.")
          }
          
          # Strings instead of numbers - interAge
          if(sum(is.character(ped$interAgeHyst) & 
                 !is.na(ped$interAgeHyst) & 
                 suppressWarnings(is.na(as.numeric(ped$interAgeHyst)))) != 0 |
             sum(is.character(ped$interAgeOoph) & 
                 !is.na(ped$interAgeOoph) & 
                 suppressWarnings(is.na(as.numeric(ped$interAgeOoph)))) != 0 |
             sum(is.character(ped$interAgeMast) & 
                 !is.na(ped$interAgeMast) & 
                 suppressWarnings(is.na(as.numeric(ped$interAgeMast)))) != 0) {
                  rlang::inform("Converting intervention ages that are strings to NA.")
                  indexHyst <- which(is.character(ped$interAgeHyst) & 
                                             !is.na(ped$interAgeHyst) & 
                                             suppressWarnings(is.na(as.numeric(ped$interAgeHyst))))
                  indexOoph <- which(is.character(ped$interAgeOoph) & 
                                             !is.na(ped$interAgeOoph) & 
                                             suppressWarnings(is.na(as.numeric(ped$interAgeOoph))))
                  indexMast <- which(is.character(ped$interAgeMast) & 
                                             !is.na(ped$interAgeMast) & 
                                             suppressWarnings(is.na(as.numeric(ped$interAgeMast))))
                  if(length(indexHyst) > 0) {
                          ped[indexHyst, "interAgeHyst"] <- NA
                  }
                  if(length(indexOoph) > 0) {
                          ped[indexOoph, "interAgeOoph"] <- NA
                  }
                  if(length(indexMast) > 0) {
                          ped[indexMast, "interAgeMast"] <- NA
                  }
          }
          
          # If age but no intervention (also covers typos for riskmod = 0 rows)
          # Doesn't include if interAge = 0 (covered above)
          if(sum((ped$riskmodHyst == 0) & !is.na(ped$interAgeHyst) & (ped$interAgeHyst != 0)) != 0 |
             sum((ped$riskmodOoph == 0) & !is.na(ped$interAgeOoph) & (ped$interAgeOoph != 0)) != 0 |
             sum((ped$riskmodMast == 0) & !is.na(ped$interAgeMast) & (ped$interAgeMast != 0)) != 0) {
                  # Handled in condense() function 
                  rlang::inform("Removing intervention ages that are missing a corresponding intervention.")
          }
          
          # If intervention but no age
          if(sum((ped$riskmodHyst == 1) & is.na(ped$interAgeHyst)) != 0 |
             sum((ped$riskmodOoph == 1) & is.na(ped$interAgeOoph)) != 0 |
             sum((ped$riskmodMast == 1) & is.na(ped$interAgeMast)) != 0) {
                  # Handled in the "Riskmods that are NA will be removed" section 
                  rlang::inform("Removing interventions that are missing a corresponding age.")
          }
          
          # Call function to condense columns
          ped <- condense(ped) 
          
          # If any intervention age is below 18, throw a warning
          if(sum(unlist(ped$interAge) < 18, na.rm=T) != 0) {
                  rlang::warn("At least one intervention age is below 18.")
          } 
  }
  
  for (i in seq(nrow(ped))) { # Iterate through individuals in pedigree
    # Extract individual's risk modifiers and ages
    riskmod <- ped$riskmod[[i]]
    interv_age <- ped$interAge[[i]]
    
    # If there is only one risk modifier or age and it is NA, remove it
    if (length(riskmod) == 1 | length(interv_age) == 1) {
      riskmod_NA <- (length(riskmod) == 1) && is.na(riskmod)
      interv_age_NA <- (length(interv_age) == 1) && is.na(interv_age)
      if (riskmod_NA | interv_age_NA) {
        if (isPrintedRiskModNA == FALSE) { 
          rlang::inform("Riskmods that are NA will be removed", 
                        level = "InvalidInputForcingNone")
          isPrintedRiskModNA = TRUE 
        }
        riskmod <- ped$riskmod[[i]] <- "None"
        interv_age <- ped$interAge[[i]] <- 1
        next()
      }
    }
    
    if (length(riskmod) != length(interv_age)) {
      # If there is a mismatch in the number of risk modifiers and ages, 
      # truncate so that they match each other in length
      retain_n <- min(length(riskmod), length(interv_age))
      msg <- sprintf(
        "ID %s has %s risk modifications but has %s intervention age. Truncating to the first %s interventions.", 
        ped$ID[i], length(riskmod), length(interv_age), retain_n)
      rlang::inform(msg, level = "RiskmodInconsistency")
      riskmod <- riskmod[1:retain_n]
      interv_age <- interv_age[1:retain_n]
    } 
    
    # If zero ages are not allowed, set them to NA
    if (!allow.age.zero) {
      if (any(interv_age == 0, na.rm = TRUE)) {
        msg <- "Forcing age = 0 to be NA. "
        rlang::inform(msg, level = "FalseAgeZero")
        # Note that fix_with_rep also replaces negative values
        interv_age <- fix_with_rep(0, NA)(interv_age)
      }
    }
    
    # Round non-integer ages to integers and issue warning (all values < 1 
    # get rounded up to 1)
    baby_idx = which(interv_age > 0 & interv_age < 1)
    if (length(baby_idx) > 0) {
      rlang::inform(sprintf(
        "ID %s has intervention age between 0 and 1; rounding up to 1.",
        paste(ped$ID[i], collapse = ",")
      ), level = "AgeRounding")
      interv_age[baby_idx] = 1
    }
    
    dec_idx = which(as.numeric(interv_age) %% 1 != 0)
    if (length(dec_idx) > 0) {
      rlang::inform(sprintf(
        "ID %s has decimal-value intervention age; rounding to nearest integer.",
        paste(ped$ID[i], collapse = ",")
      ), level = "AgeRounding")
      interv_age[dec_idx] = round(interv_age[dec_idx])
    }
    
    # Force unsupported risk modifiers and unexpected ages to NA 
    riskmod <- forcingNA_ifnot_contains(ped, "Riskmod", 
                                        c("None", unlist(RISKMODS)), 
                                        vec = riskmod, force_NA = TRUE, 
                                        force_remove = FALSE)
    
    interv_age <- forcingNA_ifnot_contains(ped, "intervention age", 1:MAXAGE, 
                                           vec = interv_age, force_NA = TRUE, 
                                           force_remove = FALSE)
    
    # Remove entries where either the risk modifier or age is NA
    if (any(is.na(riskmod) | is.na(interv_age))) {
      ifretains <- !(is.na(riskmod) | is.na(interv_age))
      riskmod <- riskmod[ifretains]
      interv_age <- interv_age[ifretains]
      which_removed <- which(!ifretains)
      msg <- sprintf(
        "ID %s contains NA for either risk modifiers or intervention age, forcing removal of modification #%s.", 
        ped$ID[i], paste0(which_removed, collapse = ","))
      rlang::inform(msg, level = "InvalidInputForcingNone")
    }
    
    # Error if the number of risk modifiers and ages is different (redundant?)
    if (length(riskmod) != length(interv_age)) {
      rlang::abort(sprintf(
        "riskmod and interAge lengths not the same for individual with ID %s.", 
        ped$ID[i]), level = "RiskmodCheck")
    }
    
    # Error if there are any NAs remaining (redundant?)
    if (any(is.na(c(riskmod, interv_age)))) {
      rlang::abort(sprintf(
        "NAs not supported for riskmod and interAge in ID %s.", 
        ped$ID[i]), level = "RiskmodCheck")
    }
    
    # Error if there is a duplicated risk modifier
    if (any(duplicated(riskmod))) {
      msg <- sprintf(
        "ID %s contains duplicated risk modifier or intervention.", 
        ped$ID[i])
      rlang::abort(msg, level = "DuplicatedRiskmod")
    }
    
    # Error if any intervention is not preventative
    for (j in 1:length(RISKMODS)) {
      rm_idx <- which(riskmod == RISKMODS[[j]])
      rm_cancer = names(RISKMODS)[[j]]
      if (length(rm_idx) != 0 && 
          rm_cancer %in% affect_cancers && 
          ped[i,paste0("isAff",rm_cancer)] == 1 &&
          (is.na(ped[i,paste0("Age", rm_cancer)]) || 
           ped[i,paste0("Age", rm_cancer)] <= interv_age[rm_idx])) {
        msg <- sprintf(
          "ID %s had an intervention after the associated cancer. Interventions must be preventative.", 
          ped$ID[i])
        rlang::abort(msg, level = "InvalidRiskmod"
        )
      }
    }
    
    # If person has no risk modifiers/ages, set riskmod to "None" and age to 1
    if (length(riskmod) == 0 | length(interv_age) == 0) {
      riskmod <- ped$riskmod[[i]] <- "None"
      interv_age <- ped$interAge[[i]] <- 1
    } 
    
    # Replace riskmod and age in pedigree with the processed version
    ped$riskmod[[i]] <- riskmod
    ped$interAge[[i]] <- interv_age
  }
  
  # ----- CHECK TUMOR BIOMARKER TESTING -----
  # Iterate through all supported tumor markers
  for (j in 1:length(MARKER_TESTING)) { 
    # Assume that if a person has tumor marker testing results, they must be 
    # affected for the related cancer
    mt_cancer = names(MARKER_TESTING)[j]
    if(mt_cancer=="CBC"){ next } # someone with BC and without CBC can have ER result
    if (mt_cancer %in% .getCancersFromFam(ped)) {
      marker_subset <- intersect(colnames(ped), MARKER_TESTING[[j]][["MARKERS"]])
      if (length(marker_subset) >= 1) {
        should_have_been_diagnosed_cancer <- apply(ped[marker_subset], 1, 
                                                   function(v) any(!is.na(v)))
        if (any(should_have_been_diagnosed_cancer & 
                !ped[[paste0("isAff", mt_cancer)]])) {
          # If the assumption is violated, ignore the testing results
          to_fix_cancer <- ped$ID[should_have_been_diagnosed_cancer & 
                                    !ped[[paste0("isAff", mt_cancer)]]]
          rlang::inform(sprintf(
            "ID %s has tumor marker testing results but is unaffected for relevant cancer, so testing will be ignored.", 
            paste0(to_fix_cancer, collapse = ",")), 
            "InvalidTests")
          ped[ped$ID %in% to_fix_cancer, marker_subset] <- NA
        }
      }
    }
  }
  
  # ----- CHECK GERMLINE TESTING -----
  # Rename germline testing columns according to default variants and print 
  # message that testing results are assumed to be for the default variants
  for (gene in names(DEFAULT_VARIANTS)) {
    idx = which(names(ped) == gene) 
    if (length(idx) == 1) {
      names(ped)[idx] = DEFAULT_VARIANTS[gene]
      rlang::inform(
        sprintf(
          "Germline testing results for %s are assumed to be for default variant %s.", 
          gene, DEFAULT_VARIANTS[gene]), 
        "GermlineTestAssumption")
    }
  }
  
  # If requested, remove proband germline testing information by replacing 
  # results with NA. 
  if (ignore.proband.germ == TRUE) {
    ped[ped$ID %in% proband, intersect(colnames(ped), DEFAULT_VARIANTS)] = NA
  }
  
  
  # ----- TYPECHECKS -----
  # Error if there are duplicated IDs (redundant?)
  if (any(duplicated(ped$ID))) {
    rlang::abort("Duplicate IDs.", level = "TypeCheck")
  }
  # Error if pedigree does not contain all cancers in model specification (redundant?)
  if (!all(cancers %in% .getCancersFromFam(ped))) {
    rlang::abort(
      "Pedigree is missing cancers in model specification.", 
      level = "TypeCheck")
  }
  # Error if sex is missing for any individual
  if (!all(ped$Sex %in% c(0, 1))) {
    rlang::abort(
      "Pedigree contains missing sex for one or more individuals.", 
      level = "TypeCheck")
  }
  # Error if any individual has a missing twin status
  if (!all(is.finite(ped$Twins))) {
    rlang::abort(
      "Pedigree contains missing twin status for one or more individuals.", 
      level = "TypeCheck")
  }
  # Error if any individual has a missing cancer affection status
  if (length(affect_cancers) > 0) {
    aff_status <- ped[, paste0("isAff", affect_cancers)]
    if (!all(aff_status == 0 | aff_status == 1)) {
      rlang::abort(
        "Pedigree contains missing affection status for one or more individuals.", 
        level = "TypeCheck")
    }
  }
  
  # Identify the columns that are not lists
  non_list_cols <- colnames(ped)[sapply(ped, class) != "list"]
  # Replace NAs and negative values with -999
  ped[non_list_cols] <- lapply(ped[non_list_cols], fix_with_rep(NA, -999))
  
  # Identify the columns that are not lists and involve ages
  non_list_age_cols = non_list_cols[grep("Age", non_list_cols)]
  # Ensure that ages don't exceed MAXAGE, censoring if necessary
  isOlderThanMaxage = ped[non_list_age_cols] > MAXAGE
  if (sum(isOlderThanMaxage) > 0) {
    ped[non_list_age_cols][isOlderThanMaxage] = MAXAGE
    rlang::inform(
      sprintf(
        "Ages that are greater than the maximum age %s will be set to the maximum age.", 
        MAXAGE), 
      "ExceedsMAXAGE")
  }
  
  # ----- CHECK CONNECTIVITY -----
  # Ensure that race and ancestry are consistent
  ped <- .checkHeredity(ped)
  
  # Remove disconnected family members
  ped <- .removeDisconnected(ped, proband)
  
  # Ensure consistent gender information
  ped <- .checkGender(ped)
  # Map the sex variable from a numeric to character vector
  ped$Sex <- .mapGenderNames(number = ped$Sex)
  
  # Check or create age bound columns for each cancer
  if (length(cancers) > 0) {
    for (cc in cancers) {
      ped <- .boundAffectionAges(ped, cc)
    }
  }
  
  # Split disconnected pedigrees
  ped_list <- split(ped, ped$famID)
  # Get list of relationships for each pedigree
  rel_ls <- lapply(ped_list, get_relation_list)
  
  # Check each pedigree's mating structure
  mapply(function(ped, rel_l) .checkMating(ped, rel_l),
         ped = ped_list, rel_l = rel_ls
  )
  
  return(list(ped_list = ped_list, proband = proband, rel_ls = rel_ls))
}


#' Check and create age bound columns for cancer affection age
#'
#' Checks for the presence and validity of lower and upper age bound columns 
#' in `ped` for a given `cancer`, (e.g. `AgeBC_lower` and `AgeBC_upper` for 
#' breast cancer). If age bound columns are not in `ped`, they will be created. 
#' These bounds will be used for age imputation. 
#' 
#' @param ped A pedigree data frame. 
#' @param cancer A short cancer name corresponding to the cancer whose 
#' affection ages should be bounded. 
#' @return A modified pedigree data frame with valid lower and upper age bound 
#' columns for `cancer`. 
#' @family check
#' @examples
#' sample_ped1 <- data.frame(isAffBC = c(0, 1, 0), AgeBC = c(-999, -999, -999))
#' sample_ped2 <- data.frame(isAffBC = c(0, 1, 0), AgeBC = c(-999, 3, -999), 
#'                           AgeBC_lower = c(-999, 1, -999))
#' sample_ped3 <- data.frame(isAffBC = c(0, 1, 0), AgeBC = c(-999, 3, -999))
#' PPP:::.boundAffectionAges(sample_ped1, "BC")
#' PPP:::.boundAffectionAges(sample_ped2, "BC")
#' PPP:::.boundAffectionAges(sample_ped3, "BC")
.boundAffectionAges <- function(ped, cancer) {
  
  # Identify affection age column for cancer
  ped_cols <- colnames(ped)
  available_age_cols <- grep(paste0("Age", cancer), ped_cols, value = TRUE)
  
  # Add upper and lower age bound columns to pedigree, necessary
  if ((length(available_age_cols) == 1) && 
      (available_age_cols == paste0("Age", cancer))) {
    # If we only have the exact age column
    ped[[paste0("Age", cancer, "_lower")]] <- ped[[paste0("Age", cancer)]]
    ped[[paste0("Age", cancer, "_upper")]] <- ped[[paste0("Age", cancer)]]
  } else {
    # If we only have the exact age column and lower bound
    if ((length(available_age_cols) == 2) && 
        (paste0("Age", cancer, "_lower") %in% ped_cols)) {
      ped[[paste0("Age", cancer, "_upper")]] <- ped[[paste0("Age", cancer)]]
    }
    # If we only have the exact age column and upper bound
    if ((length(available_age_cols) == 2) && 
        (paste0("Age", cancer, "_upper") %in% ped_cols)) {
      ped[[paste0("Age", cancer, "_lower")]] <- ped[[paste0("Age", cancer)]]
    }
  }
  
  # Identify row indices of individuals who are affected for cancer and have a 
  # non-missing cancer age
  ids2fix <- which((ped[[paste0("isAff", cancer)]] == 1) & 
                     (ped[[paste0("Age", cancer)]] != -999))
  for (row_num in ids2fix) {
    # Ensure that the lower bound is smaller than the upper bound
    if ((ped[row_num, paste0("Age", cancer, "_lower")] > 
         ped[row_num, paste0("Age", cancer, "_upper")]) & 
        (ped[row_num, paste0("Age", cancer, "_upper")] != -999)) {
      rlang::abort(
        "Affection age lower bound should be lower than upper bound", 
        level = "AgeBounds")
    }
    
    # Assign lower bound as the maximum of the lower bound and exact age
    ped[row_num, paste0("Age", cancer, "_lower")] <- 
      max(ped[row_num, paste0("Age", cancer, "_lower")], 
          ped[row_num, paste0("Age", cancer)])
    # Assign upper bound as the minimum of the lower bound and exact age
    ped[row_num, paste0("Age", cancer, "_upper")] <- 
      min(ped[row_num, paste0("Age", cancer, "_upper")], 
          ped[row_num, paste0("Age", cancer)])
  }
  
  return(ped)
}


#' Check for consistency among homozygotic twins
#'
#' Identical twins and other multiple births encoded by the `ped$Twin` column 
#' are expected to have the same parents, sexes, current ages (if all are 
#' alive), races, and ancestries. 
#' * Twin sets that include only one individual will be removed (the individual 
#' will remain in the pedigree, but will not be treated as a twin). 
#' * If `MotherID`, `FatherID`, `Sex`, or `CurAge` is known for one or more 
#' twins in the set but is missing for others, the missing value(s) will be 
#' inferred from the other twin(s). 
#' * Otherwise, if there is a conflict for `MotherID`, `FatherID`, or `Sex`, 
#' the twin label will be removed for the set. 
#' * Among twins who are currently alive, the proband's current age will take 
#' precedence in the case of `CurAge` conflicts if one (or more) of the twins 
#' is a proband. If none of the twins are probands, the average of the current 
#' ages will be assigned as `CurAge` for all members. 
#' * If the current age of any of the living twins is less than `CurAge` of a 
#' dead twin, assign the average as `CurAge` for all members. 
#' * If there is a conflict in `race` or `Ancestry`, all members of the twin 
#' set will be assigned the default values specified by the `unknown.race` and 
#' `unknown.ancestry` arguments. 
#'
#' @param ped A pedigree data frame. 
#' @param unknown.race A character string indicating the default race to use 
#' when race information is missing. The default is `PPP:::UNKNOWN_RACE`. 
#' @param unknown.ancestry A character string indicating the default ancestry 
#' to use when ancestry information is missing. The default is 
#' `PPP:::UNKNOWN_ANCESTRY`. 
#' @family check
#' @return A modified pedigree data frame with problematic twin sets resolved 
#' or removed. 
.checkTwins <- function(ped, unknown.race = UNKNOWN_RACE, 
                        unknown.ancestry = UNKNOWN_ANCESTRY) {
  
  # Default race and ancestry
  default <- c(race = unknown.race, ethnic = unknown.ancestry)
  
  # Abort if there pedigree has no twins (those with twin labels other than 0) 
  if (all(ped$Twins == 0)) {
    rlang::abort(
      "checkTwins should not be run when there are no twins in the pedigree.", 
      level="NoTwins")
  }
  # Identify (nonzero) twin labels
  twin_labels <- unique(ped$Twins[ped$Twins != 0])
  
  # Columns that must be identical but that do not have default values
  set_1 <- c("FatherID", "MotherID", "CurAge", "Sex")
  # Columns that will be changed to default values if there is a conflict
  set_2 <- c("race", "Ancestry")
  
  # Iterate through twin sets
  for (label in twin_labels) {
    # Number of twins in the set
    twin_num <- sum(ped$Twins == label)
    # Skip the twin set if there is only one twin in the set
    if (twin_num == 1) {
      rlang::inform(sprintf(
        "Twins labelled %s only have one entry, forcing it to 0", 
        label
      ), level = "MissingTwinPair")
      ped$Twins[ped$Twins == label] <- 0
      next()
    }
    
    # Identify IDs of twins in the set
    twin_ID_set <- ped$ID[ped$Twins == label]
    
    # Identify variables in set_1 where the twins disagree
    var_1 <- ped[ped$Twins == label, set_1]
    if_agree <- apply(var_1, 2, function(v) length(unique(v)) == 1)
    disagree_var <- set_1[!if_agree]
    
    # Iterate through the set_1 variables twin disagreements
    for (dis in disagree_var) {
      dis_values <- var_1[dis]
      
      if (-999 %in% dis_values) {
        # If the disagreement is because one of the twins has a missing value, 
        # infer its value from the other twin(s)
        informative_value <- unique(dis_values[dis_values != -999])[1]
        who_is_missing <- ped$ID[ped[[dis]] == -999]
        rlang::inform(sprintf(
          "One of the twins labelled %s has missing %s that can be inferred from the other twin.", 
          label, dis), level = "MissingTwinInfoModified")
        ped[ped$ID == who_is_missing, dis] <- informative_value
      } else {
        # If the conflict is a parent ID or sex, remove the twin label
        if (dis == "MotherID" | dis == "FatherID" | dis == "Sex") {
          rlang::warn(sprintf(
            "%s conflict, twin labelled %s is removed.",
            dis, label
          ), level = "TwinInfoConflict")
          ped$Twins[ped$Twins == label] <- 0
        }
        
        # If the conflict is current age, there are several viable scenarios
        if (dis == "CurAge") {
          # Identify twins who are alive in this set
          alive_twin_ID_set = ped$ID[ped$Twins == label & ped$isDead == 0]
          # Identify twins who are dead in this set
          dead_twin_ID_set = ped$ID[ped$Twins == label & ped$isDead == 1]
          
          # If at least two twins are currently alive, but have different ages
          if (length(alive_twin_ID_set) >= 2) {
            # Twins who are alive and probands in this set
            proband_alive_twin_ID_set = ped$ID[ped$Twins == label & 
                                                 ped$isDead == 0 & 
                                                 ped$isProband == 1]
            if (length(proband_alive_twin_ID_set) == 1) {
              # If one of the twins is alive and a proband, set current ages of
              # all living twins to the proband twin's age
              ped$CurAge[ped$ID %in% alive_twin_ID_set] = 
                ped$CurAge[proband_alive_twin_ID_set == ped$ID]
              rlang::warn(sprintf(
                "CurAge conflict for twin label %s; assigning CurAge of twin who is a proband to all living twins.", 
                label), level = "TwinInfoConflict")
            } else {
              # Otherwise, assign the average of the current ages to all living 
              # twins
              ped$CurAge[ped$ID %in% alive_twin_ID_set] = round(mean(ped$CurAge[ped$ID %in% alive_twin_ID_set]))
              rlang::warn(sprintf(
                "CurAge conflict for twin label %s; assigning CurAge for all living twins to be the average of their ages.", 
                label), level = "TwinInfoConflict")
            }
          } 
          
          # If the current age of the living twin(s) is less than the death age 
          # of a dead twin, assign the average of the current ages to all twins
          if (length(alive_twin_ID_set) > 0 && length(dead_twin_ID_set) > 0 &&
              max(ped$CurAge[ped$ID %in% alive_twin_ID_set]) < 
              min(ped$CurAge[ped$ID %in% dead_twin_ID_set])) {
            ped$CurAge[ped$ID %in% twin_ID_set] = 
              round(mean(ped$CurAge[ped$ID %in% twin_ID_set]))
            rlang::warn(sprintf(
              "Living twin with twin label %s is currently younger than twin who died;
              assigning CurAge to be the avreage of their ages", 
              label), level = "TwinInfoConflict")
          }
        }
      }
    }
    
    # Identify variables in set_2 where the twins disagree
    var_2 <- ped[ped$Twins == label, set_2]
    if_agree <- apply(var_2, 2, function(v) length(unique(v)) == 1)
    disagree_var <- set_2[!if_agree]
    
    # Iterate through the set_2 variables twin disagreements
    for (dis in disagree_var) {
      # Assign default values to all twins
      rlang::warn(sprintf(
        "Twins are assumed monozygotic but %s conflict for twin label %s. Forcing them to default values.", 
        dis, label), level = "TwinInfoConflict")
      ped[ped$Twins == label, dis] <- default[dis]
    }
  }
  
  return(ped)
}


#' Remove disconnected family members
#'
#' @param ped A pedigree data frame. 
#' @param proband A numeric value or vector of the unique IDs in `ped` 
#' corresponding to the proband(s). 
#' @return A modified pedigree data frame that only includes members connected 
#' to the proband(s) and has an additional `famID` column for family groupings.  
#' @family check
.removeDisconnected <- function(ped, proband) {
  # Initialize the new famID variable as a column of 1s
  ped$famID <- rep(1, nrow(ped))
  
  # Useing kinship2's makefamid function, construct family groupings from the 
  # pedigree information
  fvec <- with(ped, {
    FatherID <- as.character(FatherID)
    MotherID <- as.character(MotherID)
    FatherID[is.na(FatherID)] <- ""
    MotherID[is.na(MotherID)] <- ""
    kinship2::makefamid(ID, FatherID, MotherID)
  })
  
  # Identify the family group(s) that contains proband(s)
  proband_fam_id <- fvec[ped$ID %in% proband]
  # find out singleton proband(s)
  newped1 <- ped[which(fvec==0 & ped$isProband==1), ]
  # Subset the pedigree to remove all members who are not part of the subfamily
  newped2 <- ped[which(fvec %in% proband_fam_id[which(proband_fam_id!=0)]), ]
  # rbind two subsets
  newped <- rbind(newped1, newped2)
  
  # Store the IDs in the new and original pedigrees
  newped_ids <- newped$ID
  oriped_ids <- ped$ID
  
  # Print warning message that disconnected relatives have been removed
  if (!identical(newped_ids, oriped_ids)) {
    rlang::warn(sprintf(
      "ID %s has been removed from pedigree since they are not connected with any proband.", 
      paste0(setdiff(oriped_ids, newped_ids), collapse = ",")), 
      level = "DisconectedFamilyRemoval")
  }
  
  return(newped)
}


#' Correct inconsistencies between individual sex and cancer
#' 
#' @param ped A pedigree data frame. 
#' @param fcancers A vector of cancer short names that only apply to females. 
#' @param mcancers A vector of cancer short names that only apply to males. 
#' @param gender A character string representing the sex to fix, either 
#' `"male"` or `"female"`. The default is `"male"`. 
#' @return A modified pedigree data frame with inconsistencies between sex and 
#' cancer corrected, such that individuals who are affected for cancers that 
#' they cannot have (based on their sex) have their cancer status changed to 
#' non-affected. 
#' @family check
.fix_gender_cancer <- function(ped, fcancers, mcancers, gender = "male") {
  
  # Identify the sex we want to fix and the corresponding cancers
  gender <- match.arg(gender, c("male", "female"))
  specific_cancers <- ifelse(gender == "male", fcancers, mcancers)
  specific_cols <- paste0("isAff", specific_cancers)
  specific_sex <- ifelse(gender == "male", 1, 0)
  
  # Find row indices of problematic individuals
  row_id <- which(apply(ped[specific_cols], 1, function(v) 1 %in% v) & 
                    ped$Sex == specific_sex)
  # For each problematic individual, identify the specific cancers that they 
  # cannot have (based on their sex), turn the cancer status to non-affected, 
  # and print a warning message
  for (rid in row_id) {
    id <- ped$ID[rid]
    wrong_cancers <- specific_cancers[which(ped[specific_cols][rid, ] == 1)]
    rlang::warn(sprintf(
      "ID %s is %s but has been affected by %s; turning into non-affection.",
      id, gender, paste0(wrong_cancers, collapse = ",")
    ), level = "GenderCancerMismatch")
    ped[rid, paste0("isAff", wrong_cancers)] <- 0
  }
  
  return(ped)
}


#' Ensure consistent gender information
#'
#' First, uses the `fixParents` function from the `kinship2` package to fix the 
#' sexes of parents and add IDs of missing parents. Then, checks that cancer 
#' affection status is compatible with the sex of the individual. When 
#' individuals are affected for cancers that they cannot have (based on their 
#' sex), their cancer status is changed to non-affected and a warning message 
#' is printed. 
#' 
#' @param ped A pedigree data frame. 
#' @param gender.specific.cancers A list with two components: `Male`, which 
#' consists of a vector of short cancer names specific to males; and `Female`, 
#' which consists of a vector of short cancer names specific to females. The 
#' default is `list(Male = c("PROS"), Female = c("CER", "ENDO", "OC")`. 
#' @return A modified pedigree data frame with consistent gender information. 
#' @family check
.checkGender <- function(
  ped, 
  gender.specific.cancers = list(Male = c("PROS"), 
                                 Female = c("CER", "ENDO", "OC"))
  ) {
  
  # Save a copy of the pedigree data frame as the metadata
  metadata <- ped
  # Using kinship2's fixParents function, create a pedigree data frame that 
  # fixes the sexes of parents and adds IDs of missing parents
  fixdf <- with(ped, {
    Sex[Sex == 0] <- 2
    fdf <- kinship2::fixParents(ID, FatherID, MotherID, Sex, missid = -999)
    # Apparently, fixParents uses 0 instead of missid when adding parents with 
    # missing IDs
    fdf$momid[fdf$momid == 0] <- -999 
    fdf$dadid[fdf$dadid == 0] <- -999
    fdf$sex[fdf$sex == 2] <- 0
    colnames(fdf) <- c("ID", "MotherID", "FatherID", "Sex")
    fdf
  })
  
  # Update the pedigree by merging in the fixed genders
  ped <- merge(fixdf[c("ID", "Sex")], 
               metadata[!names(metadata) == "Sex"], all.x = TRUE)
  # Maintain the same order of individuals
  ped <- ped[ped$ID %in% metadata$ID, ]
  # Check that the IDs were not duplicated in the process
  if (any(duplicated(ped$ID))) {
    rlang::abort("Duplicate IDs.", level = "TypeCheck")
  }
  
  # Replace NAs with -999 for all non-list columns
  non_list_cols <- colnames(ped)[sapply(ped, class) != "list"]
  ped[non_list_cols] <- lapply(ped[non_list_cols], fix_with_rep(NA, -999))
  
  # Cancers in pedigree
  cancer_domain <- .getCancersFromFam(ped)
  # Only males can have these cancers
  mcancers <- intersect(cancer_domain, gender.specific.cancers$Male)
  # Only females can have these cancers
  fcancers <- intersect(cancer_domain, gender.specific.cancers$Female)
  
  # Check if there are inconsistencies between individual genders and the 
  # cancers they can get
  if (length(mcancers) > 0) {
    ped <- .fix_gender_cancer(ped, fcancers, mcancers, "female")
  }
  if (length(fcancers) > 0) {
    ped <- .fix_gender_cancer(ped, fcancers, mcancers, "male")
  }
  
  # Error if sex is missing for anyone in the pedigree
  if (any(ped$Sex == -999)) {
    rlang::abort(sprintf(
      "ID %s are missing Sex",
      paste0(ped$ID[ped$Sex == -999], collapse = ",")
    ),
    level = "InvalidNA"
    )
  }
  
  return(ped)
}


#' Ensure consistent hereditary features within trios
#'
#' Checks that race and ancestry are consistent within trios (a child and 
#' their two parents). If they are not consistent, the information will be 
#' modified accordingly and a descriptive message will be issued. 
#' 
#' @details 
#' In particular, if both parents share a feature but their child does not, the 
#' child's feature will be corrected to match their parents. If all children of 
#' a set of parents share a feature, but one of the parents does not, the 
#' the parent's feature will be corrected to match their children. 
#'
#' @param ped A pedigree data frame. 
#' @return A modified pedigree data frame with consistent hereditary features. 
#' @family check
.checkHeredity <- function(ped) {
  # Subset pedigree for the ID, MotherID, and FatherID columns
  # Only include individuals with at least one parent in the pedigree
  trio_df <- unique(ped[ped$MotherID != -999 | ped$FatherID != -999, 
                        c("ID", "MotherID", "FatherID")])
  
  if (nrow(trio_df) > 0) {
    # Create new data frame that has a row corresponding to each mother and 
    # father ID combination
    temp_df <- stats::aggregate(trio_df, 
                         by = list(trio_df$FatherID, trio_df$MotherID), 
                         FUN = c, simplify = FALSE)
    
    if (nrow(temp_df) > 0) {
      for (i in 1:nrow(temp_df)) { # Iterate through aggregated trio df
        
        # Not entirely sure why this is needed, but ensure that there are no 
        # duplicate parent IDs in each row
        temp_df$MotherID <- sapply(temp_df$MotherID, unique)
        temp_df$FatherID <- sapply(temp_df$FatherID, unique)
        
        # Identify parents in this row
        parents <- unlist(ped$ID %in% temp_df[i, c("MotherID", "FatherID")])
        # Identify children of the parents
        kids <- temp_df[i, "ID"][[1]]
        
        for (feature in c("race", "Ancestry")) { # Check hereditary features
          # Extract features from parents and children
          children_feature <- ped[[feature]][ped$ID %in% kids]
          parents_feature <- ped[[feature]][ped$ID %in% parents]
          
          if (length(unique(parents_feature) == 1) && 
              any(children_feature != parents_feature)) {
            # If the parents share a common feature, but a child's feature  
            # differs from it, correct the child's feature using their parents' 
            # feature
            who_is_diff <- kids[children_feature != unique(parents_feature)]
            rlang::inform(
              sprintf(
                "ID %s 's %s has been changed to %s to meet hereditary consistency.", 
                paste0(who_is_diff, collapse = ","), 
                feature, unique(parents_feature)), 
              level = "HeredityFix")
            ped[[feature]][ped$ID %in% who_is_diff] <- unique(parents_feature)
          } else if (length(unique(children_feature)) == 1 && 
                     any(parents_feature != children_feature)) {
            # If the children share a common feature, but a parent's feature 
            # differs from it, correct the parent's feature using their 
            # children's feature
            who_is_diff <- parents[parents_feature != unique(children_feature)]
            rlang::inform(
              sprintf(
                "ID %s 's %s has been changed to %s to meet hereditary consistency.", 
                paste0(who_is_diff, collapse = ","), 
                feature, unique(children_feature)), 
              level = "HeredityFix")
            ped[[feature]][ped$ID %in% who_is_diff] <- unique(children_feature)
            
          }
        }
      }
    }
  }
  return(ped)
}


#' Find IDs of an individual's close relatives
#' 
#' Gets IDs in `ped` corresponding to the individual and their mother, father, 
#' parents' other spouses, full siblings, spouses, and children, if any. 
#' 
#' @param ped A pedigree data frame. 
#' @param indid An individual ID. 
#' @return A list with the following 8 numeric-valued components: 
#' * `indid`: The individual's ID (the same as `indid`). 
#' * `mid`: The ID of the individual's mother. 
#' * `fid`: The ID of the individual's father.  
#' * `smid`: The IDs of the mother's other spouses (besides `fid`). 
#' * `sfid`: The IDs of the father's other spouses (besides `mid`). 
#' * `fsibid`: The IDs of the individual's full siblings. 
#' * `sid`: The IDs of the individual's spouses. 
#' * `cid`: The IDs of the individual's children. 
#' @family check
#' @seealso \code{\link{.getOtherRels}}, \code{\link{.getSpouses}}
.findIDs <- function(ped, indid) {
  # Initialize empty vectors for each ID type
  mid <- fid <- smid <- sfid <- fsibid <- cid <- sid <- numeric(0)
  
  # Identify the individual's mother and father IDs from the pedigree columns
  mid <- ped[ped$ID == indid, "MotherID"]
  fid <- ped[ped$ID == indid, "FatherID"]
  
  # Identify IDs of the parents' other spouses and the individual's full siblings
  if (mid != -999 | fid != -999) { 
    # If at least one of the parents' IDs is known
    if (mid != -999 & fid != -999) { 
      # If both parents' IDs are known
      smid <- ped[ped$MotherID == mid & ped$FatherID != fid, "FatherID"]
      sfid <- ped[ped$FatherID == fid & ped$MotherID != mid, "MotherID"]
    } else {
      if (mid != -999) { 
        # If only the mother's ID is known
        smid <- ped[ped$MotherID == mid & ped$FatherID != fid, "FatherID"]
      } else { 
        # If only the father's ID is known
        sfid <- ped[ped$FatherID == fid & ped$MotherID != mid, "MotherID"]
      }
    }
    
    fsibid <- ped[ped$FatherID == fid & ped$MotherID == mid &
                    ped$ID != indid, "ID"]
  } 
  smid <- unique(smid[smid != -999])
  sfid <- unique(sfid[sfid != -999])
  
  # Identify IDs of the individual's spouses and children
  if (ped[ped$ID == indid, "Sex"] == "Male") { # Male individual
    temp <- ped$FatherID == indid
    cid <- ped[temp, "ID"]
    sid <- ped[temp, "MotherID"]
    sid <- unique(sid[sid != -999])
  } else { # Female individual
    temp <- ped$MotherID == indid
    cid <- ped[temp, "ID"]
    sid <- ped[temp, "FatherID"]
    sid <- unique(sid[sid != -999])
  }
  
  # If the mother's ID is missing (-999), set mid to be an empty vector
  if (mid == -999) {
    mid <- numeric(0)
  }
  # If the father's ID is missing (-999), set fid to be an empty vector
  if (fid == -999) {
    fid <- numeric(0)
  }
  
  # Put all the ID vectors together into a list
  ids <- list(indid, mid, fid, smid, sfid, fsibid, sid, cid)
  names(ids) <- c("indid", "mid", "fid", "smid", "sfid", "fsibid", "sid", "cid")
  return(ids)
}


#' Find IDs of an individual's distant relatives
#' 
#' Gets IDs in `ped` corresponding to the individual's cousins, their mother's 
#' full and half siblings, and their father's full and half siblings, if any. 
#' 
#' @param ped A pedigree data frame. 
#' @param indid An individual ID. 
#' @return A list with the following 3 numeric-valued components: 
#' * `cousin`: The IDs of the individual's cousins. 
#' * `msibs`: The IDs of the individual's mother's full and half siblings.   
#' * `psibs`: The IDs of the individual's father's full and half siblings. 
#' @family check
#' @seealso \code{\link{.findIDs}}, \code{\link{.getSpouses}}
.getOtherRels <- function(ped, indid) {
  # Initialize empty vectors for each ID type
  cid <- msibs <- psibs <- NULL
  
  # Identify the individual's mother and father IDs from the pedigree columns
  mid <- ped[ped$ID == indid, "MotherID"]
  fid <- ped[ped$ID == indid, "FatherID"]
  
  # If the mother's ID is known
  if (mid != -999) { 
    # Identify IDs of the mother's parents
    mgmid <- ped[ped$ID == mid, "MotherID"]
    mgfid <- ped[ped$ID == mid, "FatherID"]
    
    # Identify IDs of the mother's full and half siblings
    # Maternal grandparents' IDs are both known
    if (mgmid != -999 & mgfid != -999) { 
      msibs <- ped[(ped$FatherID == mgfid | 
                      ped$MotherID == mgmid) & ped$ID != mid, "ID"]
    } 
    # Only maternal grandmother's ID is known
    if (mgmid != -999 & mgfid == -999) { 
      msibs <- ped[ped$MotherID == mgmid & ped$ID != mid, "ID"]
    } 
    # Only maternal grandfather's ID is known
    if (mgmid == -999 & mgfid != -999) { 
      msibs <- ped[ped$FatherID == mgfid & ped$ID != mid, "ID"]
    } 
    
    # If the mother has siblings, identify IDs of the individual's maternal 
    # cousins
    if (length(msibs) != 0) {
      for (jjj in 1:length(msibs)) {
        if (ped[ped$ID == msibs[jjj], "Sex"] == "Male") { # Sibling is male
          temp <- ped$FatherID == msibs[jjj]
          cid <- c(cid, ped[temp, "ID"])
        } else { # Sibling is female
          temp <- ped$MotherID == msibs[jjj]
          cid <- c(cid, ped[temp, "ID"])
        }
      }
    }
  }
  
  # If the father's ID is known
  if (fid != -999) {
    # Identify IDs of the father's parents
    pgmid <- ped[ped$ID == fid, "MotherID"]
    pgfid <- ped[ped$ID == fid, "FatherID"]
    
    # Identify IDs of the father's full and half siblings
    # Paternal grandparents' IDs are both known
    if (pgmid != -999 | pgfid != -999) {
      psibs <- ped[(ped$FatherID == pgfid | ped$MotherID == pgmid) & 
                     ped$ID != fid, "ID"]
    } 
    # Only paternal grandmother's ID is known
    if (pgmid != -999 & pgfid == -999) {
      psibs <- ped[ped$MotherID == pgmid & ped$ID != fid, "ID"]
    } 
    # Only paternal grandfather's ID is known
    if (pgmid == -999 & pgfid != -999) {
      psibs <- ped[ped$FatherID == pgfid & ped$ID != fid, "ID"]
    } 
    
    # If the father has siblings, identify IDs of the individual's paternal 
    # cousins
    if (length(psibs) != 0) {
      for (jjj in 1:length(psibs)) {
        if (ped[ped$ID == psibs[jjj], "Sex"] == "Male") { # Sibling is male
          temp <- ped$FatherID == psibs[jjj]
          cid <- c(cid, ped[temp, "ID"])
        } else { # Sibling is female
          temp <- ped$MotherID == psibs[jjj]
          cid <- c(cid, ped[temp, "ID"])
        }
      }
    }
  }
  
  # Put all the ID vectors together into a list
  out <- list(cid, msibs, psibs)
  names(out) <- c("cousid", "msibs", "psibs")
  return(out)
}


#' Find IDs of an individual's spouses
#' 
#' Gets IDs in `ped` corresponding to the individual's spouses, if any. 
#' 
#' @param ped A pedigree data frame. 
#' @param indid An individual ID. 
#' @return A vector of the IDs of the individual's spouses. 
#' @family check
#' @seealso \code{\link{.getOtherRels}}, \code{\link{.findIDs}}
.getSpouses = function(ped, indid) {
  if (length(indid) == 0) { # If the ID does not exist, return nothing
    sid <- numeric(0)
  } else { # Otherwise, identify IDs of the individual's spouses 
    if (ped[ped$ID == indid, "Sex"] == "Male") { # Male individual
      temp <- ped$FatherID == indid
      sid <- ped[temp, "MotherID"]
      sid <- unique(sid[sid != -999])
    } else { # Female individual
      temp <- ped$MotherID == indid
      sid <- ped[temp, "FatherID"]
      sid <- unique(sid[sid != -999])
    }
  }
  
  return(sid)
}


#' Check a pedigree's mating structure for validity
#'
#' Raises errors for invalid mating patterns, i.e. intermarriages that create 
#' "loops" in the family: 
#' * Mating between siblings or first cousins. 
#' * Mating between a child's paternal aunt/uncle and maternal aunt/uncle. 
#' * Mating between children and their parents, aunts, or uncles. 
#' * Mating between a person's relatives and the relatives of any of their 
#' spouses. 
#'
#' @param fff A pedigree data frame. 
#' @param rel_l The list of relationships in `ped` returned by a call to
#' \code{\link{get_relation_list}}. 
#' @family check
#' @return None
.checkMating <- function(fff, rel_l) {
  # Number of relatives
  psize <- nrow(fff)
  # Pedigree with only the ID, Sex, FatherID, and MotherID columns
  ped <- fff[, c("ID", "Sex", "FatherID", "MotherID")]
  
  err <- "BadMating"
  
  # Iterate through each individual in the pedigree
  for (pp in 1:psize) {
    # Extract relationship information for the individual
    tempid <- rel_l[[pp]] 
    
    # Error if the individual mated with their sibling
    if (length(intersect(tempid$sid, tempid$fsibid)) > 0) {
      if (any(intersect(tempid$sid, tempid$fsibid) != -999)) {
        rlang::abort(sprintf(
          "ID %s mated with their %s.",
          ped$ID[pp], "sibling"
        ), level = err)
      }
    }
    
    # Error if the individual mated with their first cousin
    if (length(intersect(tempid$sid, tempid$cousid)) > 0) {
      if (any(intersect(tempid$sid, tempid$cousid) != -999)) {
        rlang::abort(sprintf(
          "ID %s mated with their %s.",
          ped$ID[pp], "cousin"
        ), level = err)
      }
    }
    
    # Error if one of the individual's paternal aunts/uncles mated with a 
    # maternal aunt/uncle
    temp <- table(tempid$cousid)
    if (length(temp) != 0) {
      if (sum(temp) / length(temp) > 1) {
        rlang::abort(sprintf(
          "ID %s 's paternal aunt/uncle mated with their maternal aunt or uncle.",
          ped$ID[pp]), level = err)
      }
    }
    
    # Error if the individual mated with their mother
    if (length(intersect(tempid$sid, tempid$mid)) > 0) {
      if (any(intersect(tempid$sid, tempid$mid) != -999)) {
        rlang::abort(sprintf(
          "ID %s mated with their %s.",
          ped$ID[pp], "mother"
        ), level = err)
      }
    }
    
    # Error if the individual mated with their father
    if (length(intersect(tempid$sid, tempid$fid)) > 0) {
      if (any(intersect(tempid$sid, tempid$fid) != -999)) {
        rlang::abort(sprintf(
          "ID %s mated with his/their %s.",
          ped$ID[pp], "father"
        ), level = err)
      }
    }
    
    # Error if the individual mated with their maternal aunt or uncle
    if (length(intersect(tempid$sid, tempid$msibs)) > 0) {
      if (any(intersect(tempid$sid, tempid$msibs) != -999)) {
        rlang::abort(sprintf(
          "ID %s mated with their %s.",
          ped$ID[pp], "maternal aunt or uncle"
        ), level = err)
      }
    }
    
    # Error if the individual mated with their paternal aunt or uncle
    if (length(intersect(tempid$sid, tempid$psibs)) > 0) {
      if (any(intersect(tempid$sid, tempid$psibs) != -999)) {
        rlang::abort(sprintf(
          "ID %s mated with their %s.",
          ped$ID[pp], "paternal aunt or uncle"
        ), level = err)
      }
    }
    
    # Error if the individual's relatives are mated to any of their mate(s) 
    # relatives
    for (sid in tempid$sid) { # Iterate through individual's mates
      # Identify IDs for all of the individual's relatives
      # Includes the individual and their children, but not their mate
      pp_rels <- c(
        ped$ID[pp], tempid$mid, tempid$fid, tempid$smid,
        tempid$sfid, tempid$fsibid, tempid$cid,
        tempid$cousin, tempid$msibs, tempid$psibs, tempid$extra_spouses
      )
      
      # Extract relationship information for this mate
      tempsid <- rel_l[[which(ped$ID == sid)]]
      # Identify IDs for all of the mate's relatives
      # Does not include the individual, the mate, or their children, but 
      # includes all of the mate's other mates and children
      sid_rels <- c(
        tempsid$mid, tempsid$fid, tempsid$smid, 
        tempsid$sfid, tempsid$fsibid, 
        tempsid$cid[!(tempsid$cid %in% tempid$cid)],
        tempsid$sid[tempsid$sid != ped$ID[pp]],
        tempsid$cousin, tempsid$msibs, tempsid$psibs, tempsid$extra_spouses
      )
      
      # Check that there is no overlap between the two sets of relatives
      if (length(intersect(pp_rels, sid_rels)) > 0) {
        if (any(intersect(pp_rels, sid_rels) != -999)) {
          rlang::abort(sprintf(
            "IDs %s and %s are mated but have shared relative(s) with ID(s) %s",
            ped$ID[pp], sid, paste0(intersect(pp_rels, sid_rels), collapse = ",")
          ), level = err)
        }
      }
    }
  }
}


#' Get list of relationships from a pedigree
#' 
#' @param ped A pedigree data frame. 
#' @return The list of relationships in `ped`. Each list element corresponds to 
#' an ID in `ped` and contains a vector of relative IDs returned by 
#' \code{\link{.findIDs}} and \code{\link{.getOtherRels}}. Also includes 
#' IDs of additional spouses, returned by apply \code{\link{.getSpouses}} to 
#' various IDs returned by \code{\link{.findIDs}} and \code{\link{.getOtherRels}}. 
get_relation_list <- function(ped) {
  rel_l <- list()
  for (pp in 1:nrow(ped)) {
    rel_l[[pp]] <- c(
      .findIDs(ped, ped$ID[pp]), # IDs of close relatives
      .getOtherRels(ped, ped$ID[pp]) # IDs of distant relatives
    )
    # Spouses of mother's spouses, father's spouses, siblings, 
    # cousins, mother's siblings, and father's siblings
    rel_l[[pp]]$extra_spouses = unlist(sapply(c(rel_l[[pp]]$smid, 
                                                rel_l[[pp]]$sfid, 
                                                rel_l[[pp]]$fsibid, 
                                                rel_l[[pp]]$cousin, 
                                                rel_l[[pp]]$msibs, 
                                                rel_l[[pp]]$psibs), 
                                              .getSpouses, ped = ped))
    # Remove individual's ID from the list of extra spouses
    rel_l[[pp]]$extra_spouses = rel_l[[pp]]$extra_spouses[rel_l[[pp]]$extra_spouses 
                                                             != ped$ID[pp]]
  }
  return(rel_l)
}


#' Condenses three riskmod and three interAge columns into one riskmod and
#' one interAge column.
#' 
#' @param ped A pedigree data frame.
#' @return A pedigree data frame with two intervention columns instead of six
condense <- function(ped) {
        
        # Create new columns
        ped$riskmod <- as.list(NA)
        ped$interAge <- as.list(NA)
        
        # Populate new columns w/ info from the others
        for(i in 1:nrow(ped)) {
                
                row <- ped[i, ]
                
                # Start w/ riskmod
                hyst <- row$riskmodHyst
                ooph <- row$riskmodOoph
                mast <- row$riskmodMast
                
                # Condition 1 (Y, N, N)
                if(hyst == 1 &
                   ooph == 0 &
                   mast == 0) {
                        row$riskmod <- list("Hysterectomy")
                        row$interAge <- list(as.numeric(row$interAgeHyst))
                }
                # Condition 2 (Y, Y, N)
                if(hyst == 1 &
                   ooph == 1 &
                   mast == 0) {
                        row$riskmod <- list(c("Hysterectomy", "Oophorectomy"))
                        row$interAge <- list(as.numeric(c(row$interAgeHyst,
                                                          row$interAgeOoph)))
                }
                # Condition 3 (Y, N, Y)
                if(hyst == 1 &
                   ooph == 0 &
                   mast == 1) {
                        row$riskmod <- list(c("Hysterectomy", "Mastectomy"))
                        row$interAge <- list(as.numeric(c(row$interAgeHyst,
                                                          row$interAgeMast)))
                }
                # Condition 4 (Y, Y, Y)
                if(hyst == 1 &
                   ooph == 1 &
                   mast == 1) {
                        row$riskmod <- list(c("Hysterectomy", "Oophorectomy", 
                                              "Mastectomy"))
                        row$interAge <- list(as.numeric(c(row$interAgeHyst,
                                                          row$interAgeOoph,
                                                          row$interAgeMast)))
                }
                # Condition 5 (N, N, N)
                # Not included b/c default column value is NA, so don't need to 
                # do anything here
                
                # Condition 6 (N, Y, N)
                if(hyst == 0 &
                   ooph == 1 &
                   mast == 0) {
                        row$riskmod <- list("Oophorectomy")
                        row$interAge <- list(as.numeric(row$interAgeOoph))
                }
                # Condition 7 (N, N, Y)
                if(hyst == 0 &
                   ooph == 0 &
                   mast == 1) {
                        row$riskmod <- list("Mastectomy")
                        row$interAge <- list(as.numeric(row$interAgeMast))
                }
                # Condition 8 (N, Y, Y)
                if(hyst == 0 &
                   ooph == 1 &
                   mast == 1) {
                        row$riskmod <- list(c("Oophorectomy", "Mastectomy"))
                        row$interAge <- list(as.numeric(c(row$interAgeOoph,
                                                          row$interAgeMast)))
                }
                
                ped[i, ] <- row
                
        }
        
        # Remove old columns
        ped <- subset(ped, 
                      select = -c(riskmodHyst, interAgeHyst, riskmodOoph, 
                                  interAgeOoph, riskmodMast, 
                                  interAgeMast))
        
        return(ped)
}
