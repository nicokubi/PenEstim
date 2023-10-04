#' Check and make corrections to pedigree
#' 
#' @details 
#' First, `checkFam` checks that all required columns in pedigree are present, 
#' valid, and consistent, printing informative error/warning messages and 
#' making corrections when possible. It also removes family members that are 
#' not connected to the proband(s) and splits the pedigree into connected 
#' sub-families, if necessary. For more details, see the 
#' \code{\link{.pedStructureCheck}} documentation. 
#' 
#' Second, `checkFam` imputes missing current, cancer affection, and death 
#' ages for relatives in `ped`, subject to age-bounding information. For more 
#' details, see the \code{\link{ageImpute}} documentation. 
#' 
#' @param ped A pedigree data frame. (See the description of the `pedigree` 
#' parameter in the \code{\link{PPP}} master function documentation.)
#' @param db A model-specific database returned by \code{\link{buildDatabase}}. 
#' The default is `NULL`, in which case the database will be constructed based 
#' on the cancer information encoded in `ped` and the genes in `GENE_TYPES`.  
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
#' columns will be treated as missing and imputed. `allow.age.zero` should only 
#' be `TRUE` if individuals are confirmed to born with a cancer. Currently, the 
#' `allow.age.zero = TRUE` is not supported behavior for the `PPP` model 
#' functions. 
#' @param impute.missing.ages A logical value indicating whether missing ages 
#' should be multiply imputed. The default is `TRUE`. If set to `FALSE`, an 
#' error will be raised if any individuals have missing current or cancer 
#' affection ages. 
#' @param impute.times A numeric value indicating the number of iterations that 
#' should be run when multiply imputing missing ages. The default is `20`.
#' @param max.iter.tries A numeric value indicating the maximum number of 
#' iterations that should be attempted when multiply imputing missing ages. 
#' (Invalid iterations typically occur when an individual is missing both their 
#' current and cancer ages. Current ages for living individuals are imputed 
#' first; if a valid cancer affection cannot be found, subject to the upper 
#' bound set by the imputed current age, this imputation run will be 
#' discarded.) The default is `NULL`, in which case the maximum number of 
#' tries will be set to five times `iterations`. 
#' @param random.seed The random seed (a numeric value) to set for imputing 
#' missing ages. The default is `42L`. 
#' 
#' @return
#' A list with three components: 
#' * `lms`: A nested list where each component corresponds to the imputed ages 
#' for a connected pedigree in `ped_list`. Each component is itself a list 
#' where each component is a matrix of imputed ages corresponding to `CurAge` 
#' (if one or more individuals in `ped` is missing their current age) or 
#' `AgeX`, where `X` is a cancer for which one or more individuals in `ped` 
#' is missing the affection age. Each matrix has `impute_times` rows; the 
#' column names are the IDs of the individuals for whom ages were imputed. 
#' * `proband`: A numeric value or vector of the unique IDs in the pedigree
#' corresponding to the proband(s). If the `proband` argument was specified by 
#' the user, the returned list component will be identical. Otherwise, the IDs 
#' will correspond to individuals in the `ped` argument where 
#' `ped$isProband == 1`. 
#' * `ped_list`: A list where each element is a connected pedigree that is a 
#' subset of the relatives in `ped`. Typically, `ped_list` has length `1`, 
#' unless there are multiple probands who are related to multiple disconnected 
#' pedigrees. Relatives who are not connected to the proband(s) will be removed. 
#'
#' @family check
#' @family impute
#' @seealso \code{\link{.pedStructureCheck}}, \code{\link{ageImpute}}
#' @examples
#' # Build database for specific model configurations
#' brcadb <- buildDatabase("BRCAPRO")
#' # Run pedigree checks and age imputations
#' checkFam(ped = test_fam_1, db = brcadb)
#' @md
#' @export
checkFam <- function(ped, db = NULL, proband = NULL,
                     unknown.race = UNKNOWN_RACE, 
                     unknown.ancestry = UNKNOWN_ANCESTRY,
                     ignore.proband.germ = FALSE, 
                     allow.age.zero = FALSE, impute.missing.ages = TRUE,
                     impute.times = 20, max.iter.tries = NULL,
                     random.seed = 42L) {
  
  # Build model-specific database if db is not specified
  if (is.null(db)) {
    MS_cancers <- .getCancersFromFam(ped)
    db <- buildDatabase(
      ppd = PanelPRODatabase,
      cancers = .mapCancerNames(short = MS_cancers),
      genes = GENE_TYPES
    )
  } else {
    MS_cancers <- .mapCancerNames(long = db$MS$CANCER)
  }

  # Check pedigree structure
  struc_check <- .pedStructureCheck(ped,
    cancers = MS_cancers,
    proband = proband,
    unknown.race = unknown.race,
    unknown.ancestry = unknown.ancestry,
    ignore.proband.germ = ignore.proband.germ, 
    allow.age.zero = allow.age.zero
  )
  
  # Extract list of pedigrees, probands, and list of relationship lists
  ped_list <- struc_check$ped_list
  proband <- struc_check$proband
  rel_ls <- struc_check$rel_ls

  impute.times <- as.integer(impute.times)
  # If max.iter.tries was not specified, set it to 5*impute.times
  if (is.null(max.iter.tries) == TRUE) {
    max.iter.tries = 5*impute.times
  }
  
  if (impute.missing.ages) {
    # Save previous random state
    old_state <- get_rand_state()
    
    # Set a local random state
    set.seed(as.integer(random.seed))
    
    # Impute missing ages
    lms <- Map(function(ped, rel_l) {
      ageImpute(
        ped = ped, rel_l = rel_l, proband = proband,
        impute_times = impute.times, 
        max_iter_tries = max.iter.tries,
        penet_db = db$penet$penet_c,
        penet_cbc_db = db$penet$penet_cbc
      )
    }, ped = ped_list, rel_l = rel_ls)
    
    # Restore previous random state
    on.exit(set_rand_state(old_state))
  } else {
    # Otherwise, check that there are no missing ages
    lms <- vector("list", 1)
    ped_concat <- do.call(rbind, ped_list)
    # Error if there are missing affection ages
    for (cc in MS_cancers) {
      if (any(ped_concat[ped_concat[[paste0("isAff", cc)]] == 1, 
                         paste0("Age", cc)] == -999)) {
        rlang::abort(
          "Affection ages cannot be missing when impute.missing.ages = FALSE.", 
          level = "CannotImpute")
      }
    }
    # Error if CurAge is missing
    if (any(ped_concat[ped_concat$ID %in% proband, "CurAge"] == -999)) {
      rlang::abort(
        "Current ages of probands cannot be missing when impute.missing.ages = FALSE.", 
        level = "CannotImpute")
    }
  }

  out <- list(lms = lms, proband = proband, ped_list = ped_list)
  return(out)
}

