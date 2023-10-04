#' Subset penetrances
#'
#' Subsets penetrances for carriers from a model-specific database, calculates 
#' the penetrances for non-carriers based on SEER and carrier penetrances, 
#' applies risk modifiers if requested, and calculates the MAXAGE+1 penetrance 
#' as 1 minus the sum of the other penetrances. The resulting array stores the 
#' penetrance for each family member. 
#'
#' @param fam A checked pedigree data frame returned by \code{\link{checkFam}}. 
#' @param lm.ages A list of imputed ages returned by \code{\link{checkFam}}. 
#' @param proband A numeric value or vector of the unique IDs in `ped` for 
#' whom to return cancer penetrances. 
#' @param db A model-specific database returned by \code{\link{buildDatabase}}. 
#' The default is `NULL`, in which case the penetrances in `sub_dens` will be 
#' used directly instead of subsetted from `db`. One of `db` and `sub_dens` 
#' must be non-`NULL`. 
#' @param consider.modification A logical value indicating whether or not to 
#' consider interventions that modify the penetrance. The default is `FALSE`. 
#' @param net A logical value indicating if net or crude penetrances should be 
#' returned. The default is `TRUE`, for net penetrances. 
#' @param doc A logical value indicating whether the pentrances for death by 
#' other causes should be calculated. If `FALSE`,
#' the cancer penetrances will be calculated instead. The default is `FALSE`. 
#'
#' @details
#' The cancer carrier penetrances are obtained directly from `db`. The 
#' non-carrier cancer penetrances are derived using the SEER penetrances stored 
#' in `db` and the cancer carrier penetrances. Both the carrier and non-carrier 
#' penetrances for death by other causes are derived from the cancer 
#' penetrances and the SEER penetrances for death by other causes. The SEER 
#' penetrances are taken to be population-level, i.e., marginalized over all 
#' the genotypes. Thus, they are assumed to be weighted averages of the 
#' non-carrier penetrances and carrier penetrances, with weights based on the 
#' population-level genotype probabilities.
#'
#' @return A numeric array of subsetted cancer penetrances, possibly modified 
#' by risk modifiers, with the following dimensions:
#' 
#' * `ID`: The IDs of the probands in the family (from `ID` in `fam`). 
#' * `Cancer`: The cancers included in the model. 
#' * `Gene`: The genes included in the model. 
#' * `Age`: Ages 1-`PPP:::MAXAGE+1`.
subsetCancerPenetrance <- function(fam, lm.ages, db, proband = NULL, 
                                   consider.modification = FALSE,
                                   net = TRUE, doc = FALSE) {
  
  # Subset the net penetrance and drop levels, while preserving dimension names
  if (doc == FALSE) {
    net.crude <- ifelse(net == TRUE, "Net", "Crude")
    penet.c <- subset_array(db$penet$penet_c, axisName = "PenetType", 
                            choice = net.crude, drop = FALSE)
    penet.nc <- subset_array(db$penet$penet_nc, axisName = "PenetType", 
                             choice = net.crude, drop = FALSE)
    penet.c <- abind::adrop(penet.c, drop = 6)
    penet.nc <- abind::adrop(penet.nc, drop = 6)
  } else {
    penet <- db[["doc"]]
    # Separate lists for carriers and non-carriers
    penet.c <- penet[["doc_c"]]
    penet.nc <- penet[["doc_nc"]]
  }
  cancers <- dimnames(penet.c)$Cancer
  
  # check if CBC should be incorporated and prepare CBC penetrance array
  cbc <- ifelse(doc == FALSE & 
                  "isAffCBC" %in% colnames(fam) & 
                  !is.null(db$penet$penet_cbc), 
                TRUE, FALSE)
  if(cbc == TRUE){
    cancers <- c(cancers, "Contralateral")
    penet.cbc <- subset_array(arrayInput = db$penet$penet_cbc, 
                              axisName = "PenetType", 
                              choice = "Net", drop = FALSE)
    penet.cbc <- abind::adrop(penet.cbc, drop = 4)
    cbcDnames <- dimnames(penet.cbc)
    dim(penet.cbc) <- c(1, dim(penet.cbc))
    dimnames(penet.cbc) <- append(cbcDnames, list(Cancer = "Contralateral"), 
                                  after = 0)
  }
  
  # Check that the sex levels from the database match the pedigree
  if (any(is.na(fam$Sex))) rlang::abort("sex contains NA")
  stopifnot(all(fam$Sex %in% dimnames(penet.c)$Sex))
  
  # Subset family to only include probands, if needed
  if (!is.null(proband)) {
    fam = fam[fam$ID %in% proband,]
  }
  
  # Initialize array container for family members' penetrances
  dn_fam_pen <- list(
    ID = fam$ID,
    Cancer = cancers,
    Gene = c("noncarrier", dimnames(penet.c)$Gene),
    Age = 1:MAXAGE
  )
  fam_pen <- array(dim = lengths(dn_fam_pen), dimnames = dn_fam_pen)
  
  # Helper function to get carrier and noncarrier penetrances customized for 
  # each individual. 
  .align_fam <- function(pen, type, ped = fam, imp.ages = NULL) {
    
    # Set up information for array container based on selected type 
    # (carrier or noncarrier)
    type <- match.arg(type, c("c", "nc", "cbc")) 
    if (type == "c") {
      by_cols <- c("Race", "Sex")
      ped_cols <- c("race", "Sex")
    } else if (type == "nc") {
      ped$nc_names <- paste0("noncarrier_", ped$Ancestry, ".", ped$race)
      by_cols <- c("Gene", "Sex")
      ped_cols <- c("nc_names", "Sex")
    } else if(type == "cbc"){ # create CBC parameter data frame
      bc.ages <- first.birth.age <- had.1bc <- bc.types <- fdr.bc <- er.status <- 
        antiE <- HRp <- bDens <- bc1.t.cat <- as.numeric(rep(NA, nrow(ped)))
      CBCparams <- data.frame(ids = ped$ID, sex = ped$Sex, 
                              had.1bc, bc.ages, bc.types, fdr.bc, er.status, 
                              first.birth.age, antiE, HRp, bDens, bc1.t.cat)
      for(id in unique(ped$ID)){
        cbc.params <- getCBCParams(ped, id, penet.cbc, imp.ages)
        CBCparams[CBCparams$ids == id, 
                  c("had.1bc","bc.ages","bc.types","fdr.bc","er.status",
                    "first.birth.age","antiE","HRp","bDens","bc1.t.cat","race")] <- 
          c(cbc.params$HadBC, cbc.params$FirstBCAge, cbc.params$FirstBCType, 
            cbc.params$FDRBC, cbc.params$ER, cbc.params$FirstBirth,
            cbc.params$AntiEstrogen, cbc.params$HRPreneoplasia, 
            cbc.params$BreastDensity, cbc.params$BC1TSize, 
            cbc.params$race)
      }
      c_by_cols <- c("Sex", "FirstBCAge")
      c_cbc_cols <- c("sex", "bc.ages")
    }
    
    if(type != "cbc"){
      fam_pen <- lapply(seq(nrow(ped)), function(i) {
        subset_array(pen, by_cols, ped[i, ped_cols], drop = FALSE)
      })
    } else if(type == "cbc") {
      fam_pen <- lapply(seq(nrow(ped)), function(i) {
        
        # check if 1st BC occurred
        if(CBCparams$had.1bc[i] == 1){
          
          # noncarrier gene penetrances
          nc.cbc.vec <- CBCnoncarrierCrude(diag.age = CBCparams$bc.ages[i], 
                                           bc1.type = CBCparams$bc.types[i], 
                                           fdr.status = CBCparams$fdr.bc[i], 
                                           ER.status = CBCparams$er.status[i], 
                                           birth1.age = CBCparams$first.birth.age[i], 
                                           antiEst = CBCparams$antiE[i], 
                                           HighRiskP = CBCparams$HRp[i], 
                                           BD = CBCparams$bDens[i],
                                           bc1.t.size = CBCparams$bc1.t.cat[i],
                                           sex = CBCparams$sex[i],
                                           race = CBCparams$race[i])
          tmp.nc.cbc.arr <- array(nc.cbc.vec, 
                                  dim = c(1, 1, length(nc.cbc.vec)),
                                  dimnames = list(Cancer = "Contralateral",
                                                  Gene = "noncarrier",
                                                  Age = 1:MAXAGE))
          
          # non-CBC associated gene penetrances
          nonCBCgenes <- setdiff(dimnames(penet.c)$Gene, 
                                 unlist(CBC_GENE_VARIANT_TYPES))
          nc.cbc.arr <- tmp.nc.cbc.arr
          if(length(nonCBCgenes) > 0){
            for(ncbcg in nonCBCgenes){
              dimnames(tmp.nc.cbc.arr)[[2]] <- ncbcg
              nc.cbc.arr <- abind::abind(nc.cbc.arr, tmp.nc.cbc.arr,
                                         along = 2, use.dnns = TRUE)
            }
          }
          
          # CBC associated gene penetrances
          c.cbc.arr <- subset_array(pen, c_by_cols, 
                                    CBCparams[i, c_cbc_cols], drop = FALSE)
          c.cbc.arr <- abind::adrop(c.cbc.arr, drop = 3:4)
          
          # combine all CBC penetrances into one array
          cbc.arr <- abind::abind(nc.cbc.arr, c.cbc.arr, along = 2, use.dnns = TRUE)
          
        } else { # all 0 penetrances if no previous BC
          cbc_dnames <- list(
            Cancer = "Contralateral",
            Gene = c("noncarrier", dimnames(penet.c)$Gene),
            Age = 1:MAXAGE
          )
          cbc.arr <- array(0, dim = lengths(cbc_dnames), dimnames = cbc_dnames)
        }
        cbc.arr
      })
    }
    names(fam_pen) <- ped$ID # Assign names by unique IDs
    
    # Fill in array with penetrances
    arr <- do.call(abind::abind, c(fam_pen, along = -1))
    if (type == "c") {
      return(abind::adrop(arr, drop = c(4, 5)))
    } else if (type == "nc") {
      arr <- abind::adrop(arr, drop = c(4, 5))
      dimnames(arr)[[3]] <- c("noncarrier")
      return(arr)
    } else if (type == "cbc"){
      return(arr)
    }
  }
  
  # Fill in carrier and noncarrier penetrances for each family member
  abind::afill(fam_pen) <- .align_fam(penet.c, type = "c")
  abind::afill(fam_pen) <- .align_fam(penet.nc, type = "nc")
  if(cbc == TRUE){ 
    abind::afill(fam_pen) <- 
      .align_fam(penet.cbc, type = "cbc", imp.ages = lm.ages) 
  }
  
  # Initialize array container for family members' penetrances
  dn_sub_dens <- list(
    ID = fam$ID,
    Cancer = cancers,
    Gene = c("noncarrier", dimnames(penet.c)$Gene),
    Age = 1:(MAXAGE + 1)
  )
  sub_dens = array(dim = lengths(dn_sub_dens), dimnames = dn_sub_dens)
  
  # Modify cancer penetrances using risk modifiers
  if (consider.modification == TRUE && doc == FALSE) {
    sub_dens[,,,1:MAXAGE] <- riskmod_perfam(fam_pen, fam, riskmod = db$riskmod)
  } else {
    sub_dens[,,,1:MAXAGE] <- fam_pen
  }
  
  # Calculate the penetrances for MAXAGE+1
  sub_dens[,,,MAXAGE+1] = 1 - apply(sub_dens, c(1, 2, 3), sum, na.rm = TRUE)
  
  return(sub_dens)
}


#' Calculate penetrances
#'
#' Calculates the density and survival penetrance functions (for either cancer 
#' or for death by other causes) for each proband in the pedigree, based on 
#' the genes and cancers in the model.
#'
#' @param fam A checked pedigree data frame returned by \code{\link{checkFam}}. 
#' @param proband A numeric value or vector of the unique IDs in `ped` for 
#' whom to return cancer penetrances. 
#' @param db A model-specific database returned by \code{\link{buildDatabase}}. 
#' The default is `NULL`, in which case the penetrances in `sub_dens` will be 
#' used directly instead of subsetted from `db`. One of `db` and `sub_dens` 
#' must be non-`NULL`. 
#' @param sub_dens A numeric array of subsetted cancer penetrances returned by 
#' \code{\link{subsetCancerPenetrance}}. The default is `NULL`, in which case the 
#' penetrances will be subsetted from `db`. One of `db` and `sub_dens` must be 
#' non-`NULL`. 
#' @param PGs Possible genotypes in both list and data frame format, returned 
#' by \code{\link{.getPossibleGenotype}}. 
#' @param direct_fill_PGs A vector of genotype names with no more than 1 
#' mutation (including `"noncarrier"`), corresponding to `PGs`. 
#' @param multi_PGs A vector of genotype names with 2 or more simultaneous 
#' mutations, corresponding to `PGs`.  
#' @param multi_muts A list where each component corresponds to a genotype in 
#' `multi_PGs`, represented as a vector of gene names to indicate the mutations. 
#' mutations, corresponding to `PGs`. 
#' @param consider.modification A logical value indicating whether or not to 
#' consider interventions that modify the penetrance. The default is `FALSE`. 
#' @param net A logical value indicating if net or crude penetrances should be 
#' returned. The default is `TRUE`, for net penetrances. 
#' @param doc A logical value indicating whether the pentrances for death by 
#' other causes should be calculated. If `FALSE`,
#' the cancer penetrances will be calculated instead. The default is `FALSE`.
#' @param lm.ages A list of imputed ages returned by \code{\link{checkFam}}. 
#' Default is `NULL`; only needed when CBC is in the model and a BC age is missing.
#'
#' @return A list of density (`Dens`) and survival (`Sur`) penetrance arrays 
#' for the probands, each of which contains the following dimensions:
#' 
#' * `ID`: The IDs of the probands in the family (from `ID` in `fam`). 
#' * `cancers`: The cancers included in the model. 
#' * `genotype`: All possible genotypes, with a maximum number of simultaneous 
#' mutations based on `max_mut`. 
#' * `ages`: Ages 1-`PPP:::MAXAGE`.
calcCancerPenetrance <- function(fam, proband, db = NULL, sub_dens = NULL, 
                                 PGs, direct_fill_PGs, multi_PGs, multi_muts, 
                                 consider.modification = FALSE,
                                 net = TRUE, doc = FALSE, lm.ages = NULL) {
  
  if (is.null(db) && is.null(sub_dens)) {
    # Check that either db or sub_dens is specified
    rlang::abort("No db or sub_dens argument specified for calcCancerPenetrance.")
  } else if (is.null(sub_dens)) {
    # If sub_dens was not specified, create it from db
    sub_dens = subsetCancerPenetrance(
      fam = fam, lm.ages = lm.ages, proband = proband, db = db, 
      consider.modification = consider.modification,
      net = net, doc = doc)
  } 
  
  # Construct the final cancer penetrance density array
  cancerDens <- array(0,
    dim = c(length(proband), length(dimnames(sub_dens)[[2]]), 
            nrow(PGs$df), MAXAGE + 1),
    dimnames = list(
      ID = proband, cancers = dimnames(sub_dens)[[2]],
      genotypes = PGs$list, ages = seq(MAXAGE + 1)
    )
  )
  # Construct the final cancer penetrance survival array
  cancerSurv <- cancerDens
  
  # Fill in penetrances for single mutations directly
  cancerDens[, , direct_fill_PGs, ] <- 
    sub_dens[, , direct_fill_PGs, ]
  cancerSurv[, , direct_fill_PGs, ] <- 
    1 - drop(aperm(apply(cancerDens[, , direct_fill_PGs, , drop = FALSE], 
                         c(1, 2, 3), cumsum), c(2, 3, 4, 1)))

  # Calculate penetrances for genotypes with multiple mutations
  for (i in seq_along(multi_PGs)) {
    cancerSurv[, , multi_PGs[i], ] <- 
      drop(apply(cancerSurv[, , multi_muts[[i]], , drop = FALSE], 
                 c(1, 2, 4), prod))
    one_array <- array(1, dim = c(dim(cancerSurv)[c(1, 2)], 1))
    cancerSur_ones <- 
      abind::abind(one_array, 
                   array(cancerSurv[, , multi_PGs[i], ], 
                         dim = dim(cancerSurv)[c(1, 2, 4)]), along = 3)
    cancerDens[, , multi_PGs[i], ] <- 
      drop((-1) * aperm(apply(cancerSur_ones, c(1, 2), diff), c(2, 3, 1)))
  }
  
  return(list(Dens = cancerDens, Sur = cancerSurv))
}


#' Calculate an individual's cancer survival probabilities
#'
#' Calculates the cancer survival functions for an individual. 
#'
#' @param ID Unique ID corresponding to an ID in `sub_dens`. 
#' @param curage Individual `ID`'s (non-missing) current age. 
#' @param cancers A vector of cancer names for cancers that individual `ID` did 
#' not develop by their `curage`. 
#' @param sub_dens A numeric array of subsetted cancer penetrances returned by 
#' \code{\link{subsetCancerPenetrance}}. 
#' @param PGs Possible genotypes in both list and data frame format, returned 
#' by \code{\link{.getPossibleGenotype}}. 
#' @param direct_fill_PGs A vector of genotype names with no more than 1 
#' mutation (including `"noncarrier"`), corresponding to `PGs`. 
#' @param multi_PGs A vector of genotype names with 2 or more simultaneous 
#' mutations, corresponding to `PGs`.  
#' @param multi_muts A list where each component corresponds to a genotype in 
#' `multi_PGs`, represented as a vector of gene names to indicate the mutations. 
#' mutations, corresponding to `PGs`. 
#' 
#' @return A numeric vector with length equal to the number of possible 
#' genotypes. Each element corresponds to the probability that individual `ID` 
#' does not develop any of the cancers in `cancers` by their `curage`, given a 
#' possible genotype in `PGs`. If `cancers` is an empty vector, `1` will be 
#' returned. 
calcCancerSurv = function(ID, curage, cancers, sub_dens, PGs, 
                          direct_fill_PGs, multi_PGs, multi_muts) {
  
  # If there are no cancers, return 1 (no contribution)
  if (length(cancers) == 0) {
    return(1)
  }
  
  # Initialize the vector for storing the probability of surviving all of the 
  # cancers, for each possible genotype
  cancerSurv <- array(0,
                      dim = nrow(PGs$df),
                      dimnames = list(genotypes = PGs$list))
  
  # Calculate survivals for single mutations
  cancerSurv[direct_fill_PGs] <- apply(1 - abind::adrop(
    aperm(
      apply(sub_dens[as.character(ID), cancers, direct_fill_PGs, 1:curage, drop = FALSE], 
            c(1, 2, 3), sum), 
      c(2, 3, 1)), 
    3), 2, prod)
  
  # Calculate survivals for genotypes with multiple mutations
  if (length(multi_PGs) > 0) {
    cancerSurv[multi_PGs] <- 
      sapply(multi_muts, function(muts) {
        prod(cancerSurv[muts, drop = FALSE])
      })
  }
  
  return(cancerSurv)
}


#' Calculate an individual's cancer density probabilities
#'
#' Calculates the cancer density functions for an individual. 
#'
#' @param ID Unique ID corresponding to an ID in `sub_dens`. 
#' @param affages Individual `ID`'s (non-missing) ages of affection 
#' corresponding to `cancers`. 
#' @param cancers A vector of cancer names for cancers that individual `ID` 
#' developed at `affages`. 
#' @param sub_dens A numeric array of subsetted cancer penetrances returned by 
#' \code{\link{subsetCancerPenetrance}}. 
#' @param PGs Possible genotypes in both list and data frame format, returned 
#' by \code{\link{.getPossibleGenotype}}. 
#' @param direct_fill_PGs A vector of genotype names with no more than 1 
#' mutation (including `"noncarrier"`), corresponding to `PGs`. 
#' @param multi_PGs A vector of genotype names with 2 or more simultaneous 
#' mutations, corresponding to `PGs`.  
#' @param multi_muts A list where each component corresponds to a genotype in 
#' `multi_PGs`, represented as a vector of gene names to indicate the mutations. 
#' mutations, corresponding to `PGs`. 
#' 
#' @return A numeric vector with length equal to the number of possible 
#' genotypes. Each element corresponds to the probability that individual `ID` 
#' develops the `cancers` at their corresponding `affages`, given a 
#' possible genotype in `PGs`. If `cancers` is an empty vector, `1` will be 
#' returned. 
calcCancerDens = function(ID, affages, cancers, sub_dens, PGs, 
                          direct_fill_PGs, multi_PGs, multi_muts) {
  
  # If there are no cancers, return 1 (no contribution)
  if (length(cancers) == 0) {
    return(1)
  }
  
  # Check that the the number of affected ages matches the number of cancers
  if (length(cancers) != length(affages)) {
    rlang::abort("Number of affected ages does not match the number of cancers.")
  }
  
  # Initialize the vector for storing the probability of developing all of the 
  # cancers, for each possible genotype
  cancerDens <- array(0,
                      dim = nrow(PGs$df),
                      dimnames = list(genotypes = PGs$list))
  
  # Fill in penetrances for single mutations directly
  cancerDens[direct_fill_PGs] <- apply(sapply(1:length(cancers), function(i) {
    sub_dens[as.character(ID), cancers[i], 
                direct_fill_PGs, affages[i], 
                drop = FALSE]
  }), 1, prod)
  
  # Construct an array for storing survival probabilities for 
  # the affected ages - 1 and the affected ages 
  cancerSurv <- array(0,
                      dim = c(length(cancers), 
                              nrow(PGs$df), 2),
                      dimnames = list(
                        cancers = cancers,
                        genotypes = PGs$list, ages = c("Age-1", "Age"))
                      )
  
  # Calculate survivals for single mutations
  cancerSurv[, direct_fill_PGs,] <- 
    1 - do.call(abind::abind, c(lapply(1:length(cancers), function(i) {
      if (affages[i] > 1) { 
        # If the affected age is more than 1, calculate survivals for the 
        # the affected age - 1 and the affected age
        abind::adrop(
          aperm(
            apply(
              sub_dens[as.character(ID), cancers[i], 
                       direct_fill_PGs, 1:(affages[i]), 
                       drop = FALSE], 
              c(1, 2, 3), cumsum)[(affages[i]-1):(affages[i]),,,, drop = FALSE], 
            c(2, 3, 4, 1)), 
          1)
      } else {
        # If the affected age is 1, calculate the survival for the affected age 
        # and use an array of 1s for the affected age - 1
        abind::abind(
          array(0, dim = c(1, length(direct_fill_PGs), 1)), 
          abind::adrop(
            sub_dens[as.character(ID), cancers[i], 
                     direct_fill_PGs, affages[i], 
                     drop = FALSE], 1), 
        along = 3)
      }
    }), along = 1))

  # Calculate penetrances for genotypes with multiple mutations
  if (length(multi_PGs) > 0) {
    # Survival
    cancerSurv[, multi_PGs, ] <- 
      aperm(do.call(abind::abind, 
                    c(lapply(multi_muts, function(muts) {
                      matrix(apply(cancerSurv[, muts,, drop = FALSE], 
                                   c(1, 3), prod), ncol = 2)
                    }), along = 3)), 
            c(1, 3, 2))
    # Density
    cancerDens[multi_PGs] = 
      sapply(multi_PGs, function(muts) {
        prod((-1) * apply(cancerSurv[, muts,, drop = FALSE], c(1, 2), diff))
      })
  }
  
  return(cancerDens)
}
