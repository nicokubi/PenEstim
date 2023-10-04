#' Calculate future risk
#'
#' Calculates the future risk of cancer for the proband(s), based on the 
#' posterior carrier probabilities. 
#' 
#' @param ped A checked pedigree data frame returned by \code{\link{checkFam}}. 
#' @param proband A numeric value or vector of the unique IDs in `ped` for 
#' whom to return future risks. 
#' @param MS_cancers The cancers in the model, obtained from the 
#' model-specific database. 
#' @param cplist A list containing the cancer penetrances for the probands in 
#' `ped`, returned by \code{\link{calcCancerPenetrance}}. If `net == TRUE`, 
#' these should be net penetrances; if `net == FALSE`, these should be crude. 
#' @param doclist A list containing the penetrances for the death by other 
#' causes for `ped`, returned by \code{\link{calcCancerPenetrance}}. 
#' @param posterior_probs The genotypic distribution matrix of the family 
#' members in `proband`, returned by \code{\link{pp.peelingParing}}. 
#' @param net A logical value indicating whether to use the net cancer 
#' penetrances to estimate future risk. The default is `FALSE`, in which case 
#' the crude penetrances will be used. 
#' @param age.by The number of years between each future risk prediction. 
#' The default is `5`; if the proband's current age is 30, risk predictions 
#' will be calculated at ages 35, 40, 45, etc.
#' @param rr.bcrat A data frame with 3 columns: the IDs of the probands ("ID"),
#' the probands' BCRAT relative risks for age < 50 ("rr1"), and the probands' BCRAT
#' relative risks for age >= 50 ("rr2"). When `net = FALSE` and breast 
#' cancer is included in the model, these values will modify the future risk 
#' calculation. The default is `NULL`, in which the BCRAT relative risks will 
#' not be incorporated. 
#' @param rr.pop A data frame of race-specific for estimates of 1 /
#' (1 - population attributable risk). The default is `rr.ref`,
#' which was estimated from NHIS 2015.
#' @param cbc_penets an array of CBC carrier penetrances as returned by the 
#' penet$penet_cbc component of \code{\link{buildDatabase}}.
#' @param lm.ages imputed ages for the pedigree as returned by 
#' \code{\link{checkFam}}.
#' @return
#' A list of arrays, one for each proband, storing the future risks of 
#' developing the cancers specified by the model. 
calcFutureRisk <- function(ped, proband, MS_cancers,
                           cplist, doclist, 
                           posterior_probs, net = FALSE, age.by = 5, 
                           rr.bcrat = NULL, rr.pop = rr.ref,
                           cbc_penets, lm.ages) {
  
  # Identify the intersect of cancers in the pedigree and in the model
  cancers <- intersect(.mapCancerNames(short = .getCancersFromFam(ped)),
                       MS_cancers)
  
  # Helper function to calculate future risk for a given individual and cancer
  calc.fr.ij <- function(pen.can, pen.death, age.cur, ages.risk, probs, net, 
                         rr.nc, cbc.params) {
    
    # Initialize container for storing future risk at different ages
    fr.overall <- rep(NA, length(ages.risk))
    
    # use separate future risk calculations for CBC and non-CBC cancers
    if(length(cbc.params) == 0){
      if (net == TRUE) { # Net future risk
        sur.can <- 1 - safe_apply(pen.can[, ages = 1:age.cur], 1, sum)
        for (i in 1:length(ages.risk)) {
          try(fr.geno <- safe_apply(pen.can[, ages = (age.cur + 1):ages.risk[i]], 
                                    1, sum) / sur.can)
          fr.overall[i] <- sum(probs * fr.geno)
        }
      } else { # Crude future risk
        # P(T > t) = P(T > t, J = C) + P(T > t, J = D) =
        # \sum_{s=t+1}^MAXAGE P(T = s, J = C) +
        # \sum_{s=t+1}^MAXAGE P(T = s, J = D) +
        # P(T > t, J = D) =
        # \sum_{s=t+1}^MAXAGE P(T = s, J = C) + \sum_{s=t+1}^MAXAGE P(T = s, J = D) +
        # 1 - \sum_{s=1}^MAXAGE P(T = s, J = C) - \sum_{s=1}^MAXAGE P(T = s, J = D)
        # Here we assume that you can only get cancer from ages 1 to MAXAGE, but you
        # can die after MAXAGE
  
        # If using BRCAPRO+BCRAT, convert the penetrance to a hazard,
        # modify the hazard by the risk ratio, and convert back to a penetrance
        if (!is.null(rr.nc)) {
          sur.all.t <- 1 - cumsum(pen.can["noncarrier", , drop = FALSE]) -
            cumsum(pen.death["noncarrier", , drop = FALSE])
          haz.can.nc <- pen.can["noncarrier", , drop = FALSE] /
            c(1, sur.all.t[-length(sur.all.t)])
          # Separate relative risks for age < 50 and age >= 50
          haz.can.nc[1:49] <- 1 - (1 - haz.can.nc[1:49])^rr.nc[1]
          haz.can.nc[50:MAXAGE] <- 1 - (1 - haz.can.nc[50:MAXAGE])^rr.nc[2]
          pen.can.nc <- haz.can.nc * c(1, sur.all.t[-length(sur.all.t)])
          pen.can["noncarrier", ] <- pen.can.nc
        }
  
        # When pen.death are NAs due to sex-mismatched cancers, convert them to 0
        pen.death[is.na(pen.death)] <- 0
  
        sur.all <- 
          apply(pen.can[, ages = (age.cur + 1):MAXAGE, drop = FALSE], 1, sum) +
          apply(pen.death[, ages = (age.cur + 1):MAXAGE, drop = FALSE], 1, sum) +
          1 - 
          apply(pen.can, 1, sum) - 
          apply(pen.death, 1, sum)
        for (i in 1:length(ages.risk)) {
          fr.geno <- 
            apply(pen.can[, ages = (age.cur + 1):ages.risk[i], drop = FALSE], 1, sum) / 
            sur.all
          fr.overall[i] <- sum(probs * fr.geno)
        }
      }
      
      # CBC future risk (which always needs Net penetrances as an input)
    } else {
      
      # calculate matrix of crude CBC future risks for all genotypes and ages
      all.cbc.frs <- fr.cbc(pen.cbc.net = pen.can, 
                            pen.d.crude = pen.death, 
                            pen.bc1.net = cbc.params$net.pen.BC, 
                            age.bc1 = cbc.params$FirstBCAge, 
                            age.current = age.cur,
                            return.NA = ifelse(cbc.params$HadBC == 0, TRUE, FALSE))
      
      # find future risk at each age interval
      for (a.r in 1:length(ages.risk)) {
        fr.geno <- all.cbc.frs[ ,ages.risk[a.r]]
        fr.overall[a.r] <- sum(probs * fr.geno)
      }
    }
    return(fr.overall)
  }

  # Initialize list for storing future risk arrays for each proband
  future.risk <- stats::setNames(
    vector("list", length(proband)),
    proband
  )
  
  for (i in 1:length(proband)) {
    if (ped$isDead[ped$ID == proband[i]] == 1) { 
      # Proband is dead
      future.risk[[i]] <- paste0("ID ", proband[i], " is dead.")
    } else if (ped$CurAge[ped$ID == proband[i]] == MAXAGE) { 
      # Proband is alive and has current age == MAXAGE
      future.risk[[i]] <- paste0("ID ", proband[i], 
                                 " has current age equal to the maximum age supported, ", 
                                 MAXAGE, ".")
    } else { 
      # Proband is alive and has current age < MAXAGE
      # Proband's current age
      age.cur.ij <- ped$CurAge[ped$ID == proband[i]]
      # Generate ages at which to compute future risk
      if (age.cur.ij + age.by > MAXAGE) {
        ages.risk.ij <- MAXAGE
      } else {
        ages.risk.ij <- seq(age.cur.ij + age.by, MAXAGE, age.by)
      }
      
      # Proband's carrier probabilities
      probs.ij <- posterior_probs[row.names(posterior_probs) == proband[i], ]
      
      # Initialize container for storing proband's future risks
      future.risk[[i]] <- stats::setNames(
        data.frame(matrix(NA, length(ages.risk.ij), length(cancers) + 1)),
        c("ByAge", cancers)
      )
      future.risk[[i]]$ByAge <- ages.risk.ij

      for (j in 1:length(cancers)) { # Iterate through cancers
        cancer_short <- .mapCancerNames(long = cancers[j])
        if (ped[ped$ID == proband[i], 
                names(ped) == paste0("isAff", cancer_short)] == 1) {
          # Return NAs if proband already has the cancer
          future.risk[[i]][, j + 1] <- rep(NA, length(ages.risk.ij))
        } else if (cancer_short == "CBC" && 
                   ped[ped$ID == proband[i], "isAffBC"] != 1) {
          # Return NAs if the cancer is CBC and the proband has not developed BC
          future.risk[[i]][, j + 1] <- rep(NA, length(ages.risk.ij))
        } else { # Otherwise, compute future risk
          # Subset death by other causes
          pen.death.ij <- abind::adrop(doclist$Dens[ID = as.character(proband[i]),
                                                    cancers = cancers[j], ,
                                                    ages = 1:MAXAGE, 
                                                    drop = FALSE], 
                                       drop = c(1,2))
          
          # Subset cancer penetrances
          pen.can.ij <- abind::adrop(cplist$Dens[ID = as.character(proband[i]),
                                                 cancers = cancers[j], ,
                                                 ages = 1:MAXAGE, 
                                                 drop = FALSE],
                                     drop = c(1,2))
          
          if (net == TRUE) {
            # Skip BRCAPRO+BCRAT if using net risk
            rr.nc <- NULL
          } else {
            # Use BRCAPRO+BCRAT for breast cancer if rr.bcrat is
            # provided and the user asks for crude risk
            if (cancer_short == "BC" & !is.null(rr.bcrat)) {
              if (ped$race[ped$ID == proband[i]] == "AIAN") {
                rr.pop$noncarrier <- rr.pop$native
              } else if (ped$race[ped$ID == proband[i]] == "Asian") {
                rr.pop$noncarrier <- rr.pop$asian
              } else if (ped$race[ped$ID == proband[i]] == "Black") {
                rr.pop$noncarrier <- rr.pop$black
              } else if (ped$race[ped$ID == proband[i]] == "White") {
                rr.pop$noncarrier <- rr.pop$white
              } else if (ped$race[ped$ID == proband[i]] %in%
                c("Hispanic", "WH", "WNH")) {
                rr.pop$noncarrier <- rr.pop$hispanic
              }
              rr.nc <- rr.bcrat[rr.bcrat$ID == proband[i], ] / rr.pop$noncarrier
            } else {
              rr.nc <- NULL
            }
          }
          
          # get additional CBC parameters (age of 1st BC & BC net penetrance)
          cbc.params <- list()
          if(cancers[j] == "Contralateral"){ 
            cbc.params <- getCBCParams(ped, proband[i], cbc_penets, lm.ages)
            
            # find net penetrance of BC based on weighted average of penetrances
            # using posterior probabilities of possible genotypes
            cbc.params$net.pen.BC <- 
              colSums(abind::adrop(cplist$Dens[ID = as.character(proband[i]),
                                               cancers = "Breast", ,
                                               ages = 1:MAXAGE, 
                                               drop = FALSE],
                                   drop = c(1,2)) * probs.ij) / length(probs.ij)
            
          }
          
          # Calculate future risk for each individual
          future.risk[[i]][, j + 1] <- calc.fr.ij(pen.can.ij, pen.death.ij,
                                                  age.cur.ij, ages.risk.ij,
                                                  probs.ij,
                                                  net = net,
                                                  rr.nc = rr.nc,
                                                  cbc.params = cbc.params)
        }
      }
    }
  }
  return(future.risk)
}


#' Calculate Future Risk of Contralateral Breast Cancer (CBC) 
#' 
#' Calculates the future risk for developing cancer in the second breast at a 
#' given future age for a subject previously affected with breast cancer. This
#' calculation considers both survival from death from non-cancer causes 
#' and survival from CBC.
#' 
#' @param pen.cbc.net Matrix of net CBC penetrances where row names are 
#' genotypes and column names are ages from 1 to MAXAGE.
#' @param pen.d.crude Matrix of crude death from non-cancer cause penetrances
#' where rows names are genotypes and column names are ages from 1 to MAXAGE.
#' @param pen.bc1.net 1 row matrix or a vector of net 1st Breast Cancer (BC) 
#' penetrances from age 1 to MAXAGE. Penetrances should be based on a weighted 
#' average of possible genotypes.
#' @param age.bc1 Age at which the proband developed 1st BC; numeric value.
#' @param age.current The proband's current age; numeric value.
#' @param return.NA logical which indicates if a matrix of all NAs should be 
#' returned instead of calculated values. Useful those without a first BC.
#' 
#' @return 
#' A matrix of future risks where the row names are genotypes and the column
#' names are ages 1 to MAXAGE.
fr.cbc <- function(pen.cbc.net, pen.d.crude, pen.bc1.net, age.bc1, 
                   age.current, return.NA){
  
  # instantiate future risk storage by genotype
  cbc.frs <- matrix(0, nrow = nrow(pen.cbc.net), ncol = MAXAGE)
  rownames(cbc.frs) <- rownames(pen.cbc.net)
  colnames(cbc.frs) <- as.character(1:MAXAGE)
  
  # return all NAs (usually if never had BC)
  if(return.NA == TRUE){ 
    cbc.frs[,] <- NA
    return(cbc.frs) 
  }
  
  ## helper: convert crude death penetrance to net
  crude2net.d <- function(pen_d_crude, pen_bc1_net){
    pen_d_crude / (1 - cumsum(pen_bc1_net))
  }
  
  ## helper: calculate overall survival using non-cause-specific framework 
  # (denominator in future risk calculation)
  surv.all <- function(pen_cbc_net, pen_d_crude, pen_bc1_net, age_bc1){
    
    # survival from CBC (net)
    surv.cbc.net <- rep(1, MAXAGE)
    surv.cbc.net[(age_bc1+1):MAXAGE] <- 1 - cumsum(pen_cbc_net[(age_bc1+1):MAXAGE])
    
    # survival from death from other causes (net)
    surv.d.net <- rep(1, MAXAGE)
    pen.d.net <- crude2net.d(pen_d_crude, pen_bc1_net)
    surv.d.net[(age_bc1+1):MAXAGE] <- 1 - cumsum(pen.d.net[(age_bc1+1):MAXAGE])
    
    # survival from CBC and death from other causes (net)
    surv.all.net <- surv.cbc.net * surv.d.net
    return(surv.all.net)
  }
  
  # loop through genotypes
  for(gt in rownames(pen.cbc.net)){
    
    # subset penetrances for this genotype
    pen.cbc.net.g <- pen.cbc.net[gt,]
    pen.d.crude.g <- pen.d.crude[gt,]
    
    # calculate all cause survival for subject's lifetime
    surv.all.net.cur <- surv.all(pen.cbc.net.g, pen.d.crude.g, 
                                 pen.bc1.net, age.bc1)
    
    # net non-cancer death penetrance vector for this genotype
    pen.d.net <- crude2net.d(pen.d.crude.g, pen.bc1.net)
    
    # calc. probs. of CBC by year from one year after 1st BC to the future age
    prob.cbc.crude.fut <- rep(0, MAXAGE)
    for(i in (age.current+1):MAXAGE){
      prob.cbc.crude.fut[i] <- pen.cbc.net.g[i] * (1 - sum(pen.d.net[1:(i-1)]))
    }
    
    prob.cbc.gt <- cumsum(prob.cbc.crude.fut) / surv.all.net.cur[age.current]
    cbc.frs[gt,] <- prob.cbc.gt
  }
  return(cbc.frs)
}
