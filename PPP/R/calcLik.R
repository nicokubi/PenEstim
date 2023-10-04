#' Calculate likelihood
#'
#' Calculates the likelihood of the observed cancer history of each individual 
#' in the pedigree, given each possible genotype, based on the cancers and 
#' genes in the model.
#'
#' @param fam A checked pedigree data frame returned by \code{\link{checkFam}}. 
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
#' @return Likelihood matrix where the rows represent the relatives in `fam` 
#' and the columns represent the possible genotypes in `PGs`.  
#' @include helpers.R
calcLik <- function(fam, sub_dens, PGs, 
                    direct_fill_PGs, multi_PGs, multi_muts) {
  
  # Get cancers from pedigree
  cancer_cols <- .getCancersFromFam(fam)
  long_cancers = .mapCancerNames(short = cancer_cols)
  
  # Identify cancer affection and age columns
  isAff <- fam[, paste0("isAff", cancer_cols), drop = FALSE]
  Age <- as.matrix(fam[, paste0("Age", cancer_cols), drop = FALSE])
  
  # Calculate likelihood for each family member
  lik = t(sapply(1:nrow(fam), function(i) { 
    
    # Extract ID
    ID = fam$ID[i]
    
    # Extract current age; if missing, return a flat likelihood
    curage = fam$CurAge[i]
    if (curage == -999) {
      return(1)
    }
    
    # Identify column indices for cancers for which the person is affected and
    # and unaffected
    aff_idx <- which(isAff[i, ] == 1) # Affected
    unaff_idx <- which(isAff[i, ] == 0) # Unaffected
    
    # CBC should only be considered for subject survival if 1st BC already occurred
    if(all(c("BC", "CBC") %in% cancer_cols)){
      if(isAff$isAffBC[i] == 0){
        unaff_idx <- setdiff(unaff_idx, which(colnames(isAff) == "isAffCBC"))
      }
    }
    
    # Likelihood is the product of the survivals for all the cancers that the 
    # person didn't get and the densities for all the cancers they did get
    ccs <- calcCancerSurv(ID, curage, long_cancers[unaff_idx], 
                          sub_dens, PGs, 
                          direct_fill_PGs, multi_PGs, multi_muts)
    ccd <- calcCancerDens(ID, Age[i, aff_idx], long_cancers[aff_idx], 
                          sub_dens, PGs, 
                          direct_fill_PGs, multi_PGs, multi_muts)
    tmp.lik <- ccs * ccd
        
    return(tmp.lik)
    
  }))
  
  return(lik)
}
