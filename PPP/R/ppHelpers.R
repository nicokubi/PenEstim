#' Create a 3D array that links each individual to their parents
#'
#' @param ped A pedigree data frame. 
#' @return A 3D numeric array where each dimension is of the same length as the 
#' number of individuals in `ped`. Entry `ijk = 1` if individual `k` is the 
#' offspring of individuals `i` and `j`, `0` otherwise. If only one parent for 
#' a given child is in `ped`, they are assigned as both of the child's parents. 
#' Parents are considered missing from the pedigree when `MotherID` or 
#' `FatherID` is `-999`. 
#' @family peeling
.ppChildMatrix <- function(ped) {
  # Initialize 3D array of zeros
  C <- array(0, rep(nrow(ped), 3))
  
  # Iterate through all the individuals in ped
  for (i in 1:nrow(ped)) {
    # If both parents of individual i are in the pedigree, assign them as the
    # parents
    if (ped$FatherID[i] != -999 & ped$MotherID[i] != -999) {
      C[ped$ID == ped$MotherID[i], ped$ID == ped$FatherID[i], i] <- 1
      C[ped$ID == ped$FatherID[i], ped$ID == ped$MotherID[i], i] <- 1
      # If only one parent is in the pedigree, assign them as both parents
    } else if (ped$FatherID[i] != -999) {
      C[ped$ID == ped$FatherID[i], ped$ID == ped$FatherID[i], i] <- 1
    } else if (ped$MotherID[i] != -999) {
      C[ped$ID == ped$MotherID[i], ped$ID == ped$MotherID[i], i] <- 1
    }
  }
  return(C)
}


#' Create a matrix that links each individual to their mate(s)
#'
#' @param ped A pedigree data frame. 
#' @return A numeric matrix where each dimension is of the same length as the 
#' number of individuals in `ped`. Entry `ij = 1` if individual `i` is mated 
#' with individual `j`, `0` otherwise. If only one parent for a given child is 
#' in `ped`, the parent is assigned as their own mate. Parents are considered 
#' missing from the pedigree when `MotherID` or `FatherID` is `-999`. 
#' @family peeling
.ppEdgeMatrix <- function(ped) {
  # Initialize matrix of zeros
  E <- matrix(0, nrow = nrow(ped), ncol = nrow(ped))
  
  # Iterate through all the individuals in ped
  for (i in 1:nrow(ped)) {
    # If both parents of individual i are in the pedigree, assign them as each
    # other's mates.
    if (ped$FatherID[i] != -999 & ped$MotherID[i] != -999) {
      E[ped$ID == ped$FatherID[i], ped$ID == ped$MotherID[i]] <- 1
      E[ped$ID == ped$MotherID[i], ped$ID == ped$FatherID[i]] <- 1
      # If only one parent is in the pedigree, assign them as their own mate
    } else if (ped$FatherID[i] != -999) {
      E[ped$ID == ped$FatherID[i], ped$ID == ped$FatherID[i]] <- 1
    } else if (ped$MotherID[i] != -999) {
      E[ped$ID == ped$MotherID[i], ped$ID == ped$MotherID[i]] <- 1
    }
  }
  return(E)
}


#' Initialize matrix to store anterior information
#'
#' @param PG A matrix of possible genotypes returned by 
#' \code{\link{.ppPossibleGenotypes}}. 
#' @param ped A pedigree data frame. 
#' @param prev A matrix of genotype prevalences, where each row corresponds to 
#' an individual in `ped` and each column corresponds to a genotype in `PG`, 
#' returned by \code{\link{.calPrevalences}}. 
#' @return A numeric matrix where the rows correspond to individuals in `ped` 
#' and the columns correspond to the genotypes in `PG`. Entry `ij` stores the 
#' anterior probability for the `i`th family member and `j`th genotype. At 
#' initialization, this matrix is populated with the genotype prevalences in 
#' `prev` for each individual, i.e. the "prior" probabilities. 
#' @family peeling
.ppInitializeA <- function(PG, ped, prev) {
  # Initialize matrix of NAs
  A <- matrix(NA, nrow = nrow(ped), ncol = nrow(PG))
  
  # Populate A with the genotype prevalences in `prev` for each individual
  for (i in 1:nrow(ped)) {
    if (ped$FatherID[i] == -999 & ped$MotherID[i] == -999) {
      A[i, ] <- prev[i, ]
    }
  }
  return(A)
}


#' Initialize 3D array to store posterior information
#'
#' @param PG A matrix of possible genotypes returned by 
#' \code{\link{.ppPossibleGenotypes}}. 
#' @param E An edge matrix returned by \code{\link{.ppEdgeMatrix}}. 
#' @return A 3D numeric array where the first two dimensions are of the same 
#' length as the number of individuals in the pedigree (equivalent to the 
#' lengths of the dimensions of `E`) and the third dimension is of the same 
#' length as the number of genotypes in `PG`. Entry `ijk` stores the posterior 
#' fprobability or the `i`th and `j`th family members and the `k`th genotype in 
#' `PG`. At initialization, entries where individuals `i` and `j` are unmated 
#' are populated with `1`s. 
#' @family peeling
.ppInitializeP <- function(PG, E) {
  # Initialize 3D array of NAs
  P <- array(NA, c(nrow(E), nrow(E), nrow(PG)))
  
  # Populate P with 1s for unmated pairs
  for (i in 1:nrow(E)) {
    for (j in 1:i) {
      if (E[i, j] == 0) {
        P[i, j, ] <- 1
      } else if (E[j, i] == 0) {
        P[j, i, ] <- 1
      }
    }
  }
  return(P)
}


#' Get possible genotypes
#'
#' @param K The total number of mutations. 
#' @param M The maximum number of simultaneous mutations allowed. 
#' @param homo_idx List of indices for genes with heterozygous and homozygous 
#' mutations. Each element consists of a vector with the index pair 
#' corresponding to the heterozygous and homozygous mutations for a given gene. 
#' The default is an empty list (for no genes).
#' @param collapse A logical value indicating whether or not to collapse the
#' columns for genes with heterozygous and homozygous mutations into a single
#' column, where heterozygous is coded as `1` and homozygous is coded as `2`.
#' The default is `FALSE` for no collapsing.
#' @return A numeric matrix of possible genotypes based on `K` mutations 
#' (potentially distinguishing between homozygous and heterozygous mutations), 
#' restricted to at most `M` simultaneous mutations. Each row represents a 
#' genotype and each column represents a gene variant. If `collapse = FALSE`, 
#' the matrix is comprised of `0` and `1` values that encode the genotypes. If 
#' `collapse = TRUE`, columns for heterozygous and homozygous mutations of the 
#' same gene variant the matrix are collapsed such that heterozygous is coded 
#' as `1` and homozygous is coded as `2`. 
#' @family peeling
#' @seealso \code{\link{.getPossibleGenotype}}
.ppPossibleGenotypes <- function(K, M, homo_idx = list(), collapse = FALSE) {
  # Gene number for each gene
  G <- array(0, c(2, K))
  G[, 1] <- c(0, 1)

  if (K > 1) {
    for (i in 2:K) {
      # All the genotypes with M-1 or fewer mutations before locus i
      GG <- as.matrix(G[(G %*% rep(1, K)) < M, ])
      # If there is only one resulting genotype, transpose it
      if ((dim(GG)[2]) == 1) {
        GG <- t(GG)
      }
      # For these genotypes, it is possible to have a mutation at locus i under 
      # the restriction
      GG[, i] <- 1

      # Collect the genotypes together
      G <- rbind(G, GG)
      # There's an implication here that you can have at most two mutations at 
      # each locus: each allele is either WT(=0) or mutant(>0); the numbering 
      # of mutant types is non-informative. 
    }
  }
  
  # Translate the genotypes to the genotype matrix PG
  PG <- c()
  total_muts <- apply(G, 1, sum)
  for (i in 0:M) {
    PG <- rbind(PG, G[total_muts == i, ])
  }
  
  # Remove genotypes that include both heterozygous and homozygous mutations 
  # for the same gene variant
  if (length(homo_idx) > 0) {
    for (i in 1:length(homo_idx)) {
      drop_idx <- which(rowSums(PG[, homo_idx[[i]]]) > 1)
      if (length(drop_idx) > 0) {
        PG <- PG[-drop_idx, ]
      }
    }
  }

  # If requested, collapse the columns for gene variants with both heterozygous 
  # and homozygous mutations into a single column, where heterozygous is coded 
  # as 1 and homozygous is coded as 2
  if (collapse == TRUE && length(homo_idx) > 0) {
    # Re-code homozygous mutations as 2
    PG[, sapply(homo_idx, function(x){x[2]})] <- 
      PG[, sapply(homo_idx, function(x){x[2]})] * 2

    # Temporary matrix for storing collapsed version of PG
    PG_temp <- c()
    # Iterate through columns of PG
    for (i in 1:ncol(PG)) {
      if (!(i %in% unlist(homo_idx))) {
        # Include genes without homozygous mutations in PG
        PG_temp <- cbind(PG_temp, PG[, i])
      } else if (i %in% sapply(homo_idx, function(x){x[1]})) {
        # For genes with homozygous mutations, use the sum of the two columns
        PG_temp <- cbind(
          PG_temp,
          rowSums(PG[, homo_idx[[which(i == sapply(homo_idx, 
                                                   function(x){x[1]}))]]])
        )
      }
    }

    # Replace PG with the collapsed version
    PG <- PG_temp
  }

  return(PG)
}
