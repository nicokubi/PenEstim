#' Run peeling and paring algorithm
#'
#' Returns the posterior genotype distribution for the proband(s), based on 
#' the peeling and paring algorithm. Uses a C++ routine to speed up the main 
#' algorithm.
#' 
#' @param ped A checked pedigree data frame returned by \code{\link{checkFam}}. 
#' @param proband A numeric value or vector of the unique IDs in `ped` for 
#' whom to return posterior probabilities.
#' @param af A matrix of allele frequencies, where each row corresponds to the 
#' allele frequencies of a relative in `ped`. 
#' @param LIK A likelihood matrix returned by \code{\link{calcLik}}, possibly 
#' modified by the return values of \code{\link{germlineContrib}} and/or 
#' \code{\link{markerContrib}}. 
#' @param K The total number of mutations. 
#' @param M The maximum number of simultaneous mutations allowed. The default 
#' is `2`.  
#' @param homo_idx List of indices for genes with heterozygous and homozygous 
#' mutations. Each element consists of a vector with the index pair 
#' corresponding to the heterozygous and homozygous mutations for a given gene. 
#' The default is an empty list (for no genes). 
#'
#' @family peeling 
#' @return Matrix of carrier probabilities where the rows represent the 
#' relatives with IDs in `proband` and the columns represent the genotypes. 
#' @md
pp.peelingParing <- function(ped, proband, af, LIK, K, M = 2, 
                             homo_idx = list()) {

  # Get possible genotypes
  PG <- .ppPossibleGenotypes(K, M, homo_idx = homo_idx, collapse = TRUE)
  
  # Extract the total number of admissible genotypes under the approximation
  U <- nrow(PG)

  # Edge matrix is a binary matrix that indicates which individuals are mated
  E <- .ppEdgeMatrix(ped) # based on rows of ped
  C <- .ppChildMatrix(ped) # based on rows of ped

  # Calculate the prevalence of each genotype
  prev <- .calPrevalences(af, PG) # CPP version
  
  # Calculate the probability of inheriting each genotype based on Mendelian 
  # genetics 
  # See .calInheritanceProb documentation for examples/interpretation of 
  # inheritance probabilities
  TR <- .calInheritanceProb(PG) # CPP version
  
  # Initialise antProb and postProb for storing anterior and posterior 
  # probabilities
  antProb <- .ppInitializeA(PG, ped, prev)
  antProb[is.na(antProb)] <- 0 # Get rid of NAs
  postProb <- .ppInitializeP(PG, E)
  postProb[is.na(postProb)] <- 0 # Get rid of NAs

  # Get the indices of the probands to pass to C++
  probandIndexes <- c()
  for (i in 1:length(proband)) {
    probandIndexes <- c(probandIndexes, which(ped$ID == proband[i]) - 1)
  }

  # Extract the ID, MotherID and FatherID columns
  # This will save time and not copy the other columns of ped
  idMatrix <- as.matrix(ped[, c("ID", "MotherID", "FatherID")])

  # Normalise LIK by row to sum to 1 for accuracy
  LIK <- LIK / rowSums(LIK)
  
  # Call CPP routine
  posterior_probs <- .peelingParing(probandIndexes, idMatrix, E, C, prev, LIK,
                                    TR, antProb, postProb)

  # Column names are the genotypes
  colnames(posterior_probs) <- colnames(LIK)
  # Row names are the proband IDs
  rownames(posterior_probs) <- proband
  
  return(posterior_probs)
}
