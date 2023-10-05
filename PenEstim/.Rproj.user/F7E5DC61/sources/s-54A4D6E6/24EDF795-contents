#' Modify likelihood using germline testing results
#'
#' @param fam A checked pedigree data frame returned by \code{\link{checkFam}}. 
#' @param db A model-specific database returned by \code{\link{buildDatabase}}. 
#' @param PGs Possible genotypes in both list and data frame format, returned 
#' by `.getPossibleGenotype`. 
#' @return Germline testing contribution matrix where the rows represent the 
#' relatives in the pedigree and the columns represent all possible genotypes 
#' in `PGs`. 
#' @family modifylik
germlineContrib <- function(fam, db, PGs) {
  
  # Extract gene and germline information from database
  MS_all_gene_variants <- db$MS$ALL_GENE_VARIANTS
  GERMLINE <- db$germline
  
  # Family size
  fam_size <- nrow(fam)
  # Pedigree columns for germline testing
  gt <- fam[, intersect(MS_all_gene_variants, colnames(fam)), drop = FALSE]
  
  # Intersect of specified genes and genes in the pedigree
  genes_to_consider <- intersect(MS_all_gene_variants, colnames(gt))

  # Positive and negative test results with NAs forced to be FALSE
  pos.test <- .forceFalse(gt == 1)
  neg.test <- .forceFalse(gt == 0)

  # Helper function for likelihood modification per gene 
  mod_per_gene <- function(gene) {
    TES <- matrix(1, ncol = length(PGs$list), nrow = fam_size)
    contains <- PGs$df[[gene]] == 1
    tpos <- pos.test[, gene]
    tneg <- neg.test[, gene]

    # Assign test sensitivity and specificity for each gene
    TES[tneg, !contains] <- GERMLINE[gene, "specificity"]
    TES[tpos, contains] <- GERMLINE[gene, "sensitivity"]
    TES[tpos, !contains] <- 1 - GERMLINE[gene, "specificity"]
    TES[tneg, contains] <- 1 - GERMLINE[gene, "sensitivity"]
    return(TES)
  }
  
  # Take product over all sensitivity/specificity of testing results for each 
  # relative
  TES_all <- Reduce("*", lapply(genes_to_consider, mod_per_gene))
  
  return(TES_all)
}


#' Modify likelihood using tumor biomarker testing results
#'
#' @param fam A checked pedigree data frame returned by \code{\link{checkFam}}. 
#' @param db A model-specific database returned by \code{\link{buildDatabase}}. 
#' @param PGs Possible genotypes in both list and data frame format, returned 
#' by `.getPossibleGenotype`. 
#' @return Tumor biomarker testing contribution matrix where the rows 
#' represent the relatives in the pedigree and the columns represent all 
#' possible genotypes in `PGs`. If no marker testing information is available 
#' for the cancers/genes in the model, return 1. 
#' @family modifylik
markerContrib <- function(fam, db, PGs) {

  # Extract gene and cancer information from database
  MS_all_gene_variants <- db$MS$ALL_GENE_VARIANTS
  MS_cancers = db$MS$CANCERS

  # Construct genotype matrix
  psize <- nrow(fam)
  geno_mat <- matrix(1, nrow = psize, ncol = length(PGs$list))

  # Map cancer names from marker testing information
  mt_cancers = .mapCancerNames(short = names(MARKER_TESTING))

  # If none of the cancers in the model have marker testing information, 
  # return 1
  if (!any(MS_cancers %in% mt_cancers)) {
    return(1)
  }

  # Helper function for likelihood modification per cancer
  per_cancer <- function(cancer, possible_genes, possible_markers) {

    # Test if genes and markers in pedigree are valid
    gene_set <- intersect(possible_genes, MS_all_gene_variants)
    marker_set <- intersect(possible_markers, colnames(fam))

    # If none of the genes in the model have marker testing information, 
    # return 1
    if (length(gene_set) == 0 | length(marker_set) == 0) {
      return(1)
    }

    # Assign test sensitivity and specificity for each marker
    ped_markers <- matrix(-999, nrow = psize, ncol = length(possible_markers))
    colnames(ped_markers) <- possible_markers
    ped_markers[, marker_set] <- as.matrix(fam[, marker_set, drop = FALSE])
    ped_markers[ped_markers == 1] <- "Pos"
    ped_markers[ped_markers == 0] <- "Neg"
    ped_markers[ped_markers == -999] <- "NoTest"

    # Assign marker testing based on genotype
    genotype_info <- stats::setNames(rep(0, length(possible_genes)), possible_genes)
    M <- geno_mat
    for (i in 1:nrow(M)) {
      patient_info <- unlist(ped_markers[i, ])
      for (j in 1:ncol(M)) {
        genotype <- unlist(PGs$df[j, gene_set])
        names(genotype) <- gene_set
        genotype_info[gene_set] <- genotype
        M[i, j] <- db$biomarker[[cancer]][matrix(c(patient_info, genotype_info), 
                                                 nrow = 1, byrow = TRUE)]
      }
    }
    return(M)
  }
  
  # Modify likelihood using tumor markers, if the associated cancer is in the 
  # model
  likm = lapply(1:length(mt_cancers), function(i) {
    if (all(MS_cancers != mt_cancers[i])) {
      return(1)
    } else {
      return(likm <- per_cancer(
        mt_cancers[i], 
        MARKER_TESTING[[i]][["GENES"]],
        MARKER_TESTING[[i]][["MARKERS"]]
      ))
    }
  })
  
  return(Reduce("*", likm))
}


#' Modify likelihood for identical twins/multiple births
#'
#' The modified likelihood matrix replaces the likelihood rows corresponding to 
#' members of a given twin set with the element-wise product of these rows, 
#' while only keeping one member of the twin set in the pedigree. 
#' 
#' @param lik A likelihood matrix returned by \code{\link{calcLik}}. 
#' @param ped A checked pedigree data frame returned by \code{\link{checkFam}}. 
#' @param proband A numeric value or vector of the unique IDs in `ped` for whom 
#' to return posterior probabilities and future risks. 
#' @return A list with four components: 
#' * `lik`: Modified likelihood matrix. 
#' * `ped`: Collapsed pedigree data frame. 
#' * `proband`: Collapsed vector of probands. 
#' * `twin_labels_df`: Data frame for keeping track of twin information with 4
#' columns (`Label`, `ID`, `isKept`, and `isProband`). 
#' @family modifylik
.twinsLikMod <- function(lik, ped, proband) {
  # Initialize data frame for storing information on twins
  twin_labels_df <- data.frame()
  
  # Iterate through each set of twins
  for (label in unique(ped$Twins[ped$Twins != 0])) {
    twin_ids <- ped$ID[ped$Twins == label]
    if (length(twin_ids) == 1) {
      next
    }
    
    # Multiply twin likelihoods together
    twin_lik_prod <- apply(lik[ped$ID %in% twin_ids, ], 2, prod)
    
    # If twins includes a proband, re-order the twin_ids so the proband
    # comes first and is therefore kept in the list of proband IDs
    if (any(twin_ids %in% proband)) {
      twin_ids <- c(
        twin_ids[twin_ids %in% proband],
        twin_ids[!(twin_ids %in% proband)]
      )
    }
    # First twin (will be keeping)
    first_twin <- twin_ids[1]
    # Other twins (will be dropped)
    other_twins <- twin_ids[-1]
    
    # old proband and keep track of which twin pair they belong to
    twin_labels_df <- rbind(
      twin_labels_df,
      data.frame(
        Label = label,
        ID = twin_ids,
        isKept = c(1, rep(0, length(other_twins))),
        isProband = ifelse(twin_ids %in% proband, 1, 0)
      )
    )
    
    # Re-assign parent IDs of the offspring of the other twins
    ped$MotherID[ped$MotherID %in% other_twins] <- first_twin
    ped$FatherID[ped$FatherID %in% other_twins] <- first_twin
    
    # Drop the other twins from the likelihood, pedigree, and probands
    lik <- lik[!(ped$ID %in% other_twins), , drop = FALSE]
    ped <- ped[!(ped$ID %in% other_twins), , drop = FALSE]
    proband <- proband[!(proband %in% other_twins)]
    
    # Assign new twin likelihood
    lik[ped$ID == first_twin, ] <- twin_lik_prod
  }
  
  return(list(
    lik = lik, ped = ped, proband = proband,
    twin_labels_df = twin_labels_df
  ))
}
