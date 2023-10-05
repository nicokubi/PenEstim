#' Prepare model-specific cancer penetrances
#'
#' @param penet A numeric array of cancer penetrances in the same format as 
#' `PanelPRODatabase$Penetrance`. 
#' @param af A numeric matrix of allele frequencies in the same format as 
#' `PanelPRODatabase$AlleleFrequency`, subset to match the genes in the model 
#' specification `MS`. 
#' @param MS A model specification stored in the `MS` component of the return 
#' value of \code{\link{buildDatabase}}. This is a list of genes and cancers 
#' included in the model. 
#' @param cbc.penet A numeric array of contralateral breast cancer penetrances 
#' for carriers of any CBC associated gene in the same format as 
#' `PanelPRODatabase$Contralateral`. `NULL` by default.
#' @details All cancer penetrances, except for CBC penetrances, for carriers are 
#' extracted directly from `penet`. Noncarrier cancer penetrances, except for 
#' CBC, are calculated based on the SEER and carrier penetrances in `penet`, 
#' where the SEER probabilities are taken to be a populational weighted average 
#' of both carriers and noncarriers. CBC carrier penetrances are extracted 
#' directly from `PanelPRODatabase$Contralateral`. CBC penetrances for 
#' noncarriers and carriers of genes not associated with CBC are not included in 
#' the database but are instead estimated as needed in 
#' \code{\link{subsetCancerPenetrance}}.
#' @return A list with three components: 
#' * `penet_c`: A numeric array of carrier cancer penetrances, stratified by 
#' `Cancer`, `Gene`, `Race`, `Sex`, `Age`, and `PenetType`. 
#' * `penet_nc`: A numeric array of noncarrier cancer penetrances, stratified 
#' by `Cancer`, `Gene`, `Race`, `Sex`, `Age`, and `PenetType`. 
#' * `penet_cbc`: A numeric array of CBC carrier penetrances, stratified by 
#' `Gene`, `Sex`, `FirstBCAge`, `PenetType`, and `Age`. 
#' @include helpers.R
.customizePenet <- function(penet, af, MS, cbc.penet = NULL) {

  # Subset carrier penetrances according to the model specification
  penet <- subset_array(penet, c("Cancer", "Gene"),
    list(setdiff(MS$CANCERS,"Contralateral"), c("SEER", MS$ALL_GENE_VARIANTS)),
    drop = F
  )
  
  # Subset contralateral carrier penetrances according to the model spec
  if(!is.null(cbc.penet)){
    cbc.penet <- subset_array(cbc.penet, "Gene", 
                              intersect(MS$ALL_GENE_VARIANTS, 
                                        unlist(CBC_GENE_VARIANT_TYPES)), 
                              drop = F)
  }

  # Ancestry-race combinations
  combs <- as.vector(outer(dimnames(af)$Ancestry, 
                           dimnames(penet)$Race, paste, sep = "."))

  # Initialize container for noncarrier penetrances
  # The "Gene" and "Race" dimension names are misnomers here: "Gene" actually 
  # stores the different ancestry-race combinations, and "Race" seems to be 
  # a bit of a dummy, since there's only a single "Default" slot. However, 
  # fixing this would require carefully checking a bunch of downstream code, 
  # so I'm leaving it for now. 
  dn_penet <- dimnames(penet)
  dn_penet$Gene <- paste0("noncarrier_", combs)
  dn_penet$Race <- "Default"
  penet_NC_all <- array(dim = lengths(dn_penet), dimnames = dn_penet)
  
  # For each ancestry, check if all of the allele frequencies for genes in the 
  # model are the same as nonAJ
  is_all_nonAJ_af = apply(af == af[,"nonAJ"], 2, all)

  # Loop over each ancestry-race combination
  for (i in seq_along(combs)) {
    comb <- unlist(strsplit(combs[i], split = ".", fixed = TRUE))
    anc <- comb[1]
    race <- comb[2]
    
    # Subset allele frequencies by ancestry
    af_sub <- subset_array(af, "Ancestry", anc) 
    names(af_sub) <- dimnames(af)$Gene

    # Probability weights for heterozygous carriers (2pq term)
    weights1 <- 2 * af_sub * (1 - af_sub)
    names(weights1) = sub("_", "_hetero_", names(weights1))

    # Probability weights for homozygous carriers (p^2 term)
    weights2 <- af_sub[names(MS$HOMOZYGOUS)]^2
    if (length(MS$HOMOZYGOUS) > 0) {
      names(weights2) <- sub("_", "_homo_", names(weights2))
    }

    # Put all weights together and reorder to match model specification
    weights <- c(weights1, weights2)[MS$ALL_GENE_VARIANTS]

    # SEER penetrances: populational average of both carriers and non-carriers
    penet_SEER <- subset_array(penet, c("Gene", "Race"),
      list("SEER", race),
      drop = F
    )
    # All carrier penetrances
    penet_sub <- subset_array(penet, c("Gene", "Race"),
      list(-1, race),
      drop = F
    )
    # Calculate noncarrier penetrances
    penet_NC_all[, i, , , , ] <- 
      (abind::adrop(penet_SEER, names(dimnames(penet_SEER)) == "Gene") -
         apply(
           sweep(penet_sub, 2, as.array(weights), "*"),
           (1:length(dim(penet_sub)))[names(dimnames(penet_sub)) != "Gene"], 
           sum)) / (1 - sum(weights))
    
    if (anc != "nonAJ" && is_all_nonAJ_af[anc] == FALSE) {
      penet_NC_all[, i, , , , ] <- 
        abind::adrop(penet_SEER, names(dimnames(penet_SEER)) == "Gene")
    }
  }

  # Make sure penetrances are non-negative
  penet_NC_all[penet_NC_all < 0] <- 0

  # Combine carrier and non-carrier penetrances into a list
  return(list(penet_c = penet, penet_nc = penet_NC_all, penet_cbc = cbc.penet))
}


#' Prepare model-specific risk modifier information
#'
#' @param riskmod A numeric array of risk modifier information in the same 
#' format as `PanelPRODatabase$Riskmod`. 
#' @param MS A model specification stored in the `MS` component of the return 
#' value of \code{\link{buildDatabase}}. This is a list of genes and cancers 
#' included in the model. The default is `NULL`, in which case all information 
#' in `riskmod` will be used. 
#' @return A numeric array of risk modifier information, stratified by 
#' `Cancer`, `Gene`, `Intervention`, `Sex`, `IntervAge`, and `DataType`. 
.customizeRiskmod <- function(riskmod, MS = NULL) {

  # Subset riskmod according to the model specification
  if (!is.null(MS)) {
    riskmod <- subset_array(riskmod, c("Cancer", "Gene"),
      list(MS$CANCERS, c("noncarrier", MS$ALL_GENE_VARIANTS)),
      drop = F
    )
  }

  # Create a "None" risk modifier where all of the entries are ones
  dn_nullarray <- dimnames(riskmod)
  dn_nullarray$Intervention <- "None"
  Null_riskmod <- array(1, dim = lengths(dn_nullarray), dimnames = dn_nullarray)

  # Combine the "None" risk modifier with the other risk modifiers
  riskmod_bind <- abind::abind(Null_riskmod, riskmod, along = 3)
  names(dimnames(riskmod_bind)) <- names(dimnames(riskmod))
  return(riskmod_bind)
}


#' Prepare model-specific death by other causes
#'
#' @param doc A numeric array of probabilities of dying by other causes in the 
#' same format as `PanelPRODatabase$DOC`. 
#' @param af A numeric matrix of allele frequencies in the same format as 
#' `PanelPRODatabase$AlleleFrequency`, subset to match the genes in the model 
#' specification `MS`. 
#' @param MS A model specification stored in the `MS` component of the return 
#' value of \code{\link{buildDatabase}}. This is a list of genes and cancers 
#' included in the model. 
#' @details Probabilities for dying by other causes for carriers are extracted 
#' directly from `doc`. `"nonAJ"` noncarrier probabilities are calculated based 
#' on the SEER and carrier probabilities in `penet` and the allele frequencies 
#' in `af`, where the SEER probabilities are taken to be a populational 
#' weighted average of both carriers and noncarriers. For ancestries other than 
#' `nonAJ` whose allele frequencies differ for one or more genes in the model, 
#' the SEER probabilities are used for the noncarrier probabilities for dying 
#' by other causes.  
#' 
#' @return A list with two components: 
#' * `doc_c`: A numeric array of carrier probabilities of dying by other 
#' causes, stratified by `Cancer`, `Gene`, `Race`, `Sex`, and `Age`. 
#' * `doc_nc`: A numeric array of noncarrier probabilities of dying by other 
#' causes, stratified by `Cancer`, `Gene`, `Race`, `Sex`, and `Age`. 
.customizeDoc <- function(doc, af, MS) {
  
  # Subset carrier probabilities of dying by other causes according to the 
  # model specification
  doc <- subset_array(doc, c("Cancer", "Gene"),
    list(MS$CANCERS, c("SEER", MS$ALL_GENE_VARIANTS)),
    drop = F
  )

  # Ancestry-race combinations
  combs <- as.vector(outer(dimnames(af)$Ancestry, 
                           dimnames(doc)$Race, paste, sep = "."))

  # Initialize container for noncarrier probabilities of dying by other causes
  # The "Gene" and "Race" dimension names are misnomers here: "Gene" actually 
  # stores the different ancestry-race combinations, and "Race" seems to be 
  # a bit of a dummy, since there's only a single "Default" slot. However, 
  # fixing this would require carefully checking a bunch of downstream code, 
  # so I'm leaving it for now. 
  dn_doc <- dimnames(doc)
  dn_doc$Gene <- paste0("noncarrier_", combs)
  dn_doc$Race <- "Default"
  doc_NC_all <- array(dim = lengths(dn_doc), dimnames = dn_doc)
  
  # For each ancestry, check if all of the allele frequencies for genes in the 
  # model are the same as nonAJ
  is_all_nonAJ_af = apply(af == af[,"nonAJ"], 2, all)

  # Loop over each ancestry-race combination
  for (i in seq_along(combs)) {
    comb <- unlist(strsplit(combs[i], split = ".", fixed = TRUE))
    anc <- comb[1]
    race <- comb[2]
    
    # Subset allele frequencies by ancestry
    af_sub <- subset_array(af, "Ancestry", anc) 
    names(af_sub) <- dimnames(af)$Gene

    # Probability weights for heterozygous carriers (2pq term)
    weights1 <- 2 * af_sub * (1 - af_sub)
    names(weights1) = sub("_", "_hetero_", names(weights1))
    
    # Probability weights for homozygous carriers (p^2 term)
    weights2 <- af_sub[names(MS$HOMOZYGOUS)]^2
    if (length(MS$HOMOZYGOUS) > 0) {
      names(weights2) <- sub("_", "_homo_", names(weights2))
    }
    
    # Put all weights together and reorder to match model specification
    weights <- c(weights1, weights2)[MS$ALL_GENE_VARIANTS]

    # SEER penetrances: populational average of both carriers and non-carriers
    doc_SEER <- subset_array(doc, c("Gene", "Race"),
      list("SEER", race),
      drop = F
    )
    # All carrier probabilities of dying by other causes
    doc_sub <- subset_array(doc, c("Gene", "Race"),
      list(-1, race),
      drop = F
    )
    # Calculate noncarrier probabilities of dying by other causes
    doc_NC_all[, i, , , ] <- 
      (abind::adrop(doc_SEER, names(dimnames(doc_SEER)) == "Gene") -
         apply(
           sweep(doc_sub, 2, 1 - as.array(weights), "*"),
           (1:length(dim(doc_sub)))[names(dimnames(doc_sub)) != "Gene"], 
           sum)) / (1 - sum(weights))
    
    if (anc != "nonAJ" && is_all_nonAJ_af[anc] == FALSE) {
      doc_NC_all[, i, , , ] <- 
        abind::adrop(doc_SEER, names(dimnames(doc_SEER)) == "Gene")
    }
  }

  # Make sure death by other causes probabilities are non-negative
  doc_NC_all[doc_NC_all < 0] <- 0

  # Combine carrier and non-carrier death by other causes into a list
  return(list(doc_c = doc, doc_nc = doc_NC_all))
}


#' Build model-specific database
#'
#' Build a database of parameters for only the genes and cancers in the model. 
#' The genes and cancers are taken from the model specification `model_spec` 
#' or user-supplied vectors of `genes` and `cancers`. Allele frequencies, risk 
#' modifiers, and the net cancer/death by other causes penetrances for 
#' carriers are subsetted directly from the master database, `ppd`. The crude 
#' cancer penetrances for carriers are calculated from the net cancer 
#' penetrances for carriers and the penetrances for death by other causes. The 
#' noncarrier cancer and death by other causes penetrances are calculated from 
#' the carrier penetrances and the SEER penetrances. 
#'
#' @param model_spec One of the following character strings indicating a 
#' pre-specified model: `"BRCAPRO"`, `"BRCAPRO5"`, `"BRCAPRO6"`, `"MMRPRO"`, 
#' `"PanPRO11"`, and `"PanPRO22"`. See `PPP:::MODELPARAMS` for more 
#' information. Either `model_spec` or `genes` and `cancers` should be 
#' specified by the user. When both input types are supplied, the input from 
#' `genes` and `cancers` will take precedence. The default is `NULL`, for no 
#' model specified. 
#' @param ppd A structured list providing the parameters of the model. The 
#' default is \code{\link{PanelPRODatabase}}. Users who wish to provide their 
#' own model parameters should supply a list with the same structure as 
#' `PanelPRODatabase`. 
#' @param genes A character vector of the genes of interest. The default is 
#' `NULL`. Available options are listed in `PPP:::GENE_TYPES`. Either 
#' `model_spec` or `genes` and `cancers` should be specified by the user. When 
#' both input types are supplied, the input from `genes` and `cancers` will 
#' take precedence. It is considered valid to specify no genes in `genes` as 
#' long as at least one cancer is specified in `cancers`. 
#' @param cancers A character vector of the cancers of interest. Available 
#' options are listed in `PPP:::CANCER_TYPES`. The default is `NULL`. 
#' Either `model_spec` or `genes` and `cancers` should be specified by the 
#' user. When both input types are supplied, the input from `genes` and 
#' `cancers` will take precedence. It is considered valid to specify no cancers 
#' in `cancers` as long as at least one gene is specified in `genes`. 
#' 
#' @return A list of 7 elements containing the model parameters:
#' 
#' * `af`: Allele frequencies for the genes. Can be subsetted by ethnicity. 
#' * `penet`: Cancer penetrances, split into `penet_c` for mutation carriers,
#'  `penet_nc` for noncarriers and `penet_cbc` for contralateral breast cancer
#'  penetrances for carriers and noncarriers. `penet_c` and `penet_nc` can be 
#'  subset by cancer, gene, race, sex, age, and type of penetrance (net or 
#'  crude). `penet_cbc` is a list of two arrays, `cbc_noncarrier` for 
#'  female non-carriers and `cbc_carrier` for carriers of any sex. The CBC 
#'  non-carrier array can be subsetted by gene, age of first breast cancer, type 
#'  of first breast cancer, if first degree relatives had breast cancer, ER 
#'  result, age category at first child birth, if anti-estrogen therapy was used 
#'  to treat the first breast cancer, history of high risk preneoplasia, breast 
#'  density, net or crude penetrance, and age. The CBC carrier array can be 
#'  subsetted by gene, sex, age of first breast cancer, net or crude penetrance, 
#'  and age. 
#' * `riskmod`: Information for interventions, such as mastectomies, 
#' oophorectomies, and hysterectomies, that impact the cancer penetrances. Can 
#' be subsetted by cancer, gene, intervention, sex, age of intervention, and 
#' data type (relative risk or hazard ratio).
#' * `doc`: Probability of dying of causes other than the cancer of interest, 
#' split into `doc_c` for mutation carriers and `doc_nc` for noncarriers. Each 
#' can be subsetted by cancer, gene, race, sex, and age.
#' * `germline`: Sensitivity and specificity for germline testing, which can be 
#' used to modify the likelihood. 
#' * `biomarker`: Information on tumor biomarker testing, which can be used to 
#' modify the likelihood.  
#' * `MS`: List of genes and cancers specified by the model.
#'
#' @examples
#' brcadb1 <- buildDatabase(model_spec = "BRCAPRO")
#' brcadb2 <- buildDatabase(genes = c("BRCA1", "BRCA2"), 
#'                          cancers = c("Breast", "Ovarian"))
#' stopifnot(all.equal(brcadb1, brcadb2))
#' @export
buildDatabase <- function(model_spec = NULL, ppd = PPP::PanelPRODatabase, 
                          genes = NULL, cancers = NULL) {

  # Attempt to map the model specification to genes and cancers
  if (!is.null(model_spec) && is.null(genes) && is.null(cancers)) {
    model_spec <- match.arg(model_spec, names(MODELPARAMS))
    genes <- MODELPARAMS[[model_spec]]$GENES
    cancers <- MODELPARAMS[[model_spec]]$CANCERS
  } 
  
  # Check if we have something for genes or cancers
  if (length(genes) + length(cancers) == 0) {
    rlang::abort("Gene and cancer sets cannot both be empty. 
                 Please specifiy genes and/or cancers explicitly, or use predefined model.")
  }
  
  # Capitalize all letters for matching
  genes <- toupper(genes)
  if (length(cancers) > 0) {
    cancers <- tools::toTitleCase(cancers)
  }
  
  # Check if duplicate entries exist
  if (any(duplicated(genes))) {
    rlang::inform("Genes input contain duplicates, using unique version.",
                level = "ModelSpecDuplicated"
    )
    genes <- unique(genes)
  }
  if (any(duplicated(cancers))) {
    rlang::inform("Cancers input contain duplicates, using unique version.",
                level = "ModelSpecDuplicated"
    )
    cancers <- unique(cancers)
  }
  
  # Check if the provided names of are in the list of supported gene and cancer types
  no_match_genes <- is.na(pmatch(genes, GENE_TYPES))
  no_match_cancers <- is.na(pmatch(cancers, CANCER_TYPES))
  if (any(no_match_genes)) {
    genes_not_matched <- genes[no_match_genes]
    rlang::inform(sprintf("The genes %s do not match with existing genes, deleting the unmappables.",
                          paste0(genes_not_matched, collapse = ", ")))
    genes <- genes[!no_match_genes]
  }
  if (any(no_match_cancers)) {
    cancers_not_matched <- cancers[no_match_cancers]
    rlang::inform(sprintf("The cancers %s do not match with existing cancers, deleting the unmappables.",
                          paste0(cancers_not_matched, collapse = ", ")))
    cancers <- cancers[!no_match_cancers]
  }
  
  # Check if no genes or cancers were found. 
  if (length(genes) + length(cancers) == 0) rlang::abort("No cancers or genes found.",
                                                         level = "EmptyCancersGenes")
    
    # Identify all gene variants that map to genes in model specification
    all_gene_variants = unlist(ALL_GENE_VARIANT_TYPES[intersect(genes, names(ALL_GENE_VARIANT_TYPES))])
    alleles = unique(sub("_.*_", "_", all_gene_variants))
    # Identify genes with homozygous mutations
    homozygous_genes = HOMOZYGOUS_VARIANT_TYPES[intersect(alleles, names(HOMOZYGOUS_VARIANT_TYPES))]
    
  # Model specification
  MS <- list(CANCERS = cancers, GENES = genes, 
             ALL_GENE_VARIANTS = all_gene_variants, 
             ALLELES = alleles, 
             HOMOZYGOUS = homozygous_genes)
  
  # Message on model specification
  rlang::inform(sprintf(
    "Your model has %s cancers - %s and %s genes - %s.",
    length(MS$CANCERS),
    paste0(MS$CANCERS, collapse = ", "),
    length(MS$GENES),
    paste0(MS$ALL_GENE_VARIANTS, collapse = ", ")
  ), level = "ModelSpecification")
  
  
  # Allele frequencies
  af <- subset_array(ppd$AlleleFrequency, "Gene", 
                     alleles, drop = FALSE)
  
  # Cancer penetrances, options if model included CBC or not
  if(!"Contralateral" %in% cancers){
    penet <- .customizePenet(ppd$Penetrance, af, MS)
  } else {
    penet <- .customizePenet(ppd$Penetrance, af, MS, ppd$Contralateral)
  }
  
  # Risk modifiers 
  riskmod <- .customizeRiskmod(ppd$Riskmod, MS)

  # Death by other causes 
  doc <- .customizeDoc(ppd$DOC, af, MS)

  # Likelihood modifiers
  germline <- ppd$GermlineTesting[MS$ALL_GENE_VARIANTS, , drop = FALSE]
  biomarker <- ppd$BiomarkerTesting[intersect(MS$CANCERS, names(ppd$BiomarkerTesting))]

  # Gather output as a list 
  db <- list(
    af = af, penet = penet, riskmod = riskmod, doc = doc,
    germline = germline, biomarker = biomarker,
    MS = MS
  )
  
  return(db)
}

