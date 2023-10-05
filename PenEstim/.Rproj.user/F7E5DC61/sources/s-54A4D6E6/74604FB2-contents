#' Apply function that can safely accept a vector
#' 
#' @param arr An array of arbitrary dimensions. Can be a vector.
#' @param mar A vector giving the subscripts which the function will be applied 
#' over. Ignored when `arr` is a vector. 
#' @param FUN The function to be applied.
#' @param ... Optional arguments to FUN.
#' @return 
#' When `arr` is a vector, returns `FUN(arr, ...)`. When `arr` is an array with 
#' 2 or more dimensions, returns `apply(arr, mar, FUN, ...)`. 
safe_apply <- function(arr, mar, FUN, ...) {
  if (is.null(dim(arr))) {
    return(FUN(arr, ...))
  }
  else {
    return(apply(arr, mar, FUN, ...))
  }
}


#' Replace values with another value (`NA` by default)
#' 
#' @details 
#' The expected usage is in the context of an apply function. Note that 
#' this function also replace negative values, in addition to `na_val`. 
#' @param na_val Value to replace. 
#' @param rep Replacement value. The default is `NA`. 
#' @return A function that takes a single parameter, `x`. `x` should be 
#' specified as a vector that is compatible with `na_val` and `rep`. Evaluating 
#' the return function will return `x` with values of `na_val` replaced by 
#' `rep`. 
fix_with_rep <- function(na_val, rep = NA) {
  function(x) {
    x[x %in% na_val] <- rep
    x[x < 0] <- rep
    x
  }
}


#' Extract the cancers from the column names of a pedigree data frame
#' 
#' @param fam A pedigree data frame. 
#' @return Character vector of short (abbreviated) names of cancers extracted 
#' directly from the `isAFF*` columns of a pedigree. Does not ensure that a 
#' corresponding `Age*` column exists. 
.getCancersFromFam <- function(fam) {
  cn <- dimnames(fam)[[2]]
  cancers <- regmatches(cn, regexpr("(?<=^isAff).*", cn, perl = T))
  cancers <- cancers[!cancers %in% "Any"]
  return(cancers)
}


#' Subset array by a dimension of choice and optionally reduce that dimension
#' 
#' @param arrayInput A multi-dimensional array. 
#' @param axisName A character string giving one or more named axis dimensions 
#' of `arrayInput` of which to index. 
#' @param choice A list of indices to subset for each dimension. 
#' @param ... Additional arguments passed in to `asub`. 
#' @return The array which has been subsetted
subset_array <- function(arrayInput, axisName, choice, ...) {
  dimNames <- dimnames(arrayInput)
  axisId <- which(names(dimnames(arrayInput)) %in% axisName)
  abind::asub(arrayInput, choice, axisId, ...)
}


#' Get possible genotypes
#' 
#' Generates all possible genotypes based on a list of gene variants, 
#' restricted to a maximum number of mutations. 
#' @param gene A character vector of gene names
#' @param max_mut The maximum number of simultaneous mutations allowed. 
#' @param homo_genes A list of genes with heterozygous and homozygous mutations. 
#' Each component is named after a gene and contains a character vector of the 
#' full gene names for that gene. The default is an empty list `list()`, which 
#' indicates that no genes have both heterozygous and homozygous mutations. 
#' @return A list with two components: 
#' - `df`: A data frame of possible genotypes for the genes in `gene`, 
#' restricted to at most `max_mut` simultaneous mutations. Each row represents 
#' a genotype and each column represents a gene variant. The data frame is 
#' comprised of `0` and `1` values that encode the genotypes. 
#' - `list`: A character vector of genotype names corresponding to the rows of 
#' `df`. 
#' @seealso \code{\link{.ppPossibleGenotypes}}
.getPossibleGenotype <- function(gene, max_mut = 2, 
                                 homo_genes = list()) {
  
  # Indices for gene with heterozygous and homozygous mutations
  if (length(homo_genes) == 0) {
    # Empty list if there are no genes with heterozygous and homozygous 
    # mutations
    homo_idx <- list()
  } else {
    # Otherwise, identify the index pairs for heterozygous and homozygous 
    # mutations
    homo_idx <- lapply(homo_genes, function(g) {
      which(gene %in% g)
    })
  }
  
  # Get data frame of possible genotypes
  genodf <- .ppPossibleGenotypes(length(gene), max_mut, 
                                 homo_idx, collapse = FALSE)
  colnames(genodf) <- gene

  # Generate corresponding names for possible genotypes
  gps <- apply(genodf, 1, 
               function(v) paste0(colnames(genodf)[as.logical(v)], 
                                  collapse = "."))
  gps[1] <- "noncarrier"
  
  return(list(df = as.data.frame(genodf), list = gps))
}


#' Drop array dimensions that have only one level
#' 
#' @param arr An array. 
#' @return 
#' An array or vector with reduced dimensions from dropping dimensions with 
#' only one level. 
#' @keyword internal
adrop_one <- function(arr) {
  abind::adrop(arr, drop = which(dim(arr) == 1))
}


#' Ensure that pedigree data frame column only contains allowed values
#' 
#' Used to check categorical and logical columns in \code{\link{checkFam}}. 
#' @param ped A pedigree data frame. 
#' @param col_spec The name of a column if `ped` with a specific attribute. 
#' @param allowed A vector of the values that are allowed in `ped[[col_spec]]`. 
#' @param default_NA_rep The value to default to if `NA`s are encountered. The 
#' default is `NULL`. 
#' @param NA_Ok A logical value indicating whether `NA`s can be present in 
#' `ped[[col_spec]]`. The default is `TRUE`. 
#' @param vec A vector. The default is `NULL`. If `vec` is not `NULL`, the 
#' function will check `vec` instead of `ped[[col_spec]]`. 
#' @param force_NA A logical value indicating whether non-permissible should 
#' be replaced with `NA`s. The default is `TRUE`. 
#' @param force_remove A logical value indicating whether non-permissible 
#' values should be removed. The default is `FALSE`. 
#' @return The vector `ped[[col_spec]]` (or `vec` if `vec` is not `NULL`) with 
#' non-permissible values either removed and replaced. 
forcingNA_ifnot_contains <- function(ped, col_spec, allowed,
                                     default_NA_rep = NULL,
                                     NA_Ok = TRUE, vec = NULL,
                                     force_NA = TRUE, force_remove = FALSE) {
  if (is.null(vec)) {
    vec <- ped[[col_spec]]
  }
  contains_check <- .ifContains(vec, allowed, NA_Ok)
  if (!is.null(contains_check)) {
    if (force_NA & !force_remove) {
      if (sum(vec %in% contains_check) > 0) {
        msg <- sprintf(
          "Column %s contains %s that is not recognized. Only %s are valid inputs. Forcing default/NA.",
          col_spec,
          paste0(unique(contains_check), collapse = ","),
          paste0(allowed, collapse = ",")
        )
        rlang::inform(msg, level = "InvalidInputForcingNA")
        vec[vec %in% contains_check] <- NA
      }
    }
    if (force_remove & !force_NA) {
      if (sum(!vec %in% contains_check) > 0) {
        msg <- sprintf(
          "Column %s contains %s that is not recognized. Only %s are valid inputs. Forcing removal.",
          col_spec,
          paste0(unique(contains_check), collapse = ","),
          paste0(allowed, collapse = ",")
        )
        rlang::inform(msg, level = "InvalidInputForcingNA")
        vec <- vec[!vec %in% contains_check]
      }
    }
  }
  if (!is.null(default_NA_rep)) {
    if (sum(is.na(vec)) > 0) {
      msg <- sprintf(
        "Column %s does not support NA. Forcing to %s.",
        col_spec,
        default_NA_rep
      )
      rlang::inform(msg, level = "InvalidNA")
      vec[is.na(vec)] <- default_NA_rep
    }
  }
  return(vec)
}


#' Base function used to map values using a dictionary
#' 
#' @param dict A list with two named components that map to each other. The 
#' default is `PPP:::CANCER_NAME_MAP`, which maps short and long cancer 
#' names. 
#' @return A function where the input should be a keyword argument 
#' corresponding to one of the named components of `dict`. This will be treated 
#' as the "from" component. The return value will be the input mapped to the 
#' values of the other named component of `dict` (the "to" component).
#' @seealso \code{\link{.mapCancerNames}}, \code{\link{.mapGenderNames}} 
.mapDict <- function(dict = CANCER_NAME_MAP) {
  function(...) {
    args <- list(...)
    if (any(!names(args) %in% names(dict))) {
      rlang::abort("Some arguments cannot be found in dictionary",
        level = "MatchError"
      )
    }
    if (length(args) > 1) {
      rlang::abort("Argument list has a length larger than one",
        level = "MatchError"
      )
    }
    input_key <- names(args)[1]
    key_to_match <- names(dict)[!names(dict) %in% input_key]
    matched_indices <- match(args[[1]], dict[[input_key]])

    if (any(is.na(matched_indices))) {
      rlang::inform("Returned indices contains NA that will be removed, please double check!",
        level = "MatchWarning"
      )
      matched_indices = as.vector(na.omit(matched_indices))
    }
    return(dict[[key_to_match]][matched_indices])
  }
}


#' Convert cancer short names to long names or vice versa
#' 
#' @param ... A keyword argument that represents the "from" component of the 
#' dictionary, one of: 
#' * `short`: A character vector of short cancer names that should be mapped to 
#' long cancer names. 
#' * `long`: A character vector of long cancer names that should be mapped to  
#' short cancer names. 
#' @return A vector of the input values mapped to the "to" component of the 
#' dictionary (either `long` or `short`). 
#' @seealso \code{\link{.mapDict}}
.mapCancerNames <- .mapDict(CANCER_NAME_MAP)


#' Convert between numeric-valued sex and character-valued sex
#' 
#' @param ... A keyword argument that represents the "from" component of the 
#' dictionary, one of: 
#' * `number`: A numeric vector with values that are a subset of 
#' `c(0, 1, NA, -999)` and will be mapped to 
#' `c("Female", "Male", "All_Sexes", "All_Sexes")`. 
#' * `character`: A character vector with values that are a subset of 
#' `c("Female", "Male", "All_Sexes", "All_Sexes")` and will be mapped to 
#' `c(0, 1, NA, -999)`. 
#' @return A vector of the input values mapped to the "to" component of the 
#' dictionary (either `character` or `number`). 
#' @seealso \code{\link{.mapDict}}
.mapGenderNames <- .mapDict(list(
  number = c(0, 1, NA, -999),
  character = c("Female", "Male", rep("All_Sexes", 2))
))


#' Identifies values that are not allowed in a vector
#' 
#' @param vec A vector. 
#' @param allowed_vals A vector of the values that are allowed in `vec`. The 
#' default is `c(0, 1)`. 
#' @param NA_ok A logical value indicating whether `NA`s can be present in 
#' `vec`. The default is `TRUE`. 
#' @return Vector of values in `vec` that are not in `allowed_vals` (or `NA`, 
#' if `NA_ok = TRUE`). If no values were found, returns `NULL`. 
.ifContains <- function(vec, allowed_vals = c(0, 1),
                        NA_ok = TRUE) {
  if (NA_ok) {
    vec <- vec[!is.na(vec)]
  }
  ifin <- vec %in% allowed_vals
  if (!all(ifin)) {
    who_not_in <- vec[!ifin]
    if (length(who_not_in) != 0) {
      return(who_not_in)
    }
  }
  else {
    return(NULL)
  }
}


#' Evaluate function while ignoring `NA` and another predefined value
#' 
#' @param FUN The function to be applied. 
#' @param vec Vector over which to apply FUN. 
#' @param na_code Value to ignore in addition to `NA`. The default is `-999`. 
#' @return The output of `FUN(vec)`, while ignoring values in `c(NA, na_code)`. 
.safeGet <- function(FUN, vec, na_code = -999) {
  if (class(vec) == "list") {
    tryCatch(vec <- unlist(vec),
      error = function(c) "Input is not a unlistable list!"
    )
  }
  # If all components are NA
  if (all(is.na(vec))) {
    return(NA)
  }

  # If input is empty
  if (is.null(vec)) {
    return(NA)
  }

  # If there is NA
  vec <- vec[!is.na(vec)]

  # If there is matching NA code
  if (!is.null(na_code)) vec <- vec[vec != na_code]

  # If then there is nothing left
  if (length(vec) == 0) {
    return(NA)
  } else {
    return(FUN(vec))
  }
}


#' Force `NA` to be `FALSE`
#' 
#' @param vec A logical vector. 
#' @return `vec` with values of `NA` set to `FALSE`. 
.forceFalse <- function(vec) {
  vec[is.na(vec)] <- FALSE
  return(vec)
}


#' Get the random state
#' 
#' Uses `get0` to return `NULL` when `.Random.seed` doesn't exist, with 
#' `inherits = FALSE` to only search the global environment and not its parent. 
#' @return `.Random.seed` object. 
get_rand_state <- function() {
  get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}


#' Set a random state
#' 
#' @details 
#' Assigning `state = NULL` might lead to unwanted consequences. 
#' @param state Current state. 
#' @return None. 
set_rand_state <- function(state) {
  if (!is.null(state)) {
    assign(".Random.seed", state, envir = .GlobalEnv, inherits = FALSE)
  }
}

#' converts age of 1st child birth to one of four categories
#' 
#' @param kids.ages numeric vector of subject's children's ages; can be NULL
#' or length 0
#' @param curage the subject's current age as a numeric value
#' @return one of c("<30/nulliparous","30-39","40+","Unk")
.categorize.1st.birth.age <- function(kids.ages, curage){
  
  # has children
  if(length(kids.ages) > 0){
    
    # NA age present
    if(sum(is.na(kids.ages)) > 0 | sum(kids.ages == -999) > 0){
      fb.cat <- "Unk"
      return(fb.cat)
      
      # no NA ages present
    } else {
      first.birth.age <- curage - max(kids.ages)
      if(first.birth.age < 12){
        fb.cat <- "Unk"
        return(fb.cat)
      }
    }
    
    # no children
  } else {
    fb.cat <- "<30/nulliparous"
    return(fb.cat)
  }
  
  # check first birth age compared to diag age
  if(first.birth.age < 30){
    fb.cat <- "<30/nulliparous"
  } else if(first.birth.age < 40){
    fb.cat <- "30-39"
  } else {
    fb.cat <- "40+"
  }
  return(fb.cat)
}


#' Get CBC model parameters for a subject. 
#' 
#' These parameters can be used for subsetting a contralateral breast cancer 
#' penetrance array for carriers and/or creating noncarrier penetrances using 
#' \code{\link{CBCnoncarrierCrude}}.
#'
#' @param ped a checked pedigree data frame object returned by 
#' \code{\link{checkFam}} with at least columns `ID`, `MotherID`, `FatherID`, 
#' `Sex`, `CurAge`, `isAffBC`, `AgeBC`. Optional columns are `race`, 
#' `FirstBCType`, `ER`, `AntiEstrogen`, `HRPreneoplasia`, and `BreastDensity`. 
#' @param id numeric value, a subject ID in the pedigree for whom to find the 
#' results.
#' @param cbc.penets CBC carrier penetrance array structured like 
#' `PanelPRODatabase$Contralateral`.
#' @param imputed.ages list of imputed ages as returned by 
#' \code{\link{checkFam}}, must contain a `AgeBC` component.
#' @return named list containing:
#' * `Gene`: string, subject's genotype; either `"noncarrier"` or value from 
#' `CBC_GENE_VARIANT_TYPES`.
#' * `Sex`: string, sex of the subject, either `c("Female", "Male")`.
#' * `HadBC`: binary indicator if subject has had a primary breast cancer yet.
#' * `FirstBCAge`: numeric value, age of subject's first breast cancer.
#' * `FirstBCType`: numeric value indicating whether the subject's first breast 
#' cancer was pure invasive or mixed invasive and DCIS, pure invasive or 
#' unknown type = `3` and mixed invasive and DCIS = `2`. PPP does not handle 
#' pure DCIS breast cancers.
#' * `FDRBC`: numeric value indicating if the subject has first degree relatives 
#' with a history of breast cancer, Yes = `1` or No = `2`.
#' * `ER`: numeric value indicating subject's ER result, Unk = `3`, Pos = `2`, 
#' Neg = `1`.
#' * `FirstBirth`: numeric value indicating age category when subject's first 
#' child was born, <30/nulliparous = `1`, 30-39 = `2`, 40+ = `3`, Unk = `4`.
#' * `AntiEstrogen`: numeric value indicating if anti-estrogen therapy was used 
#' to treat the 1st breast cancer, Yes = `1`, No = `2`, Unk = `3`.
#' * `HRPreneoplasia`: numeric value indicating a history of high risk 
#' preneoplasia, Yes = `1`, No/Unk = `2`.
#' * `BreastDensity`: numeric value indicating BI-RADS breast density category. 
#' These codes are different from non-Hispanic Black women and everyone else as 
#' required by the `cbcrisk` package. For everyone else: `"a"` or entirely 
#' fat = `4`, `"b"` or scattered = `3`, `"c"` or heterogeneously dense = `2`, 
#' `"d"` / extremely dense = `1`, Unk = `5`. For women who ARE non-Hispanic 
#' Black: `"a"` (entirely fat) or `"b"` (scattered) are `1`, `"c"` (heterogeneously 
#' dense) or `"d"` (extremely dense) are `2`, Unk is `3`.
#' @details 
#' * If the subject has more than one positive CBC-linked gene, the gene 
#' selected will be the one with the highest CBC crude penetrance at age 85. If 
#' no germline testing information is available for the CBC genes, "noncarrier" 
#' status is assumed. 
#' * If age of first breast cancer has absent from the raw pedigree, a randomly 
#' selected imputed age from \code{\link{checkFam}} will be used instead.
#' * The outputs of this function are designed to be compatible 
#' with the inputs to \code{\link{CBCnoncarrierCrude}}.
getCBCParams <- function(ped, id, cbc.penets, imputed.ages){
  
  # Subject's sex
  tmp.sex <- ped$Sex[which(ped$ID == id)]
  
  # race - note that cbcrisk profile values are different for breast density, 
  # FDR Fam Hx with BC, and age at 1st BC diagnosis
  if("race" %in% colnames(ped)){
    race <- ped$race[which(ped$ID==id)]
    # Unk/other = 1, Non-hispanic Black = 2
    tmp.race <- ifelse(race == "Black", 2, 1)
    
  } else {
    tmp.race <- 1
  }
  
  # Did a 1st Breast Cancer occur?
  bc1.status <- ped$isAffBC[which(ped$ID == id)]
  if(!is.na(bc1.status) & bc1.status == 1){ # primary BC occurred, so CBC is calculable
    
    # First BC age (use random previously imputed value if missing)
    bc1.age <- ped$AgeBC[which(ped$ID == id)]
    if(bc1.age == -999){
      bc1.age <- sample(imputed.ages$AgeBC[,as.character(id)], size = 1)
    }
    
    # First BC type
    if("FirstBCType" %in% colnames(ped)){
      bc1.type <- ped$FirstBCType[which(ped$ID == id)]
      # Invasive or Unk = 3, Invasive_DCIS = 2, DCIS = 1 but DCIS not used
      bc1.type <- ifelse(bc1.type == -999 | bc1.type == "Invasive", 3, 
                         ifelse(bc1.type == "Invasive_DCIS", 2, "invalid"))
    } else {
      bc1.type <- 3
    }
    
    # ER result
    if("ER" %in% colnames(ped)){
      er.stat <- ped$ER[which(ped$ID==id)]
      # Unk = 3, Pos = 2, Neg = 1
      tmp.er <- ifelse(er.stat == -999, 3, 
                       ifelse(er.stat == 1, 2, 
                              ifelse(er.stat == 0, 1, "invalid")))
    } else {
      tmp.er <- 3
    }
    
    # Anti-Estrogen Therapy
    if("AntiEstrogen" %in% colnames(ped)){
      antiEst <- ped$AntiEstrogen[ped$ID == id]
      # Yes = 1, No = 2, Unk = 3
      antiEst <- ifelse(antiEst == -999, 3, 
                        ifelse(antiEst == 0, 2, 
                               ifelse(antiEst == 1, 1, "invalid")))
    } else {
      antiEst <- 3
    }
    
    # High Risk Preneoplasia
    if("HRPreneoplasia" %in% colnames(ped)){
      HRpre <- ped$HRPreneoplasia[ped$ID == id]
      # Yes = 1, No/Unk = 2
      HRpre <- ifelse(HRpre == -999 | HRpre == 0, 2, 
                      ifelse(HRpre == 1, 1, "invalid"))
    } else {
      HRpre <- 2
    }
    
    # 1st tumor size (only used for non-Hispanic Black women)
    if("FirstBCTumorSize" %in% colnames(ped)){
      bc1.tsize <- ped$FirstBCTumorSize[ped$ID == id]
      bc1.tsize <- dplyr::recode(bc1.tsize, "-999"=4, "T0/T1/T2"=1, "T3/T4"=2, "Tis"=3)
    } else {
      bc1.tsize <- 4
    }
    
  } else { # no CBC because primary BC has not occurred, skip calculations
    bc1.age <- bc1.type <- tmp.er <- antiEst <- HRpre <- bc1.tsize <- NA
    
  } 
  
  # 1st degree relatives (FDR) with breast cancer (Yes, No, Unk)
  fdr.ped <- ped[,c("ID","MotherID","FatherID")]
  fdr.ped[is.na(fdr.ped)] <- -999
  fdrs <- .firstDegreeRelative(fdr.ped, which(ped$ID==id))
  fdrs.ind <- fdrs$index
  fdrs.with.BC <- sum(ped$isAffBC[fdrs.ind], na.rm=T)
  # Yes = 1, No = 2, Unk = 3 (Unk not used)
  fdr.bc <- ifelse(fdrs.with.BC > 0, 1, 2)
  
  # Age at First Birth
  kids.ages <- ped$CurAge[fdrs$index[which(fdrs$lineage == "C")]]
  curage <- ped$CurAge[which(ped$ID == id)]
  tmp.1birth <- .categorize.1st.birth.age(kids.ages, curage)
  fb.cutoff <- ifelse(tmp.1birth == "30-39", 30, 
                      ifelse(tmp.1birth == "40+", 40, 88))
  # <30/nulliparous = 1, 30-39 = 2, 40+ = 3, Unk = 4
  tmp.1birth <- ifelse(tmp.1birth == "<30/nulliparous" | bc1.age <= fb.cutoff, 1,
                       ifelse(tmp.1birth == "30-39", 2,
                              ifelse(tmp.1birth == "40+", 3, 
                                     ifelse(tmp.1birth == "Unk", 4, "invalid"))))
  
  # BI-RADS Breast Density Descriptor
  if("BreastDensity" %in% colnames(ped)){
    
    # for women who are not non-Hispanic Black
    if(tmp.race == 1){
      bDensity <- ped$BreastDensity[ped$ID == id]
      # a / entirely fat = 4, b / scattered = 3, c / hetero. dense = 2, 
      # d / extremely dense = 1, Unk = 5
      bDensity <- dplyr::recode(bDensity, "-999"=5, "a"=4, "b"=3, "c"=2, "d"=1)
      bDensity <- ifelse(!bDensity %in% 1:5, "invalid", bDensity)
      
      # for non-Hispanic Black women
    } else if(tmp.race == 2){
      bDensity <- dplyr::recode(bDensity, "-999"=3, "a"=1, "b"=1, "c"=2, "d"=2)
      bDensity <- ifelse(!bDensity %in% 1:3, "invalid", bDensity)
    }
    
    # no breast density information in pedigree
  } else {
    if(tmp.race == 1){ # for women who are not non-Hispanic Black
      bDensity <- 5
    } else if(tmp.race == 2){ # for non-Hispanic Black women
      bDensity <- 3
    }
  }
  
  # CBC gene
  gene.cols <- intersect(colnames(ped), unlist(CBC_GENE_VARIANT_TYPES))
  if(length(gene.cols) > 0){
    gene.results <- ped[which(ped$ID == id), gene.cols]
    gene.results[gene.results == -999] <- 0 # noncarrier if no germline info
    if(sum(gene.results, na.rm = T) > 1){
      # if multiple pos genes, select the one with the highest penetrance at 85
      tmp.genes <- colnames(gene.results)[which(gene.results == 1)]
      tmp.gene.85.penets <- cbc.penets[tmp.genes, tmp.sex, bc1.age, "Crude", 85]
      tmp.gene <- names(which(tmp.gene.85.penets == max(tmp.gene.85.penets, na.rm=T)))
      
    } else if(sum(gene.results, na.rm = T) == 1){
      tmp.gene <- colnames(gene.results)[which(gene.results == 1)]
    } else if(sum(gene.results, na.rm = T) == 0){
      tmp.gene <- "noncarrier" 
    }
  } else { # assume noncarrier if no germline information
    tmp.gene <- "noncarrier"
  }
  
  # check for invalid parameters
  return.list <- list(Gene = tmp.gene, Sex = tmp.sex, HadBC = bc1.status, 
                      FirstBCAge = bc1.age, FirstBCType = bc1.type, FDRBC = fdr.bc, 
                      ER = tmp.er, FirstBirth = tmp.1birth, AntiEstrogen = antiEst,
                      HRPreneoplasia = HRpre, BreastDensity = bDensity, 
                      BC1TSize = bc1.tsize,
                      race = tmp.race)
  
  if(any(unlist(return.list)[!is.na(unlist(return.list))] == "invalid")){
    rlang::abort(sprintf(
      "At least one CBC parameter in getCBCParams() is 'invalid'."), 
      level = "invalidCBCParam")
  }
  return(return.list)
}


#' Contralateral Breast Cancer Penetrances
#' 
#' Get a subject's age specific crude penetrances of contralateral breast cancer, 
#' based on eight subject characteristics.
#' 
#' @param diag.age numeric value, age of 1st breast cancer diagnosis
#' @param bc1.type numeric value, type of first breast cancer, 
#' `"Invasive"` or Unk = `3`, `"Invasive_DCIS"` = `2`, `"DCIS"` = `1` 
#' (`"DCIS"` not used in PPP)
#' @param fdr.status numeric value, if first degree relatives have history of 
#' breast cancer, Yes = `1`, No = `2`, Unk = `3` (Unk not used in PPP)
#' @param ER.status numeric value, ER test result, Unk = `3`, Pos = `2`, Neg = `1`
#' @param birth1.age numeric value, age at first child birth, 
#' <30/nulliparous = `1`, 30-39 = `2`, 40+ = `3`, Unk = `4`
#' @param antiEst numeric value, indicates if anti-estrogen therapy was used to 
#' treat the first breast cancer, Yes = `1`, No = `2`, Unk = `3`
#' @param HighRiskP numeric value, indicates history of high risk preneoplasia, 
#' Yes = `1`, No/Unk = `2`
#' @param BD numeric value based on BI-RADS breast density ratings `"a"` to `"d"`. 
#' The categories differ for non-Hispanic Black women versus other women. For 
#' women who are NOT non-Hispanic Black: `"a"` or entirely fat = `4`, `"b"` or 
#' scattered = `3`, `"c"` or heterogeneously dense = `2`, `"d"` or extremely 
#' dense = `1`, Unk = `5`. For women who ARE non-Hispanic Black: `"a"` (entirely 
#' fat) or `"b"` (scattered) are `1`, `"c"` (heterogeneously dense) or `"d"` 
#' (extremely dense) are `2`, Unk is `3`.
#' @param bc1.t.size number representing a category of 1st BC tumor size, 
#' `1` for sizes T0/T1/T2, `2` for T3/T4, `3` for Tis, or `4` for unknown size.
#' @param sex string, subject's sex, one of `c("Male", "Female")`, default is 
#' `"Female"`. If `sex` is not `"Female"`, function will return an all `0`'s vector 
#' because the `cbcrisk` model is not valid for males.
#' @param race number, `2` if non-Hispanic Black, `1` otherwise.
#' @param had1stBC numeric value, indicates if subject previously diagnosed with 
#' primary breast cancer, Yes = `1`, No = `0`. If `had1stBC` is not `1` then then 
#' function will return an all `0`'s vector because the model relies on a known 
#' diagnosis age of primary breast cancer.
#' @details
#' Uses the cbcrisk function from the \href{https://personal.utdallas.edu/~sxb125731/}{cbcrisk package} 
#' for female CBC penetrances and subsets the an array of CBC noncarrier 
#' penetrances from ASK2ME for males. The inputs of this function are designed 
#' to be compatible with the outputs of the \code{\link{getCBCParams}} function.
#' @return A named vector of penetrances from age `1` to age `PPP:::MAXAGE` 
#' where the vector's names are the ages. Penetrances are `0` from age `1` 
#' through the age of first breast cancer diagnosis. Penetrances are all `0` if 
#' `had1stBC` is not `1`.
CBCnoncarrierCrude <- function(diag.age, bc1.type, fdr.status, ER.status, 
                               birth1.age, antiEst, HighRiskP, BD, bc1.t.size,
                               sex, race, had1stBC = 1){
  
  # return all 0 penetrance vector if cbcrisk not applicable
  if(had1stBC != 1){
    penets <- rep(0, PPP:::MAXAGE)
    names(penets) <- 1:(PPP:::MAXAGE)
    return(penets)
  }
  
  # subset male noncarrier array for males
  if(sex == "Male"){
    penets <- MaleCBCNoncarrierArray[as.character(diag.age), "Crude",]
    names(penets) <- 1:(PPP:::MAXAGE)
    return(penets)
  }
  
  ### FEMALES WITH PREVIOUS 1st BC
  # convert age to cbcrisk profile compatible category based on race
  if(race == 1){ # for women who are not non-Hispanic Black
    diag.age.cat <- ifelse(diag.age < 30, 1, ifelse(diag.age < 40, 2, 3))
  } else if(race == 2){ # for non-Hispanic Black women
    diag.age.cat <- ifelse(diag.age < 40, 1, 2)
  }
  
  # adjust CBC diagnosis age if not between 18 and 88 to prevent cbcrisk error
  if(diag.age > 88){
    tmp.diag.age <- 88
  } else if(diag.age < 18){
    tmp.diag.age <- 18
  } else {
    tmp.diag.age <- diag.age
  }
  
  ## cbcrisk profile
  # use different profile for Black women (CBCRisk-Black, Sajal et al., 2022)
  # than everyone else
  if(race == 1){ # everyone else
    profile <- c(diag.age.cat,
                 antiEst,
                 fdr.status,
                 HighRiskP,
                 BD,
                 ER.status,
                 bc1.type,
                 birth1.age)
  } else if(race == 2){ # non-Hispanic Black women
    profile <- c(BD,
                 fdr.status, 
                 bc1.t.size,
                 diag.age.cat)
  }
  
  # get CBC cancer risks
  cum.risks.df <- cbcrisk::cbcrisk(race = race,
                                   profile = profile, 
                                   start.age = tmp.diag.age, 
                                   pred.year = 1,
                                   print.output = F)$risk
  
  # convert cumulative risks to penetrances
  cum.risks.df$penet <- rep(0, nrow(cum.risks.df))
  for(row in nrow(cum.risks.df):1){
    if(row > 1){
      previous.risk <- cum.risks.df$'CBC risk(%)'[(row-1)] / 100
    } else {
      previous.risk <- 0
    }
    cum.risks.df$penet[row] <- cum.risks.df$'CBC risk(%)'[row] / 100 - previous.risk
  }
  penets <- cum.risks.df$penet
  
  # extrapolate penetrances if CBC diagnosis was age 89 to 94
  if(diag.age > 88){
    penet.89 <- penets[1]
    interp.old <- approx(x = c(89, PPP:::MAXAGE+1), y = c(penet.89, 0), 
                         n = (PPP:::MAXAGE+1)-89+1)$y
    names(interp.old) <- 89:(PPP:::MAXAGE+1)
    if(diag.age < PPP:::MAXAGE){
      interp.old <- interp.old[which(names(interp.old) %in% 
                                       seq(diag.age+1, PPP:::MAXAGE))]
      penets <- interp.old
    } else if(diag.age == PPP:::MAXAGE) {
      penets <- c()
    }
    
    # adjust penetrances if CBC diagnosis was less than age 89 (cbcrisk limitation)
  } else {
    
    # extrapolate penetrances if CBC diagnosis age was 0 to 17 (cbcrisk limitation)
    if(diag.age < 18){
      penet.19 <- penets[1]
      interp.yng <- approx(x = c(0, 19), y = c(0, penet.19), n = 19+1)$y
      names(interp.yng) <- 0:19
      interp.yng <- interp.yng[which(names(interp.yng) %in% seq(diag.age+1, 18))]
      penets <- c(interp.yng, penets)
    } 
    
    # extrapolate values from 89 to 95, exclusively (cbcrisk limitation)
    penet.89 <- penets[length(penets)]
    interp.old <- approx(x = c(89, PPP:::MAXAGE+1), y = c(penet.89, 0), 
                         n = (PPP:::MAXAGE+1)-89+1)$y
    interp.old <- interp.old[2:(length(interp.old)-1)]
    penets <- c(penets, interp.old)
  }
  
  # final vector of penetrances, pad w/ 0's from 1 to 1st BC
  penets <- setNames(c(rep(0, diag.age), penets), 1:(PPP:::MAXAGE))
  return(penets)
}
