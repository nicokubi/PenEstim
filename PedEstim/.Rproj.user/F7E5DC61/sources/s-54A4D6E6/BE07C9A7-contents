#' Format gene names
#'
#' Convenience function for formatting gene names into shorter, more visually 
#' pleasing names. A message will be printed if choice of formatting results in 
#' non-unique gene names. 
#' 
#' @param gene_names A character vector of gene names, assumed to be of the 
#' format GENE_HETERO/HOMO_VARIANT, where GENE is the gene, HETERO/HOMO is 
#' either `"hetero"` or `"homo"` for the respective mutations, and VARIANT is 
#' either `"anyPV"` or a specific variant.
#' @param format One of the following character strings: 
#' * `"drop_hetero_anyPV"` (default): Removes `"hetero"` and `"anyPV"` from 
#' gene names, but retains `"homo"` and specific variants. 
#' * `"drop_hetero"`: Removes "`hetero"` from gene names, but retains `"homo"` 
#' and variant information. 
#' * `"drop_anyPV"`: Removes `"anyPV"` from gene names, but retains `"homo"` 
#' and specific variants. 
#' * `"drop_PVs"`: Removes variant information from gene names, but retains 
#' `"hetero"` and `"homo"`. 
#' * `"drop_heterohomo"`: Removes `"hetero"` and `"homo"` from gene names, 
#' but retains variant information. 
#' * `"only_gene"`: Removes `"hetero"`, `"homo"`, and all variant information 
#' from gene names, retaining only the genes. 
#' * `"full":` Returns the gene names, unmodified. This is also the behavior if 
#' an unsupported `format` is specified. 
#' 
#' @return A character vector of formatted gene names. 
#' 
#' @examples 
#' gene_names = c("MUTYH_hetero_anyPV", "MUTYH_homo_anyPV", 
#'                "CHEK2_hetero_1100delC", 
#'                "MUTYH_hetero_anyPV.CHEK2_hetero_1100delC", 
#'                "MUTYH_homo_anyPV.CHEK2_hetero_1100delC")
#' PPP:::formatGeneNames(gene_names)
#' PPP:::formatGeneNames(gene_names, format = "only_gene")
formatGeneNames = function(gene_names, format="drop_hetero_anyPV") {
  if (format == "drop_hetero_anyPV") {
    gene_names = sapply(strsplit(gene_names, "\\."), function(x) {
      paste(sub("_anyPV", "", sub("_hetero", "", x)), collapse = ".")
    })
  } else if (format == "drop_hetero") {
    gene_names = sapply(strsplit(gene_names, "\\."), function(x) {
      paste(sub("_hetero", "", x), collapse = ".")
    })
  } else if (format == "drop_anyPV") {
    gene_names = sapply(strsplit(gene_names, "\\."), function(x) {
      paste(sub("_anyPV", "", x), collapse = ".")
    })
  } else if (format == "drop_PVs") {
    gene_names = sapply(strsplit(gene_names, "\\."), function(x) {
      paste(sub("_[^_]+$", "", x), collapse = ".")
    })
  } else if (format == "drop_heterohomo") {
    gene_names = sapply(strsplit(gene_names, "\\."), function(x) {
      paste(sub("_.*_", "_", x), collapse = ".")
    })
  } else if (format == "only_gene") {
    gene_names = sapply(strsplit(gene_names, "\\."), function(x) {
      paste(sub("(\\[.*?\\])", "", sub("_.*", "", x)), collapse = ".")
    })
  } else if (format != "full") {
    rlang::inform(paste(format, "is not a supported format, returning full names."))
  }
  
  if (any(duplicated(gene_names))) {
    rlang::inform("Choice of formatting results in non-unique gene names.")
  }
  
  return(gene_names)
}
