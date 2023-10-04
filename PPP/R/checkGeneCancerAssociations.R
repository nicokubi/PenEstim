#' Check for gene-cancer associations in database
#' 
#' Internal function that checks for all gene-cancer associations in a 
#' database, based on the net penetrances at age 50. The user can optionally 
#' restrict the sexes and races to consider. The function will either plot 
#' the resulting matrix of associations or return it. 
#' 
#' @param db A model-specific database returned by \code{\link{buildDatabase}}. 
#' The default is `NULL`, in which case a database of all currently-supported 
#' genes and cancers in the package will be used. 
#' @param isPlot A boolean flag indicating whether or not the function should 
#' plot the matrix of gene-cancer associations. The default is `TRUE`. If 
#' `FALSE`, the plot will be suppressed and the matrix will be returned. 
#' @param Race A character vector of races to consider when checking for 
#' associations. The default is `NULL`, in which case all races in `db` will be 
#' considered when checking for gene-cancer associations. 
#' @param Sex A character vector of sexes to consider when checking for 
#' associations. The default is `NULL`, in which case all sexes in `db` will be 
#' considered when checking for gene-cancer associations. 
#' @param shading_color The color to use for shading in gene-cancer 
#' associations when `isPlot` is `TRUE`. The default is `"dimgray"`. 
#' 
#' @return 
#' If `isPlot` is `TRUE`, the function will return a plot. If `isPlot` is 
#' `FALSE`, a boolean matrix with rows corresponding to the cancers in `db` and 
#' the columns corresponding to the genes in `db`. Associations that exist in 
#' `db` are `TRUE` in the matrix. 
#' 
#' @examples
#' # Plot all gene-cancer associations in PanelPRODatabase
#' # Not run
#' # PPP:::checkGeneCancerAssociations()
#' @md
checkGeneCancerAssociations = function(db = NULL, isPlot = TRUE, 
                                       Race = NULL, Sex = NULL, 
                                       shading_color = "dimgray") {
  
  # If no database is specified, default to all supported genes and cancers
  if (is.null(db)) {
    db = buildDatabase(genes = PPP:::GENE_TYPES, 
                       cancers = PPP:::CANCER_TYPES[-which(PPP:::CANCER_TYPES == "Contralateral")])
  }
  
  # If Race is unspecified, use all races in database
  if (is.null(Race)) {
    Race = dimnames(db$penet$penet_c)$Race
  }
  # If Sex is unspecified, use all sexes in database
  if (is.null(Sex)) {
    Sex = dimnames(db$penet$penet_c)$Sex
  }
  
  # get all ages
  Age = dimnames(db$penet$penet_c)$Age
  
  # Initialize association matrix
  assoc_mat = matrix(0, nrow = dim(db$penet$penet_c)[1], 
                     ncol = dim(db$penet$penet_c)[2], 
                     dimnames = list(Cancer = dimnames(db$penet$penet_c)$Cancer, 
                                     Gene = dimnames(db$penet$penet_c)$Gene))
  
  # Identify which column corresponds to SEER
  SEER_idx = which(colnames(assoc_mat) == "SEER")
  
  # Remove SEER from the association matrix
  assoc_mat = assoc_mat[,-SEER_idx, drop = F]
  
  # Iterate through all races and sexes to check for associations
  for (race in Race) {
    for (sex in Sex) {
      # Extract cumulative SEER penetrances 
      SEER_pen = rowSums(db$penet$penet_c[, SEER_idx, race, sex, Age, "Net", drop = F], dims = 1)
      # Compare all other cumulative penetrances to the SEER penetrances
      # cumulative gene penetrances
      gene_pen = rowSums(db$penet$penet_c[, -SEER_idx, race, sex, Age, "Net", drop = F], dims = 2)
      assoc_mat = assoc_mat + 
        apply(gene_pen, 2, 
              function(x) { x != SEER_pen })
    }
  }
  
  # Convert matrix of counts to boolean
  assoc_mat = assoc_mat > 0
  
  # Sort rows and columns alphabetically
  assoc_mat = assoc_mat[order(rownames(assoc_mat)),order(colnames(assoc_mat)), 
                        drop = F]
  
  # Generate a plot, if requested
  if (isPlot == TRUE) {
    # prep data for compatibility with image() command
    p.dat <- t(assoc_mat)
    p.dat <- p.dat[, order(colnames(p.dat), decreasing = T), drop = F]
    
    ## address edge cases
    # case 1: for all TRUE matrices, ensure color is grey, not white
    if(all(p.dat)){
      img.cols <- c(shading_color)
    } else {
      img.cols <- c("white", shading_color)
    }

    # case 2: for 1 row (1 gene) matrices
    if(nrow(p.dat) == 1){
      r.len <- 0
    } else {
      r.len <- seq(1, 0, length = nrow(p.dat))
    }
    
    # Image plot of association matrix
    # image() transposes matrix; nrow(), ncol(), etc are not intuitive
    par(mar = c(1, 9, 8, 1))
    image(p.dat, col = img.cols, xaxt = "n", yaxt = "n")
    
    # grid lines, axes and labels
    box(col = "black") # border lines
    grid(nx = nrow(p.dat), ny = ncol(p.dat), # internal lines
         col = "black", lty = "solid")
    axis(3, at = r.len, # top axis of gene labels
         labels = sort(formatGeneNames(rownames(p.dat),
                                       format = "drop_hetero_anyPV"), 
                       decreasing = T), # decreasing sorting b/c at = seq(1,0)
         las = 2, tick = FALSE, line = -0.7)
    axis(2, at = seq(0, 1, length = ncol(p.dat)), # left axis of cancer labels
         labels = colnames(p.dat),
         las = 2, tick = FALSE, line = -0.7)
    
  } else {
    
    # returns logical matrix
    return(assoc_mat)
  }
}
