#' The transmission matrix for a single genetic locus
#'
#' A function to calculate the transmission matrix for a single autosomal
#' genetic locus with an arbitrary number of alleles and unphased genotypes,
#' based on Mendel's laws of inheritance.
#'
#' @param n_alleles A positive integer, interpreted as the number of possible
#' alleles at the genetic locus.
#'
#' @param annotate A logical flag. When `FALSE` (the default), the function
#' returns a matrix suitable to be used as the `trans` argument of
#' \code{\link{pedigree_loglikelihood}}. When `TRUE`, the function annotates
#' this matrix (and converts it to a data frame) to make the output more
#' easily understood by humans.
#' 
#' @param nonviable A logical flag. When `FALSE` (the default), the function
#' runs as per usual. When 'TRUE' the function assumes that the homozygous carrier
#' is not viable. Hence only two genotypes are modeled. 
#'
#' @details When `annotate` is `FALSE`, this function returns a matrix of
#' genetic transmission probabilities, whose rows corresponding to the possible
#' joint parental genotypes and whose columns corresponding to the possible
#' offspring genotypes.  There are `ngeno = n_alleles * (n_alleles + 1) / 2` possible
#' unphased genotypes, and by choosing an order on these genotypes (which can be
#' viewed by setting `annotate` to `TRUE`, see below)
#' we can label the set of possible genotypes as `1:ngeno`.
#' Then the `(ngeno * gm + gf - ngeno, go)`th element of the outputted matrix is
#' the conditional probability that a person has genotype `go`, given that his
#' or her biological mother and father have genotypes `gm` and `gf`,
#' respectively.
#'
#' When `annotate` is `TRUE`, the function converts this matrix to a data frame,
#' adds column names giving the offspring genotype corresponding to each
#' column, and adds columns `gm` and `gf` describing the parental genotypes
#' corresponding to each row.  In this data frame, genotypes are written
#' in the usual form `1/1, 1/2, ...` for the alleles `1:n_alleles`.
#'
#' Note that if the output of this function is to be used as the `trans`
#' argument of \code{\link{pedigree_loglikelihood}} then the `annotate` option
#' must be set to `FALSE`.
#'
#' @return Either a matrix of genetic transmission probabilities suitable to be
#' used as the `trans` argument of \code{\link{pedigree_loglikelihood}}
#' (if `annotate` is `FALSE`), or a data frame that is an annotated version of
#' this matrix (if `annotate` is `TRUE`).
#'
#' @export
#'
#' @examples
#' # The transition matrix for a biallelic, autosomal locus with unphased genotypes
#' trans_monogenic(2)
#' trans_monogenic(2, annotate = TRUE)
#'

trans_monogenic2 <- function(n_alleles, annotate = FALSE, nonviable) {
    # List of possible genotypes
    eg <- cbind(rep(1:n_alleles, each = n_alleles), rep(1:n_alleles, times = n_alleles))
    f <- function(x) {
        if (x[1] > x[2]) {
            return("")
        } else {
            return(paste(x, collapse = "/"))
        }
    }
    geno.list <- apply(eg, 1, f)
    geno.list <- geno.list[geno.list != ""]
    ng <- length(geno.list)

    # Exclude non-viable genotype if nonviable = TRUE
    if (nonviable && n_alleles == 2) {
        non_viable_geno <- paste(n_alleles, n_alleles, sep = "/")
        geno.list <- geno.list[!geno.list %in% non_viable_geno]
        ng <- length(geno.list)
    }

    # Initialize the transmission matrix
    trans <- matrix(0, ng^2, ng)

    # Populate the transmission matrix
    for (gm in 1:ng) {
        for (gf in 1:ng) {
            am.list <- strsplit(geno.list[gm], "/", fixed = TRUE)[[1]]
            af.list <- strsplit(geno.list[gf], "/", fixed = TRUE)[[1]]
            for (i in 1:2) {
                for (j in 1:2) {
                    offspring_geno <- sort(c(am.list[i], af.list[j]), decreasing = FALSE)
                    go <- which(geno.list == paste(offspring_geno, collapse = "/"))
                    if (length(go) > 0) {
                        trans[ng * gm + gf - ng, go] <- trans[ng * gm + gf - ng, go] + 0.25
                    }
                }
            }
        }
    }

    # Normalize each row so that probabilities sum up to 1
    for (i in 1:nrow(trans)) {
        row_sum <- sum(trans[i, ])
        if (row_sum > 0) {
            trans[i, ] <- trans[i, ] / row_sum
        }
    }

    # Annotate the matrix if requested
    if (annotate) {
        trans <- data.frame(trans)
        colnames(trans) <- geno.list
        trans$gm <- rep(geno.list, each = ng)
        trans$gf <- rep(geno.list, times = ng)
    }

    return(trans)
}
