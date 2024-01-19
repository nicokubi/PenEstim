# Load Libraries
library(clipp)
library(PPP)
library(stats4)
library(dplyr)

dat <- prepAges(carrierProbandFamilies_cohPedigree_MSH6_COL_NL)

mhLogLikelihood_clipp(c(50,40,1,0), dat, af = 0.001,"Colorectal")

penet <- t(sapply(1:nrow(combined_df), function(i) penet.fn(i, combined_df, alpha, beta, gamma, delta, baselineRisk)))

mhLogLikelihood_clipp <- function(paras, families, af, cancer_type) {
    # Function to Calculate Lifetime Risk
    calculate_lifetime_risk <- function(cancer_type, gene = "SEER", race = "All_Races", sex = "Male", type = "Crude", data) {
        # Check if dimnames are available and correct
        if (is.null(data$Penetrance) || is.null(attr(data$Penetrance, "dimnames"))) {
            stop("Penetrance data or its dimension names are not properly defined.")
        }

        dim_names <- attr(data$Penetrance, "dimnames")
        required_dims <- c("Cancer", "Gene", "Race", "Sex", "Age", "PenetType")
        if (!all(required_dims %in% names(dim_names))) {
            stop("One or more required dimensions are missing in Penetrance data.")
        }

        # Function to safely extract index
        get_index <- function(dim_name, value) {
            idx <- which(dim_names[[dim_name]] == value)
            if (length(idx) == 0) {
                stop(paste("Value", value, "not found in dimension", dim_name))
            }
            idx
        }

        # Extracting indices for each dimension except Age
        cancer_index <- get_index("Cancer", cancer_type)
        gene_index <- get_index("Gene", gene)
        race_index <- get_index("Race", race)
        sex_index <- get_index("Sex", sex)
        type_index <- get_index("PenetType", type)

        # Subsetting Penetrance data for all ages using indices
        lifetime_risk <- data$Penetrance[cancer_index, gene_index, race_index, sex_index, , type_index]
        return(lifetime_risk)
    }

    # Penetrance Functions
    penet.fn <- function(i, data, alpha, beta, gamma, delta,baselineRisk) {
        max_age <- 94 # Assuming "inc" is defined elsewhere
        if (data$age[i] == 0) {
            penet.i <- c(1, 1) # Assuming people aged 0 years are all unaffected
        } else {
            # Find the index corresponding to the age of individual i
            age_index <- match(data$age[i], 1:94)
            # Check if the age is within the valid range of ages in your vector
            if (!is.na(age_index)) {
                # Extract the corresponding estimate
                nc.pen <- baselineRisk[age_index]
            } else {
                print("Age not found in the range of ages.")
            }
            # Weibull hazard and survival for carriers
            c.pen <- dweibull(data$age[i] - delta, alpha, beta) * gamma

            # Penetrance calculations based on genotype
            penet.i <- c(1 - nc.pen, 1 - c.pen) # if person is not affected
            if (data$aff[i] == 1) penet.i <- c(nc.pen, c.pen)
        }
        if (data$geno[i] == "1/1") penet.i[-1] <- 0
        if (data$geno[i] == "1/2") penet.i[-2] <- 0
        return(penet.i)
    }

    # Function to transform each data frame
    transform_df <- function(df) {
        df %>%
            select(
                individual = SubjectID,
                family = PedigreeID,
                mother = MotherID,
                father = FatherID,
                aff = isAffCOL,
                sex = Sex,
                age = CurAge,
                geno = MSH6
            ) %>%
            mutate(
                mother = ifelse(is.na(mother), NA, paste0("ora", sprintf("%03d", mother))),
                father = ifelse(is.na(father), NA, paste0("ora", sprintf("%03d", father))),
                geno = ifelse(is.na(geno), "", ifelse(geno == 1, "1/2", geno))
            )
    }
    
    combined_df <- do.call(rbind, lapply(families, transform_df))

    # Parameters, drawn from prior/proposal
    given_median <- paras[1]
    given_first_quartile <- paras[2]
    gamma <- paras[3]
    delta <- paras[4]

    # Recalculate the parameters
    params <- calculate_weibull_parameters(given_median, given_first_quartile, delta, gamma)
    alpha <- params$alpha
    beta <- params$beta

    # Set Initial Values
    MAF <- 0.001
    geno_freq <- c(1-MAF,MAF)
    trans <- trans_monogenic(n_alleles = 2)
    baselineRisk <- calculate_lifetime_risk(cancer_type, "SEER", data = PanelPRODatabase)

    penet <- t(sapply(1:nrow(combined_df), function(i) penet.fn(i, combined_df, alpha, beta, gamma, delta, baselineRisk)))
    loglik <- pedigree_loglikelihood(combined_df, geno_freq, trans, penet, ncores = 1)

    return(-loglik)
}
