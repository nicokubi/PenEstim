# Load Libraries
library(clipp)
library(PPP)
library(stats4)
library(dplyr)

# Data Preparation and Transformation
dat <- prepAges(carrierProbandFamilies_cohPedigree_MSH6_COL)

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

combined_df <- do.call(rbind, lapply(dat, transform_df))

# Set Initial Values
MAF <- 0.1
geno_freq <- geno_freq_monogenic(p_alleles = c(1 - MAF, MAF))
trans <- trans_monogenic(n_alleles = 2)

# Function to Calculate Lifetime Risk
calculate_lifetime_risk <- function(cancer_type, gene, race = "All_Races", sex = "Female", type = "Crude", data = PanelPRODatabase) {
    # Find indices for respective attributes
    dim_names <- attr(data$Penetrance, "dimnames")
    indices <- sapply(list(gene, cancer_type, race, sex, type), \(x) which(dim_names[[x]] == x))
    lifetime_risk <- data$Penetrance[indices]
    return(lifetime_risk)
}

# Penetrance Functions
penet.fn <- function(i, data, alpha, beta, gamma, delta) {
    max_age <- 94 # Assuming "inc" is defined elsewhere
    if (data$age[i] == 0) {
        penet.i <- c(1, 1, 1) # Assuming people aged 0 years are all unaffected
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
        penet.i <- c(1 - nc.pen, 1 - c.pen, 1 - c.pen) # if person is not affected
        if (data$aff[i] == 1) penet.i <- c(nc.pen, c.pen, c.pen)
    }
    if (data$geno[i] == "1/1") penet.i[-1] <- 0
    if (data$geno[i] == "1/2") penet.i[-2] <- 0
    if (data$geno[i] == "2/2") penet.i[-2] <- 0
    return(penet.i)
}

# Example Usage and Calculations
baselineRisk <- calculate_lifetime_risk("Breast", "SEER")
alpha1 <- 1.1
beta1 <- 20
delta1 <- 0
gamma1 <- 1

penet <- t(sapply(1:nrow(combined_df), penet.fn, combined_df, alpha = alpha1, beta = beta1, gamma = gamma1, delta = delta1))
loglik <- pedigree_loglikelihood(combined_df, geno_freq, trans, penet, ncores = 2)
