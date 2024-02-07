library(PPP)
library(parallel)
library(dplyr)
library(stats4)
library(clipp)
# Comparing PanelPRO and Clipp
# The idea is to compare the likelihood calculation with PanelPRO and clipp
# Should be somewhat in the same magniute

 #  Create the prior distributions
 prop <- makePriors(
     data = NULL,
     sample_size = NULL,
     cancer = "Breast",
     ratio = NULL,
     prior_params = prior_params_default,
     risk_proportion = risk_proportion_default
 )
# Using the approach without clipp

# Import the results from the other machine 


# Compare 
ll_pp <- out_PP$results$loglikelihood_current
ll_clipp <- out_clipp$results$loglikelihood_current

ll_pp <- out_PP$results$loglikelihood_proposal 
ll_clipp <- out_clipp$results$loglikelihood_proposal 

# Comparing the output of the likelihood calculation just for one set of parameters
# Compare the mhChain fuction with one set of parameters, once with PP and once with Clipp

# set seed
seed <- 2024

# Set the parameters
paras <- c(60,40,20,0.9)

# prep data for PP
load("/Users/nicolaskubista/Documents/Master Statistics/Master Thesis/Code/Submission/Simulated Families/simFamilies_C_1000_nocen_selected.RData")
datPP <- prepAges(simFamilies_C_1000_nocen_selected,removeProband = FALSE)

# set other inputs
n_iter <- 1
chain_id <- 1
data <- dat2
max_age <- 94
db <- PanelPRODatabase
prior_distributions <- prop
cancer_type <- "Breast"
gene_input <- "BRCA1"
max_penetrance <- 1

PanelPRODatabase$AlleleFrequency

str(datPP)

# Calculation using clipp
ll_out <- mhLogLikelihood(paras, datPP[10], max_age, db, cancer_type, gene_input)
print(ll_out)
head(datPP)

# Approach with clipp
# Using clipp
load("/Users/nicolaskubista/Documents/Master Statistics/Master Thesis/Code/Submission/Simulated Families/simFamilies_C_1000_nocen_selected.RData")
str(simFamilies_C_1000_nocen_selected)
length(simFamilies_C_1000_nocen_selected)


#  Data Prep

for (i in seq_along(simFamilies_C_1000_nocen_selected)) {
    # Add a new column "PedigreeID" with the list number
    simFamilies_C_1000_nocen_selected[[i]]$PedigreeID <- i
}
for (i in seq_along(simFamilies_C_1000_nocen_selected)) {
    if ("ID" %in% colnames(simFamilies_C_1000_nocen_selected[[i]])) {
        colnames(simFamilies_C_1000_nocen_selected[[i]])[colnames(simFamilies_C_1000_nocen_selected[[i]]) == "PedigreeID"] <- "PedigreeID"
    }
}

#  Data Prep
for (i in seq_along(simFamilies_C_1000_nocen_selected)) {
    if ("ID" %in% colnames(simFamilies_C_1000_nocen_selected[[i]])) {
        colnames(simFamilies_C_1000_nocen_selected[[i]])[colnames(simFamilies_C_1000_nocen_selected[[i]]) == "ID"] <- "SubjectID"
    }
}

# Apply the prepAges function to treat missing data
dat <- prepAges(simFamilies_C_1000_nocen_selected, removeProband = FALSE)

# Apply the transformation to adjust the format for the clipp package
dat <- do.call(rbind, lapply(dat, transformDF))
dat10 <- dat[dat$family == 10, ]
dat1

# use same af
af <- PanelPRODatabase$AlleleFrequency[3,2]
af
llClipp <- mhLogLikelihood_clipp(paras, dat10, max_age, cancer_type, db, af)
ll_out

PPP(pedigree = datPP[[1]], genes = c("BRCA1"), cancers = "Breast",
    database = PanelPRODatabase, impute.missing.ages = FALSE)$posterior.prob[[1]]

datPP
log(2.16e-30)
