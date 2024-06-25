library(clipp)
library(survival)
library(plyr) # need to load plyr before dplyr
library(truncnorm)
library(PanelPRO)
library(tidyverse)
library(stringr)
library(PedUtils)
library(survival)
library(MASS)
library(profvis)
library(survminer)
library(ggridges)
library(ggplot2)
library(dplyr)

# Data
dat <- load("/Users/nicolaskubista/Dropbox (Partners HealthCare)/CCGCRN Hispanic Cohort Data/PenEstim/Data/carrierProbandFamilies_cohPedigree_MLH1.RData")
dat <- carrierProbandFamilies_cohPedigree_MLH1


for (i in seq_along(dat)) {
  if ("ID" %in% colnames(dat[[i]])) {
    colnames(dat[[i]])[colnames(dat[[i]]) == "PedigreeID"] <- "FamilyID"
  }
}

for (i in seq_along(dat)) {
  # Add a new column "PedigreeID" with the list number
  dat[[i]]$PedigreeID <- i
}

# Change "isAffCOL" to "isAff" if "isAffCOL" is a column
for (i in seq_along(dat)) {
  if ("isAffCOL" %in% colnames(dat[[i]])) {
    colnames(dat[[i]])[colnames(dat[[i]]) == "isAffCOL"] <- "isAff"
  }
}

for (i in seq_along(dat)) {
  if ("AgeCOL" %in% colnames(dat[[i]])) {
    colnames(dat[[i]])[colnames(dat[[i]]) == "AgeCOL"] <- "Age"
  }
}

for (i in seq_along(dat)) {
  if ("MLH1" %in% colnames(dat[[i]])) {
    colnames(dat[[i]])[colnames(dat[[i]]) == "MLH1"] <- "geno"
  }
}


data <- do.call(rbind, lapply(dat, transformDF))

head(data)
data

# Exploring different priors
prior_params <- list(
  asymptote = list(g1 = 1, g2 = 1),
  threshold = list(min = 5, max = 30),
  median = list(m1 = 2, m2 = 2),
  first_quartile = list(q1 = 6, q2 = 3)
)


# Run Estimation procedure with default prior setting 
# Main Estimation for Female
system.time(out_sim_COL <- PenEstim(
    data = dat, n_chains = 1, n_iter_per_chain = 10, af = PPP::PanelPRODatabase$AlleleFrequency[paste0("MLH1", "_anyPV"), "nonAJ"],
 ageImputation = TRUE,
    prior_params = prior_params, burn_in = 0.1, median_max = TRUE
))

PenEstim()



# Prop 
prop <- makePriors(
  data = NULL,
  sample_size = NULL,
  cancer = "Breast",
  ratio = NULL,
  prior_params = prior_params,
  risk_proportion = risk_proportion_default
)


set.seed(2)
out25 <-mhChain(
  seed =1, n_iter=10, burn_in = 0.1, chain_id=1, data=data, ncores = 4,
  max_age=94, af = PPP::PanelPRODatabase$AlleleFrequency[paste0("MLH1", "_anyPV"), "nonAJ"],
  baseline_data =  baseline_data_default,
  prior_distributions = prop, removeProband = FALSE,
  median_max = TRUE, max_penetrance = 1, ageImputation = FALSE, BaselineNC  = TRUE, var = c(0.1,0.1,2,2,5,5,5,5))

out25
$out22$summary_stats
plot_traceSingle(out22$results[[1]])
plot_penetrance_sim(out_sim_4$combined_chains, prob = 0.95, max_age = 94,"Female",2.5,50,0,1)

