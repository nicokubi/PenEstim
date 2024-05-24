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
dat <- load("/Users/nicolaskubista/Dropbox (Partners HealthCare)/CCGCRN Hispanic Cohort Data/PenEstim/Data/carrierProbandFamilies_cohPedigree_MLH1_ages.RData")
dat <- carrierProbandFamilies_cohPedigree_MLH1_ages

# Data Prep
for (i in seq_along(dat)) {
  if ("ID" %in% colnames(dat[[i]])) {
    colnames(dat[[i]])[colnames(dat[[i]]) == "PedigreeID"] <- "FamilyID"
  }
}

#  Data Prep
for (i in seq_along(dat)) {
  if ("ID" %in% colnames(dat[[i]])) {
    colnames(dat[[i]])[colnames(dat[[i]]) == "ID"] <- "SubjectID"
  }
}

for (i in seq_along(dat)) {
  # Add a new column "PedigreeID" with the list number
  dat[[i]]$PedigreeID <- i
}

data <- do.call(rbind, lapply(dat, transformDF,
                              cancer_type = "Breast",
                              gene = "BRCA1"
))

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
system.time(out_sim_COL <- PenEstim_v7(
    data = dat,
    cancer_type = "Colorectal", gene_input = "MLH1", n_chains = 1, n_iter_per_chain = 10, 
    prior_params = prior_params, af = 0.1, burn_in = 0.1, median_max = TRUE, priors = prior_params
))

system.time(out_sim_COL <- PenEstim_v10(
  data = dat,
  cancer_type = "Colorectal", gene_input = "MLH1", ncores = 2, n_iter_per_chain = 5,
  prior_params = prior_params, af = 0.1, burn_in = 0.1, median_max = TRUE, priors = prior_params
))

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
out25 <-mhChain_vMWG(
  seed =1, n_iter=10, burn_in = 0.1, chain_id=1, data=data,
  max_age=94, db=PanelPRODatabase,
  prior_distributions = prop, cancer_type = "Breast", gene_input = "BRCA1", af = 0.1,
  median_max = TRUE, max_penetrance = 1, homozygote = TRUE, SeerNC = TRUE, priors = prior_params, var = c(0.1,0.1,2,2,5,5,5,5))

out25
$out22$summary_stats
plot_traceSingle(out22$results[[1]])
plot_penetrance_sim(out_sim_4$combined_chains, prob = 0.95, max_age = 94,"Female",2.5,50,0,1)

