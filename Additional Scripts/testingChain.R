# Data
dat <- simfamilies_sim4.3

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

