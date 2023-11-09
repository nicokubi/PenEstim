# Example usage:

# Simulating penetrance data
set.seed(23)  # for reproducibility

# Define the shape parameter (kappa)
alpha <- 1

# Calculate the scale parameter (lambda) for a median of 60
median_age <- 80
beta <- median_age / (log(2)^(1 / alpha))

# Generate ages from 0 to 100 in 5-year increments
ages <- seq(20, 90, by = 5)

# Calculate the cumulative distribution function (CDF) for the Weibull distribution
penetrance_probabilities <- pweibull(ages, shape = alpha, scale = beta)

# Create a data frame
penetrance_df <- data.frame(age = ages, penetrance_prob = penetrance_probabilities)

# Simulating at-risk data (assuming a decreasing number at risk with age)
total_samples <- seq(1000, 500, length.out=length(ages))
at_risk <- total_samples - cumsum(rpois(length(ages), lambda=20))  # some events occurring at each age
samples_table <- data.frame(age = ages, total_samples = total_samples, at_risk = at_risk)


print(penetrance_df)
print(samples_table)
proposal_distributions = create_distributions(penetrance_df,sample_size =1000)
print(proposal_distributions$asymptote_distribution(1))
print(proposal_distributions$shift_distribution(1))
print(proposal_distributions$median_distribution(1))
print(proposal_distributions$first_quartile_distribution(1))

# Create a list to store the results
results_list <- list()

# Run the functions 100 times
for (i in 1:10000) {
  proposal_distributions <- create_distributions()
  results_list[[i]] <- list(
    asymptote = proposal_distributions$asymptote_distribution(1),
    shift = proposal_distributions$shift_distribution(1),
    median = proposal_distributions$median_distribution(1),
    first_quartile = proposal_distributions$first_quartile_distribution(1)
  )
}

# Create 4 histograms for the results
par(mfrow = c(2, 2)) # Arrange the plots in a 2x2 grid

hist(unlist(sapply(results_list, function(x) x$asymptote)), main = "Asymptote Distribution")
hist(unlist(sapply(results_list, function(x) x$shift)), main = "Shift Distribution")
hist(unlist(sapply(results_list, function(x) x$median)), main = "Median Distribution")
hist(unlist(sapply(results_list, function(x) x$first_quartile)), main = "First Quartile Distribution")

# Reset the plotting layout
par(mfrow = c(1, 1))



# Output
library(parallel)
library(PPP)

load("/Users/nicolaskubista/Documents/Master Statistics/Master Thesis/Code/Submission/Simulated Families/simFamilies_C_1000_nocen_selected.RData")

# Example call of PenEstim with proposal_fns and proposal_params.
out7 <- PenEstim(
data=simFamilies_C_1000_nocen_selected, cancer_type = "Breast",
gene_input = "BRCA1", n_chains = 2, n_iter_per_chain = 2,
proposal_distributions = proposal_distributions, density_plots = FALSE,
trace_plots = FALSE
)
