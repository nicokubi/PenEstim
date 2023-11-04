# Example usage:

# Functions to generate initial and proposed parameter values.
proposal_dist <- list(
  asymptote = function(params) total_prob + rbeta(1, params$g1, params$g2) * (1 - params$p0),
  shift = function(params) runif(1, params$min, params$max),
  median = function(params,baseline_mid) rbeta(1, params$m1, params$m2) * baseline_mid+params$eps + params$shift,
  rbeta(1,m1,m2)*(baseline_mid+eps-shift_presroposal) + shift_proposal
  quartile = function(params) rbeta(1, params$q1, params$q2) * params$scale + params$shift
)

# Parameters for the proposal functions.
proposal_params <- list(
  asymptote = list(g1 = 9, g2 = 1),
  shift = list(min = 0, max = 25),
  median = list(m1 = 2, m2 = 2, scale = 10, shift = 5), # You should adjust scale and shift according to your data.
  quartile = list(q1 = 6, q2 = 3, scale = 10, shift = 5) # You should adjust scale and shift according to your data.

)

proposal_params <- list(
  m1 = 2, m2 = 2.0,
  g1 = 0.5, g2 = 2.0,
  eps = 5,
  shift_prior_min = 0.0, shift_prior_max = 25,
  q1 = 0.5, q2 = 2.0,
  p0 = 0.1
)

# Example for the automatic prior elictiation

# Example usage:
# Simulating penetrance data
set.seed(42)  # for reproducibility

ages <- seq(20, 80, by=5)  # ages from 20 to 80 in 5 year intervals
a <- 0.2  # controls the steepness of the curve
b <- 50  # median age of penetrance
penetrance_prob <- 1 / (1 + exp(-(a * (ages - b))))

# Add some noise to the penetrance probabilities
penetrance_prob <- penetrance_prob + rnorm(length(ages), 0, 0.05)
penetrance_prob <- pmin(pmax(penetrance_prob, 0), 1)  # Ensure values are between 0 and 1

# Create the dataframe
dataframe <- data.frame(age = ages, penetrance_prob = penetrance_prob)

# Simulating at-risk data (assuming a decreasing number at risk with age)
total_samples <- seq(1000, 500, length.out=length(ages))
at_risk <- total_samples - cumsum(rpois(length(ages), lambda=20))  # some events occurring at each age
samples_table <- data.frame(age = ages, total_samples = total_samples, at_risk = at_risk)

proposal_distributions =create_distributions()
print(proposal_distributions$asymptote_distribution(1))
print(proposal_distributions$shift_distribution(1))
print(proposal_distributions$median_distribution(1))
print(proposal_distributions$first_quartile_distribution(1))

# Output
library(parallel)
library(PPP)

# Example call of PenEstim with proposal_fns and proposal_params.
out6 <- PenEstim(simFamilies_C_1000_nocen_selected, 4, 2,
proposal_distributions = proposal_distributions, density_plots = FALSE,
trace_plots = FALSE
)