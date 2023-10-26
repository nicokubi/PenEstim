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


# Example call of PenEstim with proposal_fns and proposal_params.
out3 <- PenEstim(simFamilies_C_1000_nocen_selected, 4, 3,proposal_params=proposal_params,density_plots=FALSE, trace_plots=FALSE)
out3