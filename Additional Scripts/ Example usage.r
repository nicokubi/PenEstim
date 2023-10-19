# Example usage:

# Functions to generate initial and proposed parameter values.
proposal_dist <- list(
  asymptote = function(params) params$p0 + rbeta(1, params$g1, params$g2) * (1 - params$p0),
  shift = function(params) runif(1, params$min, params$max),
  median = function(params) rbeta(1, params$m1, params$m2) * params$scale + params$shift,
  quartile = function(params) rbeta(1, params$q1, params$q2) * params$scale + params$shift
)

# Parameters for the proposal functions.
proposal_params <- list(
  asymptote = list(p0 = 0.15, g1 = 9, g2 = 1),
  shift = list(min = 0, max = 25),
  median = list(m1 = 2, m2 = 2, scale = 10, shift = 5), # You should adjust scale and shift according to your data.
  quartile = list(q1 = 6, q2 = 3, scale = 10, shift = 5) # You should adjust scale and shift according to your data.
)

# Example call of PenEstim with proposal_fns and proposal_params.
results <- PenEstim(data, 3, 1000, 200, 94, proposal_fns, proposal_params)
