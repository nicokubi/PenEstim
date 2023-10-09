# Thesis: Methods and Software for Penetrance Estimation in Complex Family-based Studies
Work related to the Master Thesis at the BayesMendel Lab at Harvard University under the supervision of Prof. Dr. Giovanni Parmigiani (Harvard University), Dr. Danielle Braun (Harvard University), and Dr. Markus Kalisch (ETH Zurich), 

## PenEstim
The PenEstim Package works as an extension to the PanelPRO R Package developed by the BayesMendel Lab. The package is used to estimate age-specific penetrance for complex family-based data in a format compatible with PanelPRO.

Additional Notes:
- Dependencies: `PPP` (adapted version of PanelPRO), `stats`. PPP needs to be installed before loading the PenEstim package. 
- Applies independent Metropolis Hastings algorithm to draw from the posterior distribution to estimate the parameters of a four-parameter Weibull distribution
- Assumes that Proposal = Prior distribution
- Parameters of the prior/proposal can be adapted from default values
 

## Simulated Families
The following three types of families were used for the simulation study. All families were simulated without censoring and with different parameters for the four-parameter Weibull distribution: 
- simFamilies_A_475_nocen (n=475)
- simFamilies_B_1000_nocen (n=488)
- simFamilies_C_1000_nocen_selected (n=60)

## Additional Scripts
Additional Scripts used as part of this thesis: 
- `LikEval.R` was used for the single-parameter likelihood evaluation used in the application of the MLE approach.
- `MCMCplots.R` was used to generate plots to analyze the results from the MCMC approach.
- `peelingParing_adjusted.cpp` contains the adjusted cpp routine for the likelihood evaluation based on the Peeling and Paring algorithm. 
- `Plot_Families.R` was used to generate plots to display the empirical penetrance curves from the simulated families. 
- `Sim_Families.R` was used to generate the Families for the simulation study. 



