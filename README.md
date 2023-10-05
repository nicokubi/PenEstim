# hvd-thesis
Work related to the Master Thesis at Harvard and DFCI 

## PenEstim
The PenEstim Package works as an extension to the PanelPRO R Package developed by the BayesMendel Lab. The pacakge is used to estimate age-specific penetrance for complex family-based data in a format compatible with PanelPRO.

Additional Notes:
- Dependencies: `PPP` (adapted version of PanelPRO), `stats`
- Applies independent Metropolis Hastings algorithm to draw from the posterior distribution to estimate the parameters of a four-parameter Weibull distibution
- Assumes that Proposal = Prior distribution
- Parameters of the prior/proposal can be adapted from default values
 

## Simulated Families
The following three types of families were used for the simulation study. All families were simulated without censoring and with different parameters for the four-parameter Weibull distribution: 
- simFamilies_A_475_nocen (n=475)
- simFamilies_B_1000_nocen (n=488)
- simFamilies_C_1000_nocen_selected (n=60)

## Additional Scripts
Additiona Scripts used as part of this thesis: 
- `LikEval.R` was used for the single-parameter likelihood evaluation used in the application of the MLE approach.
- `MCMCplots.R` was used to generate plots to analyze the results from the MCMC approach.


