# PPP 1.1.0

## Major changes
- Manual installation of the `cbcrisk` v2.0 package is no longer required when installing `PPP`. Instead, the first time the user runs the `PPP()` function (or one of its model specific wrapper functions such as `BRCAPRO()`), `cbcrisk` will be automatically installed from GitHub, if it is not already.

## Minor changes
- The intervention columns `riskmod` and `interAge`, which stored information on the prophylactic surgeries (mastectomies, oophorectomies, and hysterectomies) in list type format, made sharing .csv files of PPP pedigrees difficult due to the fact that .csv files cannot easily handle list type columns. Therefore, PPP now also accepts pedigrees with six different numeric columns to replace the riskmod and interAge list type columns. The naming conventions for the six new columns are `riskmodX` (3 columns) and `interAgeX` (3 columns), where "X" is one of the 3 abbreviated surgery names `"Mast"`, `"Ooph"`, or `"Hyst"`. See the README for information on how to encode these columns.

## Bug fixes
- None

# PPP 1.0.0

## Major changes
- Contralateral breast cancer (CBC) has been added to PPP 1.0.0. This requires the user to install the non-CRAN R package cbcrisk manually, prior to installing PPP 1.0.0. Due to the new dependency, if a user updates to PPP 1.0.0, their calls to PPP functions may result in errors if they do not install cbcrisk. Additional and optional columns related to breast cancer are now accepted in the input pedigree to support the addition of CBC. See the README file for information on how to install cbcrisk and the additional columns.
- Updated ATM and PALB2 female breast cancer penetrances based on work from Lakshika Ruberu and Swati Biswas.

## Minor changes
- Relatives in the input pedigree who are disconnected are now removed from the pedigree before analysis and a warning message informs the user of which relatives were disconnected.

## Bug fixes
- A small but fatal bug for the edge case where somebody wants to get net future risk for the proband's CurAge+1. This bug unintentionally dropped dimensions in the cancer penetrance when subsetting for a single age like that results in some mischief.


# PPP 0.3.0

## Major changes
- APC, BMPR1A, and cervical cancer are not supported and have been removed from `PanelPRODatabase`. `PPP-22` replaces `PPP-24` as the largest model. 
- Our age imputation procedure imputes missing current ages, followed by missing cancer affectation ages. If a missing cancer age cannot be sampled such that it satisfies its age bounds and the upper bound set by the current age, this imputation run will be discarded and re-attempted, with a warning message issued. If the desired number of runs (`iterations`, default `20`) cannot be imputed within `max.iter.tries` (default `5*iterations`) attempts, an error will be raised. 
- Missing current ages are now imputed based on age distributions from the literature, instead of truncated normal distributions. 
- When any of the allele frequencies for an ancestry other than `"nonAJ"` differ from the `"nonAJ"` allele frequencies (for a given model specification), the noncarrier penetrances for this ancestry gets set to the SEER penetrances for each cancer in the built database. 

## Minor changes
- Documentation notes that Asian-Pacific Islanders should be classified as `"Asian"` in the `race` column of the user pedigree. 
- Ages between `0` and `1` will be rounded up to `1`; other decimal ages will be rounded to the nearest integer. In both cases, an `AgeRounding` message will be printed. 
- The visualization function `visRisk` now has an option to return the resulting Plotly object instead of plotting it. Checks for probands without carrier probability or future risk estimates (due to the model not specifying genes or cancers) are also implemented. 
- A `random.seed` can now be specified for `checkFam`, such that the age imputations will match those in `PPP` when the same random seed is used. 
- New internal function `checkGeneCancerAssociations` automates checking for and visualizing gene-cancer associations. 
- Master function documentation clarifies that ancestry information is used to select allele frequencies for BRCA1 and BRCA2, specifically. 
- The `use.mult.variants` argument has been removed from the master function, since we currently do not support multiple variants. The documentation also notes that we output both heterozygous and heterozygous carrier probabilities when the database information is available (currently only MUTYH; references for this added). 
- The `checkMating` `checkFam` check that none of a given individual's relatives overlap with the individual's spouse's relatives has been expanded to include more spouses (i.e. relatives by marriage) in the two relative groups. This still doesn't cover all possible loop scenarios, but it helps catch more of them. 
- Removed the deprecated `normalize` option from the peeling-paring function, since we now always normalize the carrier probabilities and do not allow the user to specify otherwise. 
- Added clarification (with examples) on how inheritance probabilities should be interpreted to the `.calInheritanceProb` documentation. 
- When ages were imputed, warning messages will be printed that give each proband's lower and upper bounds for the probability of carrying any pathogenic variant and encourage the user to provide more age information when possible. 

## Bug fixes
- Fixed `buildDatabase` error that occurs when no cancers are passed in and the function attempts to coerce the (empty) `cancers` vector to title-case. 
- Fixed `ageImpute` bugs involving inconsistent encoding of the `Sex` variable and relative information slots. 
- If the lower and upper age bounds get swapped during `ageImpute` (i.e. because the lower bound was greater than the upper bound), the new bounds are now ensured to be valid (i.e. the new upper bound must be below `MAXAGE`). 
- `checkFam` now raises an error when a `MotherID` or `FatherID` in the pedigree has an invalid sex. 
- Interpretation of biomarker testing information in the documentation for `PanelPRODatabase` has been corrected; they should be the probabilities of the marker given genotype, not the other way around. 
- `checkFam` check that drops biomarker testing for those unaffected for the relevant cancer has been fixed. 
- Bug in cancer penetrance for ovarian cancer and BRCA2 has been fixed. 
- Fixed some small typos in documentation and code comments. 


# PPP 0.2.0

## Major changes
- The PPP-24 model specification has been updated to include 18 instead of 11 cancers: Brain, Breast, Cervical, Colorectal, Endometrial, Gastric, Hepatobiliary, Kidney, Leukemia, Melanoma, Ovarian, Osteosarcoma, Pancreas, Prostate, Small Intestine, Soft Tissue Sarcoma, Thyroid, and Urinary Bladder. 
- Eight additional test pedigrees have been added: `test_fam_5` through `test_fam_12`. 
- `PanelPRODatabase` has been updated to remove the association between male breast cancer and BRCA1. 
- When an individual in the pedigree has an age of cancer diagnosis that is greater than their current age, `checkFam` now raises an error and suggests that if the user cannot resolve the age inputs, the maximum cancer age can be entered as the current age. 

## Minor changes
- When one or more of the carrier probabilities is `0`, germline testing results were incorporated, and the default sensitivities/specifities of `1` were used, a warning will be issued that suggests a more accurate estimate may be achieved with user-specified sensitivities and specificities.
- The new `ignore.proband.germ` argument in `PPP` allows users the option to ignore proband germline testing results (if provided) from being considered by the model. 
- When the current age of one or more living probands is equal to `PPP:::MAXAGE`, no future risk probabilities will be calculated for these individuals. 

## Bug fixes
- A bug that scrambled the model parameters for the breast tumor markers has been fixed. 
- Imputed ages can no longer exceed the legal range of ages (between 1 and `PPP:::MAXAGE`). 


# PPP 0.1.0

* This is the first release of PPP.
