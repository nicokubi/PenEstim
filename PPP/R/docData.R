#' Built-in PPP database object
#'
#' Built-in parameters for `PPP` models. This object stores the default 
#' values, but the user can set their own parameter values or create their own 
#' database object with the same format. 
#' 
#' @details 
#' A list of 7 elements: 
#' * `Penetrance`: An array that provides the age-specific probabilities of 
#' developing a certain cancer, given genotype, race, and sex. The array has 6 
#' levels corresponding to the cancer, gene, race, sex, age, and type of 
#' penetrance (net or crude). See `str(PanelPRODatabase$Penetrance)` for the 
#' structure of the array. For example, 
#' `PanelPRODatabase$Penetrance["Breast", "BRCA1_hetero_anyPV", "Asian", "Female", "40", "Net"]`
#' is the net probability that an Asian woman with any heterozygous pathogenic 
#' BRCA1 variant develops breast cancer at age 40. For now, mutations in all 
#' genes are considered as heterozygous mutations except for mutations in MUTYH, 
#' which are allowed to be either heterozygous or homozygous.
#' * `AlleleFrequency`: An array that provides the minor allele frequencies for 
#' mutations in each gene variant and ancestry type (Ashkenazi Jewish, 
#' non-Ashkenazi Jewish, and Italian). The array has 2 levels corresponding to 
#' gene and ancestry. See `str(PanelPRODatabase$AlleleFrequency)` for the 
#' structure of the array. For example, 
#' `PanelPRODatabase$AlleleFrequency["BRCA1_anyPV", "AJ"]` is the allele 
#' frequency of any pathogenic BRCA1 variant for an Ashkenazi Jewish individual.
#' * `DOC`: An array that provides the age-specific probabilities of death by 
#' other causes (other than the specified cancer), given genotype, race, and 
#' sex. The array has 5 levels corresponding to the cancer, gene, race, sex, 
#' and age. See `str(PanelPRODatabase$DOC)` for the structure of the array. For 
#' example, 
#' `PanelPRODatabase$DOC["Breast", "BRCA1_hetero_anyPV", "Asian", "Female", "40"]` 
#' is the probability that a 40-year-old Asian women with any heterozygous 
#' pathogenic BRCA1 variant who is currently alive dies from a cause other than 
#' breast cancer. This is used for estimating crude future risk.
#' * `Riskmod`: An array of parameters (either relative risks or hazard ratios) 
#' that modify risk at the cancer penetrance level. Risk modifiers are assumed 
#' to be independent of each other. The array has 6 levels corresponding to the 
#' cancer, gene, type of intervention (mastectomy, hysterectomy, oophorectomy), 
#' sex, age of intervention, and type of parameter (relative risk ratio or 
#' hazard ratio). See `str(PanelPRODatabase$Riskmod)` for the structure of the 
#' array. For example, 
#' `PanelPRODatabase$Riskmod["Breast", "BRCA1_hetero_anyPV", "Mastectomy", "Female", "40", "HR"]`
#' is the ratio of the hazard of breast cancer for a woman with any 
#' heterozygous pathogenic BRCA1 variant who had a mastectomy at age 40 
#' compared to the hazard for such a woman who did not have a mastectomy at age 40.
#' * `GermlineTesting`: An array that provides the sensitivities and 
#' specificities for germline testing, used to modify the likelihood. The array 
#' has 2 levels corresponding to the gene and measure (sensitivity or 
#' specificity). See `str(PanelPRODatabase$GermlineTesting)` for the structure 
#' of the array. For example, 
#' `PanelPRODatabase$GermlineTesting["BRCA1_hetero_anyPV", "sensitivity"]` 
#' is the sensitivitiy for detecting BRCA1. 
#' * `BiomarkerTesting`: A list of parameters for tumor biomarker testing, 
#' where each component is an array corresponding to a different cancer. Each 
#' array represents a configuration of the tumor marker testing results (either 
#' "NoTest", "Neg", or "Pos") and the genotypes (0 or 1 for each gene that is 
#' related to the cancer). The parameters represent the probability of the 
#' tumor marker configuration given genotype. See 
#' `str(PanelPRODatabase$BiomarkerTesting)` for the structure of the list and 
#' arrays. For example, 
#' `PanelPRODatabase$BiomarkerTesting$Breast[ER="Pos", CK5.6="Pos", CK14="Pos", PR="Neg", HER2="Neg", BRCA1_hetero_anyPV="1", BRCA2_hetero_anyPV="0"]`
#' is the probability of being a BRCA1 pathogenic variant carrier but not a 
#' BRCA2 pathogenic variant carrier given that the individual tested positive 
#' for ER, CK5.6, and CK14 and negative for PR and HER2.
#' * `Contralateral`: An array that provides the age-specific probabilities of 
#' developing contralateral breast cancer for carriers of pathogenic germline 
#' mutations of CBC associated genes and who have already had a first primary 
#' breast cancer in one breast, given genotype, sex, and age of the first breast 
#' cancer diagnosis. The array has 5 levels corresponding to the gene, sex, 
#' age at first breast cancer, type of penetrance (net or crude), and age. 
#' See `str(PanelPRODatabase$Contralateral)` for the structure of the array. 
#' For example, 
#' `PanelPRODatabase$Penetrance["BRCA1_hetero_anyPV", "Female", "40", "Net", "50"]`
#' is the net probability that a woman with any heterozygous pathogenic 
#' BRCA1 variant, who had their first breast cancer at age 40, develops 
#' contralateral breast cancer at age 50. For now, mutations in the CBC 
#' associated genes are considered as heterozygous mutations.
#' 
#' @seealso \code{\link{MaleCBCNoncarrierArray}}
#' 
#' @references 
#' Antoniou, A. C., et al. "A comprehensive model for familial breast cancer incorporating BRCA1, BRCA2 and other genes." British journal of cancer 86.1 (2002): 76-83.
#' 
#' Begg, Colin B., et al. "Lifetime risk of melanoma in CDKN2A mutation carriers in a population-based sample." Journal of the National Cancer Institute 97.20 (2005): 1507-1515.
#' 
#' Berwick, Marianne, et al. "The prevalence of CDKN2A germ-line mutations and relative risk for cutaneous malignant melanoma: an international population-based study." Cancer Epidemiology and Prevention Biomarkers 15.8 (2006): 1520-1525.
#' 
#' Bishop, D. Timothy, et al. "Geographical variation in the penetrance of CDKN2A mutations for melanoma." Journal of the National Cancer Institute 94.12 (2002): 894-903.
#' 
#' Braun, Danielle, et al. "A clinical decision support tool to predict cancer risk for commonly tested cancer-related germline mutations." Journal of genetic counseling 27.5 (2018): 1187-1199.
#' 
#' Chen, Jinbo, et al. "Penetrance of Breast and Ovarian Cancer in Women Who Carry a BRCA1/2 Mutation and Do Not Use Risk-Reducing Salpingo-Oophorectomy: An Updated Meta-Analysis." JNCI cancer spectrum 4.4 (2020): pkaa029.
#' 
#' Chen, Sining, et al. "BayesMendel: an R environment for Mendelian risk prediction." Statistical applications in genetics and molecular biology 3.1 (2004).
#' 
#' Chen, Sining, et al. "Prediction of germline mutations and cancer risk in the Lynch syndrome." Jama 296.12 (2006): 1479-1487.
#' 
#' Chowdhury, M., et al. “A Model for Individualized Risk Prediction of Contralateral Breast Cancer”, Breast Cancer Research and Treatment 161 (2017): 153-160.
#' 
#' Chu, Haitao, Sining Chen, and Thomas A. Louis. "Random effects models in a meta-analysis of the accuracy of two diagnostic tests without a gold standard." (2007).
#' 
#' DevCan: Probability of Developing or Dying of Cancer Software, Version 6.7 Surveillance Research Program, Statistical Methodology and Applications, National Cancer Institute, 2012.  http://surveillance.cancer.gov/devcan/
#' 
#' Felton, K. E. A., D. M. Gilchrist, and S. E. Andrew. "Constitutive deficiency in DNA mismatch repair: is it time for Lynch III?." Clinical genetics 71.6 (2007): 499-500.
#' 
#' Guo, Yonghai, et al. "Risk of ipsilateral breast tumor recurrence and contralateral breast cancer in patients with and without TP53 variant in a large series of breast cancer patients." The Breast 65 (2022): 55-60.
#' 
#' Katki, Hormuzd A. "Incorporating medical interventions into carrier probability estimation for genetic counseling." BMC medical genetics 8.1 (2007): 1-11.
#' 
#' Kuchenbaecker, Karoline B., et al. "Risks of breast, ovarian, and contralateral breast cancer for BRCA1 and BRCA2 mutation carriers." Jama 317.23 (2017): 2402-2416.
#' 
#' Lakhani, Sunil R., et al. "The pathology of familial breast cancer: predictive value of immunohistochemical markers estrogen receptor, progesterone receptor, HER-2, and p53 in patients with mutations in BRCA1 and BRCA2." Journal of clinical oncology 20.9 (2002): 2310-2318.
#' 
#' Lakhani, Sunil R., et al. "Prediction of BRCA1 status in patients with breast cancer using estrogen receptor and basal phenotype." Clinical cancer research 11.14 (2005): 5175-5180.
#' 
#' Lee, Andrew J., et al. "Incorporating truncating variants in PALB2, CHEK2, and ATM into the BOADICEA breast cancer risk model." Genetics in Medicine 18.12 (2016): 1190-1198.
#' 
#' Rosenthal, Eric T., et al. "Clinical testing with a panel of 25 genes associated with increased cancer risk results in a significant increase in clinically significant findings across a broad range of cancer histories." Cancer genetics 218 (2017): 58-68.
#' 
#' Sajal, I.H., et al. "CBCRisk-Black: A Personalized Contralateral Breast Cancer Risk Prediction Model for Black Women", Breast Cancer Research and Treatment 194 (2022): 179-186.
#' 
#' Theodoratou, Evropi, et al. "A large-scale meta-analysis to refine colorectal cancer risk estimates associated with MUTYH variants." British journal of cancer 103.12 (2010): 1875-1884.
#' 
#' Wang, Cathy, et al. "Penetrance of colorectal cancer among mismatch repair gene mutation carriers: a meta-analysis." JNCI cancer spectrum 4.5 (2020): pkaa027.
#' 
#' Wang, Wenyi, et al. "Estimating CDKN2A carrier probability and personalizing cancer risk assessments in hereditary melanoma using MelaPRO." Cancer research 70.2 (2010): 552-559.
#' 
#' Win, Aung Ko, et al. "Risk of extracolonic cancers for people with biallelic and monoallelic mutations in MUTYH." International journal of cancer 139.7 (2016): 1557-1563.
#' 
"PanelPRODatabase"

#' MaleCBCNoncarrierArray
#'
#' An array which stores age-specific conditional penetrances for male contralateral 
#' breast cancer (CBC) for individuals without a pathogenic germline mutation.
#' 
#' @details
#' The array is stratified by age at first breast cancer (`FirstBCAge`), 
#' net or crude penetrance type (`PenetType`), and current age (`Age`). The age 
#' ranges for `FirstBCType` and `Age` are from `1` to `94` and the choices for 
#' `Penet Type` are one of `c("Net", "Crude")`.
#' 
#' @seealso \code{\link{PanelPRODatabase}}
#' 
#' @references 
#' Braun D, Yang J, Griffin M, Parmigiani G, Hughes KS. "A Clinical Decision Support Tool to Predict Cancer Risk for Commonly Tested Cancer -Related Germline Mutations." J Genet Couns. 2018 Mar 2.
"MaleCBCNoncarrierArray"

#' test_fam_1
#'
#' A simple but complete example of pedigree data frame with 19 relatives and 
#' family history of breast and ovarian cancer. It also includes germline 
#' testing results for BRCA1 and BRCA2, as well as tumor biomarker testing for 
#' ER, CK5.6, CK14, PR, and HER2. See the \code{\link{PPP}} documentation 
#' for the columns that should be included in a pedigree. 
#'
#' @seealso \code{\link{test_fam_2}}, \code{\link{test_fam_3}}, 
#' \code{\link{test_fam_4}}, \code{\link{test_fam_5}},
#' \code{\link{test_fam_6}}, \code{\link{test_fam_7}},
#' \code{\link{test_fam_7}}, \code{\link{test_fam_9}},
#' \code{\link{test_fam_8}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_1"

#' test_fam_2
#'
#' An example of a pedigree with a more complex family structure. It contains 25 
#' relatives and family history of endometrial, pancreatic, and small intestine 
#' cancer. See the \code{\link{PPP}} documentation for the columns that 
#' should be included in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_3}}, 
#' \code{\link{test_fam_4}}, \code{\link{test_fam_5}},
#' \code{\link{test_fam_6}}, \code{\link{test_fam_7}},
#' \code{\link{test_fam_7}}, \code{\link{test_fam_9}},
#' \code{\link{test_fam_8}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_2"

#' test_fam_3
#' 
#' An example of a pedigree with 50 relatives and family history of 11 cancers, 
#' corresponding to \code{\link{PanelPRO11}}. See the \code{\link{PPP}} 
#' documentation for the columns that should be included in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_4}}, \code{\link{test_fam_5}},
#' \code{\link{test_fam_6}}, \code{\link{test_fam_7}},
#' \code{\link{test_fam_8}}, \code{\link{test_fam_9}},
#' \code{\link{test_fam_10}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_3"

#' test_fam_4
#' 
#' An example of a pedigree with 9 relatives and family history of 7 cancers: 
#' breast, brain, colorectal, prostate, endometrial, small intestine, and 
#' ovarian. It contains no missing information, resulting in exact estimates 
#' for carrier probabilities and future risks. See the \code{\link{PPP}}
#' documentation for the columns that should be included in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_3}}, \code{\link{test_fam_5}},
#' \code{\link{test_fam_6}}, \code{\link{test_fam_7}},
#' \code{\link{test_fam_8}}, \code{\link{test_fam_9}},
#' \code{\link{test_fam_10}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_4"

#' test_fam_5
#' 
#' An example of a pedigree with 17 relatives and family history of breast 
#' cancer and melanoma. It contains no missing information, resulting in exact 
#' estimates for carrier probabilities and future risks. See the 
#' \code{\link{PPP}} documentation for the columns that should be included 
#' in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_3}}, \code{\link{test_fam_4}},
#' \code{\link{test_fam_6}}, \code{\link{test_fam_7}},
#' \code{\link{test_fam_8}}, \code{\link{test_fam_9}},
#' \code{\link{test_fam_10}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_5"

#' test_fam_6
#' 
#' An example of a pedigree with 19 relatives and family history of 4 cancers: 
#' breast, endometrial, pancreatic and melanoma. It contains no missing 
#' information, resulting in exact estimates for carrier probabilities and 
#' future risks. See the \code{\link{PPP}} documentation for the columns 
#' that should be included in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_3}}, \code{\link{test_fam_4}},
#' \code{\link{test_fam_5}}, \code{\link{test_fam_7}},
#' \code{\link{test_fam_8}}, \code{\link{test_fam_9}},
#' \code{\link{test_fam_10}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_6"

#' test_fam_7
#' 
#' An example of a pedigree with 19 relatives and family history of colorectal 
#' and endometrial cancers. It contains no missing information, resulting in 
#' exact estimates for carrier probabilities and future risks. See the 
#' \code{\link{PPP}} documentation for the columns that should be included 
#' in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_3}}, \code{\link{test_fam_4}},
#' \code{\link{test_fam_5}}, \code{\link{test_fam_6}},
#' \code{\link{test_fam_8}}, \code{\link{test_fam_9}},
#' \code{\link{test_fam_10}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_7"

#' test_fam_8
#' 
#' An example of a pedigree with 20 relatives and family history of colorectal 
#' and prostate cancers. It contains no missing information, resulting in exact 
#' estimates for carrier probabilities and future risks. See the 
#' \code{\link{PPP}} documentation for the columns that should be included 
#' in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_3}}, \code{\link{test_fam_4}},
#' \code{\link{test_fam_5}}, \code{\link{test_fam_6}},
#' \code{\link{test_fam_1}}, \code{\link{test_fam_9}},
#' \code{\link{test_fam_10}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_8"

#' test_fam_9
#' 
#' An example of a pedigree with 19 relatives and family history of breast and 
#' prostate cancers. It contains no missing information, resulting in exact 
#' estimates for carrier probabilities and future risks. See the 
#' \code{\link{PPP}} documentation for the columns that should be included 
#' in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_3}}, \code{\link{test_fam_4}},
#' \code{\link{test_fam_5}}, \code{\link{test_fam_6}},
#' \code{\link{test_fam_7}}, \code{\link{test_fam_8}},
#' \code{\link{test_fam_10}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_9"

#' test_fam_10
#' 
#' An example of a pedigree with 21 relatives and family history of breast and 
#' colorectal cancers. It contains no missing information, resulting in exact 
#' estimates for carrier probabilities and future risks. See the 
#' \code{\link{PPP}} documentation for the columns that should be included 
#' in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_3}}, \code{\link{test_fam_4}},
#' \code{\link{test_fam_5}}, \code{\link{test_fam_6}},
#' \code{\link{test_fam_7}}, \code{\link{test_fam_8}},
#' \code{\link{test_fam_9}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_10"

#' test_fam_11
#' 
#' An example of a pedigree with 16 relatives and family history of 5 cancers: 
#' breast, endometrial, ovarian, pancreatic and small intestine. It contains no 
#' missing information, resulting in exact estimates for carrier probabilities 
#' and future risks. See the \code{\link{PPP}} documentation for the 
#' columns that should be included in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_3}}, \code{\link{test_fam_4}},
#' \code{\link{test_fam_5}}, \code{\link{test_fam_6}},
#' \code{\link{test_fam_7}}, \code{\link{test_fam_8}},
#' \code{\link{test_fam_9}}, \code{\link{test_fam_10}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_11"

#' test_fam_12
#' 
#' An example of a pedigree with 21 relatives and family history of 17 cancers, 
#' corresponding to \code{\link{PanelPRO22}}. It contains no missing 
#' information, resulting in exact estimates for carrier probabilities and 
#' future risks. See the \code{\link{PPP}} documentation for the columns 
#' that should be included in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_4}}, \code{\link{test_fam_5}},
#' \code{\link{test_fam_6}}, \code{\link{test_fam_7}},
#' \code{\link{test_fam_8}}, \code{\link{test_fam_9}},
#' \code{\link{test_fam_10}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{test_fam_13}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_12"

#' test_fam_13
#' 
#' This pedigree is a copy of test_fam_5 but with the proband affected with
#' breast cancer at a known age. This allows contralateral breast cancer
#' to be incorporated into the model which results in future risks for this 
#' cancer. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_3}}, \code{\link{test_fam_4}},
#' \code{\link{test_fam_5}}, \code{\link{test_fam_6}},
#' \code{\link{test_fam_7}}, \code{\link{test_fam_8}},
#' \code{\link{test_fam_9}}, \code{\link{test_fam_10}},
#' \code{\link{test_fam_11}}, \code{\link{test_fam_12}},
#' \code{\link{test_fam_14}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_13"

#' test_fam_14
#' 
#' This pedigree is a copy of test_fam_1 but with the proband affected with
#' breast cancer at a known age. This allows contralateral breast cancer
#' to be incorporated into the model which results in future risks for this 
#' cancer. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_3}}, \code{\link{test_fam_4}},
#' \code{\link{test_fam_5}}, \code{\link{test_fam_6}},
#' \code{\link{test_fam_7}}, \code{\link{test_fam_8}},
#' \code{\link{test_fam_9}}, \code{\link{test_fam_10}},
#' \code{\link{test_fam_11}}, \code{\link{test_fam_12}},
#' \code{\link{test_fam_13}}, \code{\link{err_fam_1}}
#' @family family
"test_fam_14"

#' err_fam_1
#'
#' An example of a pedigree which has 'loops' (intermarriages). PPP does 
#' not currently support such pedigrees, so passing it into 
#' \code{\link{PPP}} or \code{\link{checkFam}} will result in an error. It 
#' contains 10 relatives and family history of breast and pancreatic cancer. 
#' See the \code{\link{PPP}} documentation for the columns that should be 
#' included in a pedigree. 
#' @seealso \code{\link{test_fam_1}}, \code{\link{test_fam_2}}, 
#' \code{\link{test_fam_4}}, \code{\link{test_fam_5}},
#' \code{\link{test_fam_6}}, \code{\link{test_fam_7}},
#' \code{\link{test_fam_8}}, \code{\link{test_fam_9}},
#' \code{\link{test_fam_10}}, \code{\link{test_fam_11}},
#' \code{\link{test_fam_12}}, \code{\link{err_fam_1}}
#' @family family
"err_fam_1"


#' Race-specific estimates of 1 / (1 - population attributable risk)
#'
#' @format Data frame of race-specific estimates of 
#' 1 / (1 - population attributable risk), estimated from NHIS 2015. 
#' This is used as the default value for the `rr.pop` argument when running 
#' BRCA-BCRAT. BRCA-BCRAT is run when either the `rr.bcrat` or `bcra.vars` 
#' argument is specified in a call to \code{\link{PPP}}.
#' 
"rr.ref"


#' Coefficients in BCRAT relative risk model
#' 
#' @format Race-specific coefficients from Gail et al. (1989) for White women, 
#' Gail et al. (2007) for Black women, Matsuno et al. (2011) for Asian women, 
#' and Banegas et al. (2017) for Hispanic women, saved as a data frame with 6 
#' rows and 7 columns:
#' * `beta1`: Coefficient for number of benign breast biopsies. 
#' * `beta2`: Coefficient for age at menarche. 
#' * `beta3`: Coefficient for age at first live birth. 
#' * `beta4`: Coefficient for number of female first-degree relatives with 
#' breast cancer. 
#' * `beta5`: Coefficient for interaction between number of benign breast 
#' biopsies and age>=50. 
#' * `beta6`: Coefficient for interaction between age at first live birth and 
#' number of female first-degree relatives with breast cancer. 
#' * `Race`: `"White"`, `"Black"`, `"Hispanic"`, `"FHispanic"`, `"Asian"`, or 
#' `"Other"`; `"Hispanic"` corresponds to U.S.-born Hispanic, while 
#' `"FHispanic"` corresponds to foreign-born Hispanic. 
#'
#' @references
#'
#' Gail MH, Brinton LA, Byar DP, Corle DK, Green SB, Schairer C, Mulvihill JJ. 
#' Projecting individualized probabilities of developing breast cancer for 
#' white females who are being examined annually. JNCI: Journal of the National 
#' Cancer Institute. 1989 Dec 20;81(24):1879-86.
#'
#' Gail MH, Costantino JP, Pee D, Bondy M, Newman L, Selvan M, Anderson GL, 
#' Malone KE, Marchbanks PA, McCaskill-Stevens W, Norman SA. Projecting 
#' individualized absolute invasive breast cancer risk in African American 
#' women. Journal of the National Cancer Institute. 2007 Dec 5;99(23):1782-92.
#'
#' Matsuno RK, Costantino JP, Ziegler RG, Anderson GL, Li H, Pee D, Gail MH. 
#' Projecting individualized absolute invasive breast cancer risk in Asian and 
#' Pacific Islander American women. Journal of the National Cancer Institute. 
#' 2011 Jun 22;103(12):951-61.
#'
#' Banegas MP, John EM, Slattery ML, Gomez SL, Yu M, LaCroix AZ, Pee D, 
#' Chlebowski RT, Hines LM, Thompson CA, Gail MH. Projecting individualized 
#' absolute invasive breast cancer risk in US Hispanic women. Journal of the 
#' National Cancer Institute. 2017 Feb 1;109(2):djw215.
#' 
"log.rr.bcrat"


#' Empirical age distributions used for age imputation
#'
#' A list containing empirical distributions of female and male ages at first 
#' birth; female and male birth ages; and sibling age differences. The values 
#' were obtained from the 2018 Natality Data from the National Vital Statistics 
#' System of the National Center for Health Statistics. Note that the maximum 
#' for the female ages is 50 because the data source coded ages above 50 as 50. 
#'
#' @details 
#' A list of 5 elements: 
#' * `female_first_birth`: Ages at first birth for females
#' * `male_first_birth`: Ages at first birth for males
#' * `female_birth`: Birth ages for females
#' * `male_birth`: Birth ages for males
#' * `sib_diff`: Age differences between pairs of siblings
#' 
#' @source \url{https://www.nber.org/research/data/vital-statistics-natality-birth-data}
"AgeDistribution"
