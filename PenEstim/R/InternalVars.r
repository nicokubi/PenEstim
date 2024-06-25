#' Default Baseline Data
#'
#' This dataset contains age-specific cancer penetrance rates for both females and males.
#' As an example, the data is derived from the SEER program for breast cancer in females (for both males and females) 
#'
#' @format A data frame with 94 rows and 3 variables:
#' \describe{
#'   \item{Age}{Age in years (1 to 94)}
#'   \item{Female}{Penetrance rate for females}
#'   \item{Male}{Penetrance rate for males}
#' }
#'
baseline_data_default <- data.frame(
  Age = 1:94,
  Female = c(
    2.8e-07, 9e-08, 1e-08, 3e-08, 5e-08, 7e-08, 9e-08, 7e-08, 5e-08, 4e-08, 2e-08, 4e-08, 3.2e-07, 6.3e-07, 9.5e-07, 1.27e-06,
    1.94e-06, 4.73e-06, 7.88e-06, 1.102e-05, 1.417e-05, 1.916e-05, 3.517e-05, 5.303e-05, 7.088e-05, 8.873e-05, 0.00010961,
    1.486e-04, 1.906e-04, 0.00023261, 0.00027461, 0.00032101, 0.00039385, 0.00047109, 0.00054832, 0.00062554, 7.166e-04,
    0.00089081, 0.00107885, 0.00126684, 0.00145479, 0.00163827, 0.00179526, 0.00194779, 0.00210026, 0.00225264, 0.00239638,
    0.00248859, 0.00257215, 0.00265562, 0.00273901, 0.00281841, 0.00287439, 0.00292639, 0.00297831, 0.00303014, 0.00309486,
    0.00323735, 0.00339265, 0.00354775, 0.00370266, 0.00385852, 0.00402115, 0.0041847, 0.00434799, 4.511e-03, 0.00466219,
    0.00474384, 0.00481372, 0.00488337, 0.0049528, 0.00500362, 0.00494413, 0.00486620, 0.00478821, 0.00471016, 0.00462597,
    0.00450511, 0.00437818, 0.00425130, 0.00412450, 0.00399836, 0.00387577, 0.00375380, 0.00363187, 0.00351002, 0.00338186,
    0.00321541, 0.00304280, 0.00287050, 0.00269853, 0.00253244, 0.00239963, 0.00227262
  ),
  Male = c(
    2.8e-07, 9e-08, 1e-08, 3e-08, 5e-08, 7e-08, 9e-08, 7e-08, 5e-08, 4e-08, 2e-08, 4e-08, 3.2e-07, 6.3e-07, 9.5e-07, 1.27e-06,
    1.94e-06, 4.73e-06, 7.88e-06, 1.102e-05, 1.417e-05, 1.916e-05, 3.517e-05, 5.303e-05, 7.088e-05, 8.873e-05, 0.00010961,
    1.486e-04, 1.906e-04, 0.00023261, 0.00027461, 0.00032101, 0.00039385, 0.00047109, 0.00054832, 0.00062554, 7.166e-04,
    0.00089081, 0.00107885, 0.00126684, 0.00145479, 0.00163827, 0.00179526, 0.00194779, 0.00210026, 0.00225264, 0.00239638,
    0.00248859, 0.00257215, 0.00265562, 0.00273901, 0.00281841, 0.00287439, 0.00292639, 0.00297831, 0.00303014, 0.00309486,
    0.00323735, 0.00339265, 0.00354775, 0.00370266, 0.00385852, 0.00402115, 0.0041847, 0.00434799, 4.511e-03, 0.00466219,
    0.00474384, 0.00481372, 0.00488337, 0.0049528, 0.00500362, 0.00494413, 0.00486620, 0.00478821, 0.00471016, 0.00462597,
    0.00450511, 0.00437818, 0.00425130, 0.00412450, 0.00399836, 0.00387577, 0.00375380, 0.00363187, 0.00351002, 0.00338186,
    0.00321541, 0.00304280, 0.00287050, 0.00269853, 0.00253244, 0.00239963, 0.00227262
  )
)

# Unexported internal parameters
#' Currently supported cancer types
CANCER_TYPES <- c(
    "Brain", "Breast", "Colorectal", "Endometrial",
    "Gastric", "Kidney", "Leukemia", "Melanoma", "Ovarian",
    "Osteosarcoma", "Pancreas", "Prostate", "Small Intestine",
    "Soft Tissue Sarcoma", "Thyroid", "Urinary Bladder",
    "Hepatobiliary", "Contralateral"
)

#' Mapping of short and long cancer names
CANCER_NAME_MAP <- list(short = c(
    "BRA", "BC", "COL", "ENDO",
    "GAS", "KID", "LEUK", "MELA", "OC",
    "OST", "PANC", "PROS", "SI",
    "STS", "THY", "UB",
    "HEP", "CBC"
), long = CANCER_TYPES)

#' Currently supported gene types
GENE_TYPES <- c(
    "ATM", "BARD1", "BRCA1", "BRCA2",
    "BRIP1", "CDH1", "CDK4", "CDKN2A", "CHEK2", "EPCAM",
    "MLH1", "MSH2", "MSH6", "MUTYH", "NBN", "PALB2",
    "PMS2", "PTEN", "RAD51C", "RAD51D", "STK11", "TP53"
)

#' Mapping of genes to supported variants
ALL_GENE_VARIANT_TYPES <- list(
    ATM = "ATM_hetero_anyPV",
    BARD1 = "BARD1_hetero_anyPV",
    BRCA1 = "BRCA1_hetero_anyPV",
    BRCA2 = "BRCA2_hetero_anyPV",
    BRIP1 = "BRIP1_hetero_anyPV",
    CDH1 = "CDH1_hetero_anyPV",
    CDK4 = "CDK4_hetero_anyPV",
    CDKN2A = "CDKN2A[P16]_hetero_anyPV",
    CHEK2 = "CHEK2_hetero_1100delC",
    EPCAM = "EPCAM_hetero_anyPV",
    MLH1 = "MLH1_hetero_anyPV",
    MSH2 = "MSH2_hetero_anyPV",
    MSH6 = "MSH6_hetero_anyPV",
    MUTYH = c(
        "MUTYH_hetero_anyPV",
        "MUTYH_homo_anyPV"
    ),
    NBN = "NBN_hetero_657del5",
    PALB2 = "PALB2_hetero_anyPV",
    PMS2 = "PMS2_hetero_anyPV",
    PTEN = "PTEN_hetero_anyPV",
    RAD51C = "RAD51C_hetero_anyPV",
    RAD51D = "RAD51D_hetero_anyPV",
    STK11 = "STK11_hetero_anyPV",
    TP53 = "TP53_hetero_anyPV"
)
