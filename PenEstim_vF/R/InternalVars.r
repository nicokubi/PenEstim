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
