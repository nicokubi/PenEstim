# Unexported internal parameters

#' Currently supported race types
RACE_TYPES <- c("All_Races", "AIAN", "Asian", "Black", 
                "White", "Hispanic", "WH", "WNH")

#' Currently supported ancestry types
ANCESTRY_TYPES <- c("AJ", "nonAJ", "Italian")

#' Currently supported sex types
SEX_TYPES <- c("Female", "Male")

#' Currently supported penetrance types
PENET_TYPES <- c("Net", "Crude")

#' Maximum age
MAXAGE <- 94


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
  "HEP", "CBC"), long = CANCER_TYPES)


#' Currently supported gene types
GENE_TYPES <- c(
  "ATM", "BARD1", "BRCA1", "BRCA2", 
  "BRIP1", "CDH1", "CDK4", "CDKN2A", "CHEK2", "EPCAM", 
  "MLH1", "MSH2", "MSH6", "MUTYH", "NBN", "PALB2",
  "PMS2", "PTEN", "RAD51C", "RAD51D", "STK11", "TP53"
)

#' Mapping of genes to supported variants
ALL_GENE_VARIANT_TYPES <- list(ATM = "ATM_hetero_anyPV", 
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
                               MUTYH = c("MUTYH_hetero_anyPV", 
                                         "MUTYH_homo_anyPV"), 
                               NBN = "NBN_hetero_657del5", 
                               PALB2 = "PALB2_hetero_anyPV",
                               PMS2 = "PMS2_hetero_anyPV", 
                               PTEN = "PTEN_hetero_anyPV", 
                               RAD51C = "RAD51C_hetero_anyPV", 
                               RAD51D = "RAD51D_hetero_anyPV", 
                               STK11 = "STK11_hetero_anyPV", 
                               TP53 = "TP53_hetero_anyPV")

#' Default variants to use
DEFAULT_VARIANTS <- c(ATM = "ATM_hetero_anyPV", 
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
                      MUTYH = "MUTYH_hetero_anyPV", 
                      NBN = "NBN_hetero_657del5", 
                      PALB2 = "PALB2_hetero_anyPV",
                      PMS2 = "PMS2_hetero_anyPV", 
                      PTEN = "PTEN_hetero_anyPV", 
                      RAD51C = "RAD51C_hetero_anyPV", 
                      RAD51D = "RAD51D_hetero_anyPV", 
                      STK11 = "STK11_hetero_anyPV", 
                      TP53 = "TP53_hetero_anyPV")

#' Genes with homozygous and heterozygous mutation variants
HOMOZYGOUS_VARIANT_TYPES <- list(
  "MUTYH_anyPV" = c("MUTYH_hetero_anyPV", "MUTYH_homo_anyPV"))

#' Genes with contralateral breast cancer associations
CBC_GENE_TYPES <- c("ATM", "BRCA1", "BRCA2", "CHEK2", "PALB2", "TP53")
CBC_GENE_VARIANT_TYPES <- list(ATM = "ATM_hetero_anyPV",
                               BRCA1 = "BRCA1_hetero_anyPV", 
                               BRCA2 = "BRCA2_hetero_anyPV",
                               CHEK2 = "CHEK2_hetero_1100delC",
                               PALB2 = "PALB2_hetero_anyPV",
                               TP53 = "TP53_hetero_anyPV")


#' Default race when missing from pedigree
UNKNOWN_RACE <- "All_Races"

#' Default ancestry when missing from pedigree
UNKNOWN_ANCESTRY <- "nonAJ"


#' Currently supported risk modifiers
RISKMODS <- list(BC = c("Mastectomy"), 
                 CBC = c("Mastectomy"),
                 ENDO = c("Hysterectomy"), 
                 OC = c("Oophorectomy"))

#' Currently supported marker testing
MARKER_TESTING <- list(
  BC = list(GENES = c("BRCA1_hetero_anyPV", "BRCA2_hetero_anyPV"), 
            MARKERS = c("ER", "CK5.6", "CK14", "PR", "HER2")),
  COL = list(GENES = c("MLH1_hetero_anyPV", "MSH2_hetero_anyPV", 
                       "MSH6_hetero_anyPV", "PMS2_hetero_anyPV"), 
             MARKERS = "MSI")
)


#' Model specification of genes and cancers for pre-defined models
MODELPARAMS <- list(
  BRCAPRO = list(
    GENES = c("BRCA1", "BRCA2"),
    CANCERS = c("Breast", "Ovarian", "Contralateral")
  ),
  BRCAPRO5 = list(
    GENES = c("BRCA1", "BRCA2", "ATM", "CHEK2", "PALB2"),
    CANCERS = c("Breast", "Ovarian", "Contralateral")
  ),
  BRCAPRO6 = list(
    GENES = c("BRCA1", "BRCA2", "MLH1", "MSH2", "MSH6", "CDKN2A"), 
    CANCERS = c("Breast", "Ovarian", "Contralateral")
  ),
  MMRPRO = list(
    GENES = c("MLH1", "MSH2", "MSH6"),
    CANCERS = c("Colorectal", "Endometrial")
  ),
  PanPRO11 = list(
    GENES = c(
      "BRCA1", "BRCA2", "ATM", "PALB2", "CHEK2",
      "EPCAM", "PMS2", "MLH1", "MSH2", "MSH6", "CDKN2A"
    ),
    CANCERS = c(
      "Brain", "Breast", "Colorectal", "Endometrial",
      "Gastric", "Kidney", "Melanoma", "Ovarian", "Pancreas", 
      "Prostate", "Small Intestine", "Contralateral"
    )
  ),
  PanPRO22 = list(
    GENES = c(
      "ATM", "BARD1", "BRCA1",
      "BRCA2", "BRIP1", "CDH1", "CDK4", "CDKN2A",
      "CHEK2", "EPCAM", "MLH1", "MSH2", "MSH6",
      "MUTYH", "NBN", "PALB2", "PMS2", "PTEN",
      "RAD51C", "RAD51D", "STK11", "TP53"
    ),
    CANCERS = c(
      "Brain", "Breast", "Colorectal", "Endometrial", "Gastric", "Kidney", 
      "Leukemia", "Melanoma", "Ovarian", "Osteosarcoma", "Pancreas", "Prostate",
      "Small Intestine", "Soft Tissue Sarcoma", "Thyroid", "Urinary Bladder", 
      "Hepatobiliary", "Contralateral"
    )
  )
)

# Global variables
utils::globalVariables(c("PanelPRODatabase", "log.rr.bcrat", "rr.ref",
                         "ByAge", "genes", "Cancer",
                         "estimate", "lower", "upper"))
