#' Calculate BCRAT relative risks
#'
#' @param df A data frame of BCRAT covariates, where each row corresponds to 
#' a proband in the pedigree. The expected columns are: 
#' * `ID`: Numeric ID for each proband. Should not contain duplicated 
#' entries. 
#' * `Age`: Initial age (should be >=20 and <90). 
#' * `NBiops`: Number of benign breast biopsies (`99` if unknown). 
#' * `Hyp`: A logical value indicating whether the individual has had a breast 
#' biopsy with atypical hyperplasia (`0` if no, `1` if yes, `99` if unknown). 
#' * `AgeMen`: Age at menarche (`99` if unknown)
#' * `AgeFLB`: Age at first live birth (`98` if nulliparous, `99` if unknown). 
#' * `NumRel`: The number of female first-degree relatives with breast cancer 
#' (`99` if unknown). 
#' * `Race`: A character string indicating race (`"White"`, `"Black"`, 
#' `"Hispanic"`, `"FHispanic"`, `"Asian"`, or `"Other"`; `"Hispanic"` 
#' corresponds to U.S.-born Hispanic while `"FHispanic"` corresponds to 
#' foreign-born Hispanic). 
#' 
#' @family bcrat
#' @seealso \code{\link{check.inputs.bcrat}}, \code{\link{recode.inputs.bcrat}}
#' @return A data frame of BCRAT relative risk estimates with three columns:
#' * `ID`: Numeric ID for each individual. 
#' * `rr1`: Individual's BCRAT relative risk estimate for ages < 50. 
#' * `rr2`: Individual's BCRAT relative risk estimate for ages >= 50. 
#' 
#' @details 
#' The model for White women is described in Gail et al. (1989) and Costantino 
#' et al. (1999); the model for African-American women is described in Gail et 
#' al. (2007); the model for Asian women is described in Matsuno et al. (2011); 
#' and the model for Hispanic women is described in Banegas et al. (2017).
#'
#' @references
#' Gail MH, Brinton LA, Byar DP, Corle DK, Green SB, Schairer C, Mulvihill JJ. 
#' Projecting individualized probabilities of developing breast cancer for 
#' white females who are being examined annually. JNCI: Journal of the National 
#' Cancer Institute. 1989 Dec 20;81(24):1879-86.
#'
#' Costantino JP, Gail MH, Pee D, Anderson S, Redmond CK, Benichou J, Wieand 
#' HS. Validation studies for models projecting the risk of invasive and total 
#' breast cancer incidence. Journal of the National Cancer Institute. 1999 Sep 
#' 15;91(18):1541-8.
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
#' @examples
#' bcrat.input <- data.frame(
#'   ID = 1:6, Age = 45:50, NBiops = c(0, 0, 0, 99, 1, 3),
#'   Hyp = c(99, 99, 99, 99, 0, 1), AgeMen = c(11, 12, 13, 14, 99, 15),
#'   AgeFLB = c(20, 22, 98, 99, 24, 30), NumRel = c(0, 1, 2, 0, 0, 0),
#'   Race = c("White", "White", "Black", "Hispanic", "Asian", "Other")
#' )
#' PPP:::calc.rr.bcrat(bcrat.input)
calc.rr.bcrat <- function(df) {
  
  # Check data frame inputs
  df <- check.inputs.bcrat(df)

  # Recode covariates to get categorical variables and Hyp multiplicative 
  # factor for logistic regression model
  df <- recode.inputs.bcrat(df)

  # Set up race-specific coefficients for each individual
  df[, paste0("beta", 1:6)] <- NA
  df[which(df$Race %in% "White"), paste0("beta", 1:6)] <- 
    rbind(log.rr.bcrat[which(log.rr.bcrat$Race %in% "White"), 
                       paste0("beta", 1:6)][rep(1, length(which(df$Race %in% "White"))), ])
  df[which(df$Race %in% "Black"), paste0("beta", 1:6)] <- 
    rbind(log.rr.bcrat[which(log.rr.bcrat$Race %in% "Black"), 
                       paste0("beta", 1:6)][rep(1, length(which(df$Race %in% "Black"))), ])
  df[which(df$Race %in% "Hispanic"), paste0("beta", 1:6)] <- 
    rbind(log.rr.bcrat[which(log.rr.bcrat$Race %in% "Hispanic"), 
                       paste0("beta", 1:6)][rep(1, length(which(df$Race %in% "Hispanic"))), ])
  df[which(df$Race %in% "FHispanic"), paste0("beta", 1:6)] <- 
    rbind(log.rr.bcrat[which(log.rr.bcrat$Race %in% "FHispanic"), 
                       paste0("beta", 1:6)][rep(1, length(which(df$Race %in% "FHispanic"))), ])
  df[which(df$Race %in% "Other"), paste0("beta", 1:6)] <- 
    rbind(log.rr.bcrat[which(log.rr.bcrat$Race %in% "Other"), 
                       paste0("beta", 1:6)][rep(1, length(which(df$Race %in% "Other"))), ])
  df[which(df$Race %in% "Asian"), paste0("beta", 1:6)] <- 
    rbind(log.rr.bcrat[which(log.rr.bcrat$Race %in% "Asian"), 
                       paste0("beta", 1:6)][rep(1, length(which(df$Race %in% "Asian"))), ])

  # Calculate relative risks
  # Age < 50
  df$rr1 <- exp(df$beta1 * df$NBiopsCat + df$beta2 * df$AgeMenCat + 
                  df$beta3 * df$AgeFLBCat + df$beta4 * df$NumRelCat + 
                  df$beta6 * df$AgeFLBCat * df$NumRelCat + log(df$HypMF))
  # Age >- 50
  df$rr2 <- df$rr1 * exp(df$beta5 * df$NBiopsCat)
  
  # Set relative risk estimates with errors (except "Age error") to NA
  df[which(!df$error %in% c(NA, "Age error")), c("rr1", "rr2")] <- NA

  return(df[, c("ID", "rr1", "rr2")])
}


#' Recode BCRAT covariates
#'
#' @param df A dataframe of BCRAT covariates. (See the description of the `df` 
#' parameter in the \code{\link{calc.rr.bcrat}} documentation.)
#' 
#' @return Same data frame as input, but with five additional columns: 
#' * `NBiopsCat`: Categorical variable for number of benign breast biopsies 
#' (`NBiops`). 
#' * `AgeMenCat`: Categorical variable for age at menarche (`AgeMen`). 
#' * `AgeFLBCat`Categorical variable for age at first live birth (`AgeFLB`). 
#' * `NumRelCat`: Categorical variable for the number of female first-degree 
#' relatives with breast cancer (`NumRel`). 
#' * `HypMF`: Multiplicative factor for breast biopsies with atypical 
#' hyperplasia (`Hyp`). 
#' @seealso \code{\link{calc.rr.bcrat}}
#' @family bcrat
recode.inputs.bcrat <- function(df) {

  # Initialize new categorical variable columns
  df[, paste0(c("NBiops", "AgeMen", "AgeFLB", "NumRel"), "Cat")] <- NA

  # Categorical variable for number of benign breast biopsies (0, 1, 2)
  # 0, 99, NA -> category 0
  df$NBiopsCat[which(df$NBiops %in% c(0, 99, NA))] <- 0
  # 1 -> category 1
  df$NBiopsCat[which(df$NBiops %in% 1)] <- 1
  # >=2, not 99 -> category 2
  df$NBiopsCat[which(df$NBiops > 1 & df$NBiops < 99)] <- 2
  # US- and foreign-born Hispanic models group category 2 with category 1
  df$NBiopsCat[which(df$NBiops > 1 & df$NBiops < 99 & 
                       df$Race %in% c("Hispanic", "FHispanic"))] <- 1

  # Categorical variable for age at menarche (0, 1, 2)
  # >14, 99, NA -> category 0
  df$AgeMenCat[which(df$AgeMen %in% c(99, NA) | df$AgeMen >= 14)] <- 0
  # [12, 14) -> category 1
  df$AgeMenCat[which(df$AgeMen >= 12 & df$AgeMen < 14)] <- 1
  # <12 -> category 2
  df$AgeMenCat[which(df$AgeMen < 12)] <- 2
  # African-American model groups codes <12 -> category 1 (category 2 grouped 
  # with category 1)
  df$AgeMenCat[which(df$AgeMen < 12 & df$Race %in% "Black")] <- 1
  # US-born Hispanic-American model recodes all values to the same category
  df$AgeMenCat[which(df$Race %in% "Hispanic")] <- 0
  
  # Categorical variable for age at first live birth (0, 1, 2, 3)
  # <20, 99, NA -> category 0
  df$AgeFLBCat[which(df$AgeFLB %in% c(99, NA) | df$AgeFLB < 20)] <- 0
  # [20, 25) -> category 1
  df$AgeFLBCat[which(df$AgeFLB >= 20 & df$AgeFLB < 25)] <- 1
  # [25, 30), 98 -> category 2
  df$AgeFLBCat[which((df$AgeFLB >= 25 & df$AgeFLB < 30) | df$AgeFLB %in% 98)] <- 2
  # >30, not 98 or 99 -> category 3
  df$AgeFLBCat[which(df$AgeFLB >= 30 & df$AgeFLB < 98)] <- 3
  # African-American model recodes all values to the same category
  df$AgeFLBCat[which(df$Race %in% "Black")] <- 0
  # US- and foreign-born Hispanic-American models code [20, 30) -> category 1
  df$AgeFLBCat[which(df$Race %in% c("Hispanic", "FHispanic") & 
                       df$AgeFLB >= 20 & 
                       df$AgeFLB < 30)] <- 1
  # US- and foreign-born Hispanic-American models code 
  # >30, not 98 or 99 -> category 2
  df$AgeFLBCat[which(df$Race %in% c("Hispanic", "FHispanic") & 
                       (df$AgeFLB >= 30 | df$AgeFLB %in% 98) & 
                       df$AgeFLB < 99)] <- 2
  
  # Categorical variable for number of first-degree relatives with breast 
  # cancer (0, 1, 2)
  # 0, 99, NA -> category 0
  df$NumRelCat[which(df$NumRel %in% c(99, NA, 0))] <- 0
  # 1 -> category 1
  df$NumRelCat[which(df$NumRel %in% 1)] <- 1
  # >=2, not 99 -> category 2
  df$NumRelCat[which(df$NumRel >= 2 & df$NumRel < 99)] <- 2
  # Asian-American and US- and foreign-born Hispanic-American models code 
  # >=2, not 99 -> category 1 (category 2 grouped with category 1)
  df$NumRelCat[which(df$NumRel >= 2 & df$NumRel < 99 & 
                       df$Race %in% c("Hispanic", "FHispanic", "Asian"))] <- 1

  # Multiplicative factor for breast biopsies with atypical hyperplasia
  df$HypMF <- NA
  # No benign breast biosopies or missing/unknown atypical hyperplasia status
  df$HypMF[which(df$NBiopsCat == 0 | df$Hyp %in% c(NA, 99))] <- 1
  # Nonzero number of biosopies and no biopsies with atypical hyperplasia 
  df$HypMF[which(df$NBiopsCat > 0 & df$Hyp == 0)] <- 0.93
  # Nonzero number of biosopies and biopsy with atypical hyperplasia 
  df$HypMF[which(df$NBiopsCat > 0 & df$Hyp == 1)] <- 1.82

  return(df)
}


#' Check validity of BCRAT covariates
#'
#' @param df A dataframe of BCRAT covariates. (See the description of the `df` 
#' parameter in the \code{\link{calc.rr.bcrat}} documentation.)
#'
#' @return Same data frame as input, but with an additional column named 
#' `error` containing error messages for rows with invalid covariate values. 
#' `NA` indicates no error. 
#' @seealso \code{\link{calc.rr.bcrat}}
#' @family bcrat
check.inputs.bcrat <- function(df) {
  
  # Initialize error column
  df$error <- NA

  # Error if age is out of bounds
  ind.age <- which(df$Age < 20 | df$Age >= 90)
  df$error[ind.age] <- "Age error"

  # Error if an unsupported race is specified
  ind.race <- which(!df$Race %in% c("White", "Black", "Hispanic", 
                                    "FHispanic", "Other", "Asian"))
  df$error[ind.race] <- "Race error"

  # Error if biopsy atypical hyperplasia status is invalid or incompatible with 
  # number of benign breast biopsies
  # If number of biopsies is 0 or unknown (99), should have either unknown (99) 
  # or no biopsies with atypical (0) hyperplasia
  # Otherwise, atypical hyperplasia status should be 0, 1, or 99
  ind.biop <- which((df$NBiops %in% c(0, 99) & !df$Hyp %in% c(0, 99, NA)) | 
                      (df$NBiops > 0 & !df$Hyp %in% c(99, NA, 0, 1)))
  df$error[ind.biop] <- "Hyp error"

  # Error if age at menarche is invalid. Ignore errors for US-born Hispanic 
  # women (not currently in model).
  ind.agemen <- which((df$AgeMen < 0 | (df$AgeMen > df$Age & df$AgeMen < 99)) & 
                        !df$Race %in% "Hispanic")
  df$error[ind.agemen] <- "AgeMen error"

  # Error if age at first live birth is invalid or occurs before age of 
  # menarche. Ignore errors for black women (not currently in model).
  ind.ageflb <- which((df$AgeFLB < 0 | (df$AgeFLB > df$Age & df$AgeFLB < 98) | 
                         (df$AgeFLB < df$AgeMen & df$AgeMen < 99)) & 
                        !df$Race %in% "Black")
  df$error[ind.ageflb] <- "AgeFLB error"

  return(df)
}
