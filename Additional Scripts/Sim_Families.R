# Simulation
library(plyr) #need to load plyr before dplyr
library(truncnorm)
library(PanelPRO)
library(tidyverse)
library(stringr)
library(PedUtils)
library(survival)
library(MASS)
library(profvis)

# set seed 
set.seed(2024)

# Generate Families ----

families_sim1= list()
probandIDS = c()
probandBRCA1Status = c()
probandAffectionStatus = c()


# Set number of families to be generates
numberFamilies <- 2000
# age
age <- seq(1, 94, 1)  # same as in the panelpro db

# use same BRCA1 frequency (not estimated)
BRCA1freq <- 0.1
PanelPRODatabase$AlleleFrequency[["BRCA1_anyPV",1]] <- BRCA1freq
PanelPRODatabase$AlleleFrequency[["BRCA1_anyPV",2]] <- BRCA1freq
PanelPRODatabase$AlleleFrequency[["BRCA1_anyPV",3]] <- BRCA1freq

# Given mean and variance
alpha<- 2.5 # Replace with the value you have
beta  <- 50 # Replace with the value you have
gamma <- 0.9
delta <-  20

median_true <- delta + beta * (-log((gamma-0.5)/gamma))^(1/alpha)
quartile_true= delta + beta * (-log((gamma-0.25)/gamma))^(1/alpha)
print(median_true)
print(quartile_true)

# Now use alpha and beta in your simulation
penetrance.mod <- dweibull(age - delta, alpha, beta) * gamma


# For now focus on just one vector of penetrance estimates
gene <- "BRCA1_hetero_anyPV"
cancer <- "Breast"
race <- "All_Races"
sex <- "Female"
type <- "Net"

# Find the indices for the resp. attributes 

dim_names <- attr(PanelPRODatabase$Penetrance, "dimnames")
gene_index <- which(dim_names$Gene == gene)
cancer_index <- which(dim_names$Cancer == cancer)
race_index <- which(dim_names$Race == race)
sex_index <- which(dim_names$Sex == sex)
type_index <- which(dim_names$PenetType == type)

# Overwrite the penetetrane for all races, geneders and races
PanelPRODatabase$Penetrance[cancer_index, gene_index, race_index, sex_index , ,]<- penetrance.mod

for(i in 1:numberFamilies){
  # Cancers
  cancers = "Breast"
  # Genes
  genes = "BRCA1"
  #family members
  # Paternal aunts, paternal uncles
  nSibsPatern =floor(rtruncnorm(n=2, mean=0, 1))
  # Maternal aunts, maternal uncles
  nSibsMatern = floor(rtruncnorm(n=2, mean=0, 1))
  # Sisters and brothers
  nSibs = floor(rtruncnorm(n=2, mean=1, 1))
  # We make the assumption that the number of sons and daughters for the
  # proband and all siblings, is the same. Nieces and nephews of the proband
  # are not sampled separately1
  nGrandchild = floor(rtruncnorm(n=2, mean=0, 1))
  nChild = floor(rtruncnorm(n=2, mean=1, 1))
  
  # Simulate family using `PedUtils` code
  fam = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nChild,
                      PanelPRODatabase, genes, cancers,
                      includeGeno = TRUE, includeBiomarkers = FALSE, 
                      censoring = TRUE)
  
  famDF = as.data.frame(fam)
  proband = famDF %>% filter(isProband==1)
  probandIDS = c(probandIDS, proband$ID)
  probandBRCA1Status = c(probandBRCA1Status, proband$BRCA1)
  probandAffectionStatus = c(probandAffectionStatus, proband$isAffBC)
  families_sim1[[i]] = famDF
  
}

save(families_sim1, file = "families_sim1.RData")
# Filter families with affected probands
carrierProbandfamilies_sim1 <- Filter(function(fam) 
  any(fam$isProband == 1 & fam$BRCA1 == 1 & fam$Sex == 0), families_sim1)


# hide genotype information for everyone but the proband
simFamiliesGeno <- function(fams) {
  
  simFamilies <- list()
  # Assuming you have the original list of families stored as `original_families`
  # and a vector of proband IDs called `proband_ids`
  for (i in 1:length(fams)) {
    family <- fams[[i]]
    proband <-  family %>% filter(isProband==1)
    family$BRCA1 <- ifelse(family$ID == proband$ID, family$BRCA1, NA)
    simFamilies[[i]] <- family
  }
  return(simFamilies)
  
}

# Subset the family list
simfamilies_sim1 <- simFamiliesGeno(carrierProbandfamilies_sim1)

save(simfamilies_sim1, file = "simfamilies_sim1.Rdata")

# Function to select families 

selectfam <- function(input_families, n) {
  selectedIndices <- sample(length(input_families), n)
  selectedFamilies <- list()
  
  for (i in selectedIndices) {
    selectedFamilies[[length(selectedFamilies) + 1]] <- input_families[[i]]
  }
  
  new_name <- paste(deparse(substitute(input_families)), "_selected", sep = "")
  assign(new_name, selectedFamilies, envir = .GlobalEnv)
  
  save(list = new_name, file = paste0(new_name, ".RData"))
}
selectfam(simfamilies_sim1,60)

simFamilies_C_1000_nocen_selected60 <- simFamilies_C_1000_nocen_selected300

save(simfamilies_sim1_selected, file = "simfamilies_sim1_selected.Rdata")

