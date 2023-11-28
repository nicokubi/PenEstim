#Load required library 
library(dplyr)

# Read the original dataset
data <- read.csv("/path/to/your/first/dataset.csv")
data <- families

transformed_data <- data %>%
    select(
        ID = ID,
        Sex = Gender,
        MotherID = MotherID,
        FatherID = FatherID,
        isProband = isProband,
        CurAge = AgeDeath,
        isDead = Death,
        AgeBC = AgeBreast,
        isAffBC = AffectedBreast,
        BRCA1 = BRCA1,
        Twins = Twins,
        Ancestry = ethnic,
        race = race
    ) %>%
    mutate(
        riskmod = NA, # Add 'riskmod' with NA or empty string as needed
        interAge = NA # Add 'interAge' with NA or empty string as needed
    )

# Include 'PedigreeID' in the transformed data for grouping purpose
transformed_data_with_pedigree <- cbind(transformed_data, data["PedigreeID"])

# Step 2: Aggregation by Family/Pedigree
# Group by 'PedigreeID' to create family lists
families_with_pedigree <- transformed_data_with_pedigree %>%
    group_by(PedigreeID) %>%
    do(family_data = .[names(transformed_data)]) %>%
    ungroup()

# Step 3: Creating a List of Lists
# Convert the data frame of lists to a list of lists
list_of_families_with_pedigree <- split(families_with_pedigree$family_data, families_with_pedigree$PedigreeID)
list_of_families_with_pedigree <- lapply(list_of_families_with_pedigree, function(x) x %>% as.data.frame())
cohPedigree <- list_of_families_with_pedigree

# check structure of Pedigree file
str(cohPedigree)

# save fiel
save(cohPedigree, file = "cohPedigree.RData")

# Filter families with affected probands for BRCA1
carrierProbandFamilies_cohPedigree <- Filter(function(fam) {
    any(fam$isProband == 1 & fam$BRCA1 == 1 & fam$Sex == 0)
}, cohPedigree)
str(carrierProbandFamilies_cohPedigree)

print(carrierProbandFamilies_cohPedigree[[2]])

# Select a subset of the families 
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
selectfam(carrierProbandFamilies_cohPedigree, 200)
str(carrierProbandFamilies_cohPedigree_selected)


