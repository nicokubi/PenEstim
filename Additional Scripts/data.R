#Load required library 
library(dplyr)

# Read the original dataset
#data <- read.csv("/path/to/your/first/dataset.csv")
data <- bm.clean

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
        AgeOC = AgeOvary,
        AgePANC = AgePancreas,
        isAffBC = AffectedBreast,
        isAffOC = AffectedOvary,
        isAffPANC = AffectedPancreas,
        BRCA1 = BRCA1,
        MLH1 = MLH1,
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

# Step 4: Select the MLH 1 families 
MLH1Pedigree <- list_of_families_with_pedigree
# Add code to take only the families given that the proband has the mutation
MLH1Pedigree

# save fiel
#save(cohPedigree, file = "cohPedigree.RData")

# Filter families with affected probands for BRCA1
carrierProbandFamilies_cohPedigree_MLH1 <- Filter(function(fam) {
    any(fam$isProband == 1 & fam$MLH1 == 1)
}, cohPedigree)
str(carrierProbandFamilies_cohPedigree_MLH1)

print(carrierProbandFamilies_cohPedigree_MLH1[[2]])
length(carrierProbandFamilies_cohPedigree_MLH1)

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
selectfam(carrierProbandFamilies_cohPedigree, 100)
str(carrierProbandFamilies_cohPedigree_selected)
carrierProbandFamilies_cohPedigree_selected

age_vectors <- lapply(carrierProbandFamilies_cohPedigree_selected, function(sublist) {
    # Extracting the 'CurAge' vector from each sublist
    sublist[["CurAge"]]
})

# Step 2: Combine these vectors into one
all_ages <- unlist(age_vectors)

# Step 3: Analyze the distribution of ages

# Basic summary statistics
summary(all_ages)
str(carrierProbandFamilies_cohPedigree)


# Function to prepare ages
prepAges <- function(data) {
    for (i in seq_along(data)) {
        # For 'AgeBC', 'AgeOC', and 'AgePANC', set NA ages to 1 and corresponding 'isAff' to 0
        columns_to_check <- c("BC", "OC", "PANC")
        for (col_suffix in columns_to_check) {
            age_col <- paste0("Age", col_suffix)
            aff_col <- paste0("isAff", col_suffix)

            # Set age_col to 1 and corresponding 'isAff' to 0 where NA
            na_rows <- is.na(data[[i]][[age_col]])
            data[[i]][[age_col]][na_rows] <- 1
            data[[i]][[aff_col]][na_rows] <- 0
        }

        # Set 'CurAge' to the maximum cancer affection age where applicable
        cancer_ages <- data[[i]][c("AgeBC", "AgeOC", "AgePANC")]
        max_cancer_age <- apply(cancer_ages, 1, max, na.rm = TRUE)

        # If CurAge is NA, set it to max cancer age, else set it to the maximum of CurAge and max cancer age
        data[[i]]$CurAge <- ifelse(is.na(data[[i]]$CurAge), max_cancer_age, pmax(data[[i]]$CurAge, max_cancer_age))
    }
    return(data)
}

test_fam_1
test = prepAges(carrierProbandFamilies_cohPedigree_MLH1)
test
