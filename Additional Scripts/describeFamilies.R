describeFamilies <- function(fams, cancer, gene) {
  # Validate and translate cancer type
  if (!cancer %in% CANCER_NAME_MAP$long) {
    stop("Cancer type '", cancer, "' is not supported. Please choose from the supported list.")
  }
  
  # Get the abbreviation for the cancer type
  cancer_index <- which(CANCER_NAME_MAP$long == cancer)
  cancer_short <- CANCER_NAME_MAP$short[cancer_index]
  
  # Construct column names based on cancer type and gene input
  aff_col_name <- paste0("isAff", cancer_short)
  age_col_name <- paste0("Age", cancer_short)
  
  # Initialize counters and vectors for storing data
  affectedFamilies <- 0 
  affectedProbands <- 0
  multipleAffectedFamilies <- 0
  famSizes <- c()
  geneFamilies <- 0
  geneProbands <- 0 
  curAgesMale <- c()
  curAgesFemale <- c()
  cancerAges <- c()
  pbCancerAges <- c()
  pbCurAges <- c()
  pbCancerAgesFemale <- c()
  affectedFemaleProbands <- 0
  geneFemaleProbands <- 0 
  cancerAgesFemale <- c()
  cancerAgesMale <- c()
  total_individuals <- 0
  gene_individuals <- 0
  affected_individuals <- 0
  affected <- NA
  
  # Initialize counters for NA values
  na_curAge_male <- 0
  na_curAge_female <- 0
  na_cancerAge_male <- 0
  na_cancerAge_female <- 0
  
  # Process each family
  for (i in 1:length(fams)) {
    f <- fams[[i]]
    sizeOfFamily <- nrow(f)
    total_individuals <- total_individuals + sizeOfFamily
    famSizes <- c(famSizes, sizeOfFamily)
    
    curAge_male <- f %>% filter(Sex == 1) %>% pull(CurAge)
    curAge_female <- f %>% filter(Sex == 0) %>% pull(CurAge)
    
    # Count NA values for current age
    na_curAge_male <- na_curAge_male + sum(is.na(curAge_male))
    na_curAge_female <- na_curAge_female + sum(is.na(curAge_female))
    
    curAgesMale <- c(curAgesMale, curAge_male)
    curAgesFemale <- c(curAgesFemale, curAge_female)
    
    ff <- fams[[i]] %>% filter(!!sym(gene) == 1, !!sym(aff_col_name) == 1)
    fa <- fams[[i]] %>% filter(!!sym(aff_col_name) == 1)
    gene_individuals <- gene_individuals + nrow(ff)
    affected_individuals <- affected_individuals + nrow(fa)
    
    cancerAge_female <- fa %>% filter(Sex == 0) %>% pull(!!sym(age_col_name))
    cancerAge_male <- fa %>% filter(Sex == 1) %>% pull(!!sym(age_col_name))
    
    # Count NA values for cancer age
    na_cancerAge_male <- na_cancerAge_male + sum(is.na(cancerAge_male))
    na_cancerAge_female <- na_cancerAge_female + sum(is.na(cancerAge_female))
    
    cancerAgesFemale <- c(cancerAgesFemale, cancerAge_female)
    cancerAgesMale <- c(cancerAgesMale, cancerAge_male)
    
    if (nrow(fa) > 0) {
      affectedFamilies <- affectedFamilies + 1
    }
    
    if (nrow(fa) > 1) {
      multipleAffectedFamilies <- multipleAffectedFamilies + 1
    }
    
    geneFams <- fams[[i]] %>% filter(isProband == 0, !!sym(gene) == 1)
    if (nrow(geneFams) > 0) {
      geneFamilies <- geneFamilies + 1
    }
    
    pb <- fams[[i]] %>% filter(isProband == 1)
    pbfemale <- fams[[i]] %>% filter(isProband == 1, Sex == 0)
    
    if (nrow(pbfemale) > 0) {
      pbCancerAgesFemale <- c(pbCancerAgesFemale, pbfemale[[age_col_name]])
      if (any(pbfemale[[aff_col_name]] == 1)) {
        affectedFemaleProbands <- affectedFemaleProbands + 1
      }
      if (!is.na(pbfemale[[gene]])) {
        if (any(pbfemale[[gene]] == 1)) {
          geneFemaleProbands <- geneFemaleProbands + 1
        }
      }
    }
    
    pbCancerAges <- c(pbCancerAges, pb[[age_col_name]])
    pbCurAges <- c(pbCurAges, pb$CurAge)
    
    if (any(pb[[aff_col_name]] == 1)) {
      affectedProbands <- affectedProbands + 1
    }
    
    if (!is.na(pb[[gene]])) {
      if (any(pb[[gene]] == 1)) {
        geneProbands <- geneProbands + 1
      }
    }
  }
  
  print(paste0("Number of families: ", length(fams)))
  print(paste0("Total Individuals: ", total_individuals))
  print(paste0("Average family size: ", mean(famSizes)))
  print("Summary of family sizes")
  print(summary(famSizes))
  print("Summary of Current Age for all family members")
  print(summary(c(curAgesMale, curAgesFemale)))
  print(paste0("Number of individuals with NA Current Age (Male): ", na_curAge_male))
  print(paste0("Number of individuals with NA Current Age (Female): ", na_curAge_female))
  print("Summary of Current Age for probands")
  print(summary(pbCurAges))
  
  print(paste0("Number of affected individuals:", affected_individuals))
  print(paste0("Number of families with affected individuals (", cancer, "): ", affectedFamilies))
  print(paste0("Number of families with more than 1 affected individual (", cancer, "): ", multipleAffectedFamilies))
  print(paste0("Number of families with affected probands (", cancer, "): ", affectedProbands))
  print(paste0("Number of families with affected female probands (", cancer, "): ", affectedFemaleProbands))
  print(paste0("Number of families with relatives with PV (and not the proband) (", gene, "): ", geneFamilies))
  print(paste0("Number of probands with PV (", gene, "): ", geneProbands))
  print(paste0("Number of female probands with PV (", gene, "): ", geneFemaleProbands))
  print(paste0("Number of individuals with PV (", gene, "): ", gene_individuals))
  
  print("Summary of Cancer Age of (female) probands")
  print(summary(pbCancerAgesFemale))
  print("Summary of Cancer Age in affected (non-proband) individuals")
  print(summary(pbCancerAges))
  print(paste0("Number of individuals with NA Cancer Age (Male): ", na_cancerAge_male))
  print(paste0("Number of individuals with NA Cancer Age (Female): ", na_cancerAge_female))
  
  # Plot Cancer Age Distribution for Males and Females
  cancer_age_data <- data.frame(Age = c(cancerAgesMale, cancerAgesFemale),
                                Sex = c(rep("Male", length(cancerAgesMale)), rep("Female", length(cancerAgesFemale))))
  
  # Add NA counts as a separate category
  cancer_age_data <- rbind(cancer_age_data, data.frame(Age = rep(NA, na_cancerAge_male), Sex = rep("Male", na_cancerAge_male)))
  cancer_age_data <- rbind(cancer_age_data, data.frame(Age = rep(NA, na_cancerAge_female), Sex = rep("Female", na_cancerAge_female)))
  
  p1 <- ggplot(cancer_age_data, aes(x = Age, fill = Sex)) +
    geom_histogram(binwidth = 5, position = "dodge", color = "white", alpha = 0.8, na.rm = TRUE) +
    geom_text(aes(x = 85, y = 30, label = paste("NA Males:", na_cancerAge_male)), color = "blue") +
    geom_text(aes(x = 85, y = 28, label = paste("NA Females:", na_cancerAge_female)), color = "red") +
    geom_text(aes(x = 85, y = 32, label = paste("Total Affected:", affected_individuals)), color = "black") +
    theme_minimal() +
    labs(title = "Cancer Age Distribution", x = "Age", y = "Count") +
    xlim(0, 100) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  print(p1)
  
  # Plot Current Age Distribution for Males and Females
  cur_age_data <- data.frame(Age = c(curAgesMale, curAgesFemale),
                             Sex = c(rep("Male", length(curAgesMale)), rep("Female", length(curAgesFemale))))
  
  p2 <- ggplot(cur_age_data, aes(x = Age, fill = Sex)) +
    geom_histogram(binwidth = 5, position = "dodge", color = "white", alpha = 0.8) +
    geom_text(aes(x = 85, y = 120, label = paste("NA Males:", na_curAge_male)), color = "blue") +
    geom_text(aes(x = 85, y = 115, label = paste("NA Females:", na_curAge_female)), color = "red") +
    geom_text(aes(x = 85, y = 125, label = paste("Total Individuals:", total_individuals)), color = "black") +
    theme_minimal() +
    labs(title = "Current Age Distribution", x = "Age", y = "Count") +
    xlim(0, 100) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  print(p2)
}
