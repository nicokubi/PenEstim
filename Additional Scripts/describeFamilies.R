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
  
  # Process each family
  for (i in 1:length(fams)) {
    f <- fams[[i]]
    sizeOfFamily <- nrow(f)
    total_individuals <- total_individuals + sizeOfFamily
    famSizes <- c(famSizes, sizeOfFamily)
    curAgesMale <- c(curAgesMale, f %>% filter(Sex == 1) %>% pull(CurAge))
    curAgesFemale <- c(curAgesFemale, f %>% filter(Sex == 0) %>% pull(CurAge))
    
    ff <- fams[[i]] %>% filter(!!sym(gene) == 1, !!sym(aff_col_name) == 1)
    fa <- fams[[i]] %>% filter(!!sym(aff_col_name) == 1)
    gene_individuals <- gene_individuals + nrow(ff)
    affected_individuals <- affected_individuals + nrow(fa)
    
    cancerAgesFemale <- c(cancerAgesFemale, ff %>% filter(Sex == 0) %>% pull(!!sym(age_col_name)))
    cancerAgesMale <- c(cancerAgesMale, ff %>% filter(Sex == 1) %>% pull(!!sym(age_col_name)))
    
    if (nrow(fa) > 0) {
      affectedFamilies <- affectedFamilies + 1
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
  print("Summary of Current Age for probands")
  print(summary(pbCurAges))
  
  print(paste0("Number of affected individuals:", affected_individuals))
  print(paste0("Number of families with affected individuals (", cancer, "): ", affectedFamilies))
  print(paste0("Number of families with affected probands (", cancer, "): ", affectedProbands))
  print(paste0("Number of families with affected female probands (", cancer, "): ", affectedFemaleProbands))
  print(paste0("Number of families with relatives with PV (and not the proband) (", gene, "): ", geneFamilies))
  print(paste0("Number of probands with PV (", gene, "): ", geneProbands))
  print(paste0("Number of female probands with PV (", gene, "): ", geneFemaleProbands))
  print(paste0("Number of individuals with PV (", gene, "): ", gene_individuals))
  
  print("Summary of Cancer Age of (female) probands")
  print(summary(pbCancerAgesFemale))
  print("Summary of Cancer Age in affected (non-proband) individuals")
  print(summary(cancerAges))
  
  # Plot Cancer Age Distribution for Males and Females
  cancer_age_data <- data.frame(Age = c(cancerAgesMale, cancerAgesFemale),
                                Sex = c(rep("Male", length(cancerAgesMale)), rep("Female", length(cancerAgesFemale))))
  
  p1 <- ggplot(cancer_age_data, aes(x = Age, fill = Sex)) +
    geom_histogram(binwidth = 5, position = "dodge", color = "white", alpha = 0.8) +
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
