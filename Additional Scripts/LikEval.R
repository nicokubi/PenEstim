# Likelihood Evaluation
library(PPP)
library(numDeriv)
library(plyr) 
library(truncnorm)
library(PanelPRO)
library(tidyverse)
library(stringr)
library(survival)
library(MASS)
library(profvis)

LogLikelihood <- function(para, families) {
  # Set up the penetrance as a function of parameters alpha and beta
  alpha <- 1.75
  beta <- para
  gamma <- 0
  delta <- 60
  age <- seq(1, 94, 1)  # same as in the panelpro db
  
  
  # use same BRCA1 frequency (not estimated)
  BRCA1freq <- 0.9
  PanelPRODatabase$AlleleFrequency[["BRCA1_anyPV",1]] <- BRCA1freq
  PanelPRODatabase$AlleleFrequency[["BRCA1_anyPV",2]] <- BRCA1freq
  PanelPRODatabase$AlleleFrequency[["BRCA1_anyPV",3]] <- BRCA1freq
  
  
  # Weibull Penetrance Function for Females
  penetrance.mod.f <- dweibull(age-delta,alpha, beta) * (1-gamma)
  
  # For now focus on just one vector of penetrance estimates
  gene <- "BRCA1_hetero_anyPV"
  cancer <- "Breast"
  race <- "All_Races"
  female <- "Female"
  male <- "Male"
  type <- "Net"
  
  # Find the indices for the resp. attributes 
  dim_names <- attr(PanelPRODatabase$Penetrance, "dimnames")
  gene_index <- which(dim_names$Gene == gene)
  cancer_index <- which(dim_names$Cancer == cancer)
  race_index <- which(dim_names$Race == race)
  sex_index <- which(dim_names$Sex == female)
  type_index <- which(dim_names$PenetType == type)
  
  
  # overwrite the penetrance in the PanelPro Database 
  PanelPRODatabase$Penetrance[cancer_index, gene_index, race_index, sex_index, ,] <- penetrance.mod.f  
  
  # Male Penetrance
  sex_index <- which(dim_names$Sex == male)
  
  # Male penetrance funciton 
  penetrance.mod.m <- 0
  # overwrite the penetrance in the PanelPro Database 
  PanelPRODatabase$Penetrance[cancer_index, gene_index, race_index, sex_index, ,] <- penetrance.mod.m  
  
  
  # Storing the estimates
  log_likelihood <- 0
  
  for (i in 1:length(families)) {
    data <- families[[i]]  # Get the data for the current family
    
    # Access the posterior probabilities (not normalized) and estimates for BRCA1 gene in Breast cancer
    postprobs <- PPP(data, genes = c("BRCA1"), cancers = "Breast", database= PanelPRODatabase)$posterior.prob[[1]]  # check bug
    estimate <- postprobs[postprobs$genes=="BRCA1_hetero_anyPV","estimate"]
    
    if (is.nan(estimate) || estimate <=0) {
      # Handle NaN or zero probabilities by adding a small value and/or penalizing
      estimate <- 1e-28
      ll <- log(estimate)
      log_likelihood <- log_likelihood + ll
      cat("NaN or zero encountered:", alpha, beta, gamma, delta, ll, log_likelihood, estimate, "\n")
    } else {
      ll <- log(estimate)
      log_likelihood <- log_likelihood + ll
      cat(alpha, beta, gamma, delta, ll, log_likelihood, estimate, "\n")
    }
  }
  
  # return the negative log-likelihood so we can apply minimization 
  return(log_likelihood)
}


# Define the function that evaluates the likelihood for different values of beta
loglikelihood_evaluator <- function(paras, families) {
  loglikelihoods <- numeric(length(paras))
  
  for (i in seq_along(paras)) {
    # Call the NegLogLikelihood function for each beta value
    # Here, e are assuming that the families data are passed as a second argument
    # If it is stored in the global environment, you can remove the second argument
    loglikelihoods[i] <- LogLikelihood(c(paras[i]), families)
  }
  
  return(loglikelihoods)
}

# Define a range of beta values for evaluation
betas <- seq(12,18,1)

# Evaluate the likelihoods for different beta values
loglikelihoods <- loglikelihood_evaluator(betas, simFamilies_A_475_nocen)
results_df_SW3.5sb<- data.frame(beta = betas, logLikelihood = loglikelihoods)


# Basic R plotting
# Set up a 4x1 layout for the plots
# Create a PDF file with A4 dimensions
pdf("plots.pdf", width = 8.27, height = 11.69)



# Define a layout with 4 rows and 1 column

# Plot for alpha
# Set the outer margin to ensure square plots
par(oma = c(0, 0, 0, 0))
# Set the margin to ensure square plots
par(mar = c(5, 5, 4, 2) + 0.1)

# Create a 2x2 layout of plots with square dimensions
par(mfrow = c(2, 2))

# Families A
# Plot for alpha
plot(results_mle_A_famA$alpha, results_mle_A_famA$logLikelihood, type = "l", col = "blue", 
     xlab = "Alpha", ylab = "Log-likelihood", 
     main = "Shape Parameter (Alpha)")
abline(v = 2.5, col = "red", lty = 2)  

# Plot for beta
plot(results_mle_B_famA$beta, results_mle_B_famA$logLikelihood, type = "l", col = "blue", 
     xlab = "Beta", ylab = "Log-likelihood", 
     main = "Scale Parameter (Beta)")
abline(v = 14, col = "red", lty = 2)  

# Plot for gamma
plot(results_mle_G_famA$gamma, results_mle_G_famA$logLikelihood, type = "l", col = "blue", 
     xlab = "Gamma", ylab = "Log-likelihood", 
     main = "Asymptote Parameter (Gamma)")
abline(v = 0.8, col = "red", lty = 2)  

# Plot for delta
plot(results_mle_D_famA$delta, results_mle_D_famA$logLikelihood, type = "l", col = "blue", 
     xlab = "Delta", ylab = "Log-likelihood", 
     main = "Shift Parameter (Delta)")
abline(v = 60, col = "red", lty = 2)  


# Families B
# Plot for alpha
plot(results_mle_A_famB$alpha, results_mle_A_famB$logLikelihood, type = "l", col = "blue", 
     xlab = "Alpha", ylab = "Log-likelihood", 
     main = "Shape Parameter (Alpha)")
abline(v = 2.5, col = "red", lty = 2)  

# Plot for beta
plot(results_mle_B_famB$beta, results_mle_B_famB$logLikelihood, type = "l", col = "blue", 
     xlab = "Beta", ylab = "Log-likelihood", 
     main = "Scale Parameter (Beta)")
abline(v = 14, col = "red", lty = 2)  

# Plot for gamma
plot(results_mle_G_famB$gamma, results_mle_G_famB$logLikelihood, type = "l", col = "blue", 
     xlab = "Gamma", ylab = "Log-likelihood", 
     main = "Asymptote Parameter (Gamma)")
abline(v = 0.8, col = "red", lty = 2)  

# Plot for delta
plot(results_mle_D_famB$delta, results_mle_D_famB$logLikelihood, type = "l", col = "blue", 
     xlab = "Delta", ylab = "Log-likelihood", 
     main = "Shift Parameter (Delta)")
abline(v = 20, col = "red", lty = 2)  

#Families C
# Plot for alpha
plot(results_mle_A_famC$alpha, results_mle_A_famC$logLikelihood, type = "l", col = "blue", 
     xlab = "Alpha", ylab = "Log-likelihood", 
     main = "Shape Parameter (Alpha)")
abline(v = 2.5, col = "red", lty = 2)  

# Plot for beta
plot(results_mle_B_famC$beta, results_mle_B_famC$logLikelihood, type = "l", col = "blue", 
     xlab = "Beta", ylab = "Log-likelihood", 
     main = "Scale Parameter (Beta)")
abline(v = 50, col = "red", lty = 2)  

# Plot for gamma
plot(results_mle_G_famC$gamma, results_mle_G_famC$logLikelihood, type = "l", col = "blue", 
     xlab = "Gamma", ylab = "Log-likelihood", 
     main = "Asymptote Parameter (Gamma)")
abline(v = 0.9, col = "red", lty = 2)  

# Plot for delta
plot(results_mle_D_famC$delta, results_mle_D_famC$logLikelihood, type = "l", col = "blue", 
     xlab = "Delta", ylab = "Log-likelihood", 
     main = "Shift Parameter (Delta)")
abline(v = 20, col = "red", lty = 2)  




