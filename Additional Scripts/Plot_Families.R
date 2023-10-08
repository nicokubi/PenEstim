library(tidyverse)
library(stringr)
library(MASS)
library(survival)

# modify the penetrance for simulation
alpha <- 1
beta <- 50
gamma <- 0
delta <- 25
age <- seq(1,94,1) # same as in the panelpro db


prob_plot_penetrance_curves <- function(fams) {
  
  # Initialize the data frames before calling the function
  allcarriers <- data.frame()
  allnoncarriers <- data.frame()
  
  # Create an empty plot
  plot(0, 0, type = "n", xlim = c(0, 100), ylim = c(0, 1), xlab = "Age", 
       ylab = "Penetrance",
       main = "Estimated Penetrance Curves (for females)")
  
  # Initialize colors and line types
  colors <- c("red", "blue")  
  line_types <- c(1, 2)  
  
  # Iterate through each family and aggregate the data
  for (i in 1:length(fams)) {
    f <- fams[[i]]
    
    # Filter only probands
    f_proband <- f
    carriers <- f_proband %>% filter(BRCA1 == 1) %>% filter(Sex == 0)
    noncarriers <- f_proband %>% filter(BRCA1 == 0) %>% filter(Sex == 0)
    
    # Append the data to the combined data frames
    allcarriers <- rbind(allcarriers, carriers)
    allnoncarriers <- rbind(allnoncarriers, noncarriers)
  }
  
  # Convert the cancer status variable to survival objects for each genotype
  surv_obj_1 <- Surv(time = allcarriers$AgeBC, event = allcarriers$isAffBC)
  #surv_obj_0 <- Surv(time = allnoncarriers$AgeBC, 
                     #event = allnoncarriers$isAffBC)
  
  # Estimate the Kaplan-Meier survival curves for each genotype
  km_fit_1 <- survfit(surv_obj_1 ~ 1, data = allcarriers)
  #km_fit_0 <- survfit(surv_obj_0 ~ 1, data = allnoncarriers)
  
  # Calculate the penetrance functions for carriers and non-carriers
  penetrance_1 <- 1 - km_fit_1$surv
  #penetrance_0 <- 1 - km_fit_0$surv
  
  # Plot the aggregated penetrance curves
  lines(km_fit_1$time, penetrance_1, type = "l", col = colors[1], 
        lty = line_types[1])
  #lines(km_fit_0$time, penetrance_0, type = "l", col = colors[2], 
   #     lty = line_types[2])
  
  # Return the aggregated data frames
  return(list(allcarriers = allcarriers, allnoncarriers = allnoncarriers))
}

# Assuming 'families' is the list of family data frames
# Plot the empirical penetrance curve
par(mfrow = c(1, 1))
result <- prob_plot_penetrance_curves(sim)

# Assuming 'age', 'delta', 'alpha', 'beta', and 'gamma' are defined
lines(age, pweibull(age-60,2.5,14) * (1), type = "l", 
      lwd = 1, col = "black")
# Add legend
legend("topleft", legend = c("Empirical Distribution - Carrier", 
                             "True Carrier Weibull Distribution"), 
       lty = c(1, 1), lwd = c(1, 1), col = c("red", "black"))


result <- prob_plot_penetrance_curves(simFamilies_B_1000_nocen)

# Assuming 'age', 'delta', 'alpha', 'beta', and 'gamma' are defined
lines(age, pweibull(age-20,2.5,14) * (0.8), type = "l", 
      lwd = 1, col = "black")
# Add legend
legend("topleft", legend = c("Empirical Distribution - Carrier", 
                             "True Carrier Weibull Distribution"), 
       lty = c(1, 1), lwd = c(1, 1), col = c("red", "black"))

result <- prob_plot_penetrance_curves(simFamilies_C_1000)

# Assuming 'age', 'delta', 'alpha', 'beta', and 'gamma' are defined
lines(age, pweibull(age-20,2.5,50) * (0.9), type = "l", 
      lwd = 1, col = "black")
# Add legend
legend("topleft", legend = c("Empirical Distribution - Carrier", 
                             "True Carrier Weibull Distribution"), 
       lty = c(1, 1), lwd = c(1, 1), col = c("red", "black"))
     
