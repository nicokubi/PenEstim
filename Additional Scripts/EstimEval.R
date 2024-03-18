#  Plotting of the estimated penetrance and the data-generating curve for simulation studies

plot_penetrance_sim <- function(data, prob, max_age, data_gen_alpha, data_gen_beta, data_gen_threshold, data_gen_asymptote) {
    if (prob <= 0 || prob >= 1) {
        stop("prob must be between 0 and 1")
    }

    params <- calculate_weibull_parameters(
        data$median_results,
        data$first_quartile_results,
        data$threshold_results,
        data$asymptote_results
    )

    alphas <- params$alpha
    betas <- params$beta
    asymptotes <- data$asymptote_results
    thresholds <- data$threshold_results

    x_values <- seq(0, max_age, length.out = max_age + 1)
    distributions <- vector("list", length(alphas))

    for (i in seq_along(alphas)) {
        if (validate_weibull_parameters(
            data$first_quartile_results[i],
            data$median_results[i],
            data$threshold_results[i],
            data$asymptote_results[i]
        )) {
            distributions[[i]] <- pweibull(x_values - thresholds[i],
                shape = alphas[i], scale = betas[i]
            ) * asymptotes[i]
        } else {
            distributions[[i]] <- rep(NA, length(x_values))
        }
    }

    distributions_matrix <- do.call(cbind, distributions)
    ci_lower <- apply(distributions_matrix, 1, function(x) quantile(x, probs = (1 - prob) / 2, na.rm = TRUE))
    ci_upper <- apply(distributions_matrix, 1, function(x) quantile(x, probs = 1 - (1 - prob) / 2, na.rm = TRUE))
    mean_density <- rowMeans(distributions_matrix, na.rm = TRUE)

    # Data-generating curve using the given parameters
    data_generating_curve <- pweibull(x_values - data_gen_threshold,
        shape = data_gen_alpha, scale = data_gen_beta
    ) * data_gen_asymptote

    par(mfrow = c(1, 1))

    plot(x_values, mean_density,
        type = "l", col = "blue", ylim = c(min(ci_lower), max(ci_upper)), xlim = c(1, max_age),
        xlab = "Age", ylab = "Estimated Penetrance", main = "Penetrance Curve with Credible Interval"
    )
    lines(x_values, ci_lower, col = "red", lty = 2)
    lines(x_values, ci_upper, col = "red", lty = 2)
    lines(x_values, data_generating_curve, col = "black", lty = 1)
    polygon(c(x_values, rev(x_values)), c(ci_lower, rev(ci_upper)), col = rgb(1, 0, 0, 0.1), border = NA)

    legend("topleft",
        legend = c("Estimated Penetrance", "Credible Interval", "Data-Generating Penetrance"),
        col = c("blue", "red", "black"), lty = c(1, 2, 1), cex = 0.8
    )
}

# Calculating the required penetrance data 

calculate_penetrance_data <- function(data, prob, max_age, data_gen_alpha, data_gen_beta, data_gen_threshold, data_gen_asymptote) {
    if (prob <= 0 || prob >= 1) {
        stop("prob must be between 0 and 1")
    }

    params <- calculate_weibull_parameters(
        data$median_results,
        data$first_quartile_results,
        data$threshold_results,
        data$asymptote_results
    )

    alphas <- params$alpha
    betas <- params$beta
    thresholds <- data$threshold_results
    asymptotes <- data$asymptote_results

    x_values <- seq(0, max_age, length.out = max_age + 1)
    distributions <- vector("list", length(alphas))

    for (i in seq_along(alphas)) {
        distributions[[i]] <- pweibull(x_values - thresholds[i], shape = alphas[i], scale = betas[i]) * asymptotes[i]
    }

    distributions_matrix <- do.call(cbind, distributions)
    mean_density <- rowMeans(distributions_matrix, na.rm = TRUE)

    ci_lower <- apply(distributions_matrix, 1, function(x) quantile(x, probs = (1 - prob) / 2, na.rm = TRUE))
    ci_upper <- apply(distributions_matrix, 1, function(x) quantile(x, probs = 1 - (1 - prob) / 2, na.rm = TRUE))

    data_generating_curve <- pweibull(x_values - data_gen_threshold, shape = data_gen_alpha, scale = data_gen_beta) * data_gen_asymptote

    return(list(
        mean_density = mean_density,
        data_generating_curve = data_generating_curve,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        x_values = x_values
    ))
}

# Calculating MSE for evaluation
calculate_mse <- function(estimated_curve, true_curve) {
    mse <- mean((estimated_curve - true_curve)^2)
    return(mse)
}

# Calculate 95% Credible Interval Coverage
calculate_ci_coverage <- function(true_curve, ci_lower, ci_upper) {
    coverage <- mean(true_curve >= ci_lower & true_curve <= ci_upper)
    return(coverage)
}
