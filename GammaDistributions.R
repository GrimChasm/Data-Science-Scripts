# Simulating Gamma Distributions
# Simulate N = 10^5 samples from gamma(n, nλ) for n = 20 and n = 100

# Set random seed for reproducibility
set.seed(123)

# Parameters
rate <- 5  # λ = 5
N <- 10^5  # Number of samples
sample_sizes <- c(20, 100)  # n values to test

# Create a list to store the simulated data
gamma_samples <- list()

# Generate N samples from gamma(n, nλ) for each sample size
for (n in sample_sizes) {
  # For gamma(shape, rate), we need:
  # shape = n
  # rate = nλ
  shape <- n
  scaled_rate <- n * rate
  
  # Generate N random variables from the gamma distribution
  gamma_samples[[as.character(n)]] <- rgamma(N, shape = shape, rate = scaled_rate)
}

# Calculate theoretical properties
theoretical_properties <- data.frame(
  n = sample_sizes,
  mean = 1/rate,  # Mean = 1/λ for all n
  sd = 1/(rate * sqrt(sample_sizes)),  # SD = 1/(λ√n)
  skewness = 2/sqrt(sample_sizes),  # Skewness = 2/√shape
  kurtosis_excess = 6/sample_sizes  # Kurtosis excess = 6/shape
)

print("Theoretical properties of gamma(n, nλ) distributions:")
print(theoretical_properties)

# Safely create histograms with theoretical curves overlaid
tryCatch({
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
  
  for (i in 1:length(sample_sizes)) {
    n <- sample_sizes[i]
    samples <- gamma_samples[[as.character(n)]]
    
    # Calculate empirical statistics
    emp_mean <- mean(samples)
    emp_sd <- sd(samples)
    
    # Create histogram
    hist(samples, 
         breaks = 100, 
         prob = TRUE, 
         main = paste("Histogram of gamma(", n, ", ", n*rate, ")", sep=""),
         xlab = "Value", 
         xlim = c(0.05, 0.35),
         col = "lightblue",
         border = "white")
    
    # Add gamma density curve
    x <- seq(0.05, 0.35, length.out = 1000)
    lines(x, dgamma(x, shape = n, rate = n*rate), 
          col = "darkblue", lwd = 2)
    
    # Add normal approximation curve
    lines(x, dnorm(x, mean = 1/rate, sd = 1/(rate*sqrt(n))), 
          col = "red", lwd = 2, lty = 2)
    
    # Add vertical line at mean
    abline(v = 1/rate, col = "darkgreen", lwd = 2)
    
    # Add legend
    legend("topright", 
           legend = c("Gamma density", "Normal approximation", "Mean (1/λ)"),
           col = c("darkblue", "red", "darkgreen"), 
           lwd = 2, 
           lty = c(1, 2, 1))
    
    # Calculate text y-position more reliably
    y_values <- dgamma(x, shape = n, rate = n*rate)
    max_y <- max(y_values, na.rm = TRUE)
    text_x <- 0.3
    text_y <- 0.9 * max_y
    
    # Fix newlines in text
    text(text_x, text_y, 
         paste("Mean: ", round(emp_mean, 5), 
               "\nSD: ", round(emp_sd, 5),
               "\nSkewness: ", round(2/sqrt(n), 5),
               "\nExcess Kurtosis: ", round(6/n, 5)),
         adj = 1)
  }
  
  # Reset layout
  par(mfrow = c(1, 1))
  
  # Create QQ plots to check for normality
  par(mfrow = c(1, 2))
  
  for (i in 1:length(sample_sizes)) {
    n <- sample_sizes[i]
    samples <- gamma_samples[[as.character(n)]]
    
    # Create QQ plot
    qqnorm(samples, 
           main = paste("Q-Q Plot for n =", n),
           pch = 19, 
           col = "blue")
    qqline(samples, col = "red", lwd = 2)
  }
  
  # Reset layout
  par(mfrow = c(1, 1))
  
  # Create a comparative plot showing both distributions
  plot(density(gamma_samples[["20"]]), 
       main = "Comparison of Gamma Densities", 
       xlab = "Value",
       ylim = c(0, 25),
       xlim = c(0.05, 0.35),
       col = "blue",
       lwd = 2)
  
  lines(density(gamma_samples[["100"]]), 
        col = "red",
        lwd = 2)
  
  # Add the theoretical densities
  x <- seq(0.05, 0.35, length.out = 1000)
  lines(x, dgamma(x, shape = 20, rate = 20*rate), 
        col = "darkblue", lwd = 2, lty = 2)
  lines(x, dgamma(x, shape = 100, rate = 100*rate), 
        col = "darkred", lwd = 2, lty = 2)
  
  # Add vertical line at mean
  abline(v = 1/rate, col = "black", lwd = 2, lty = 3)
  
  # Add legend
  legend("topright", 
         legend = c("n = 20 (empirical)", "n = 100 (empirical)", 
                    "n = 20 (theoretical)", "n = 100 (theoretical)", "Mean (1/λ)"),
         col = c("blue", "red", "darkblue", "darkred", "black"), 
         lwd = 2, 
         lty = c(1, 1, 2, 2, 3))
}, error = function(e) {
  cat("Error in plotting:", conditionMessage(e), "\n")
})

# Calculate and print quantiles for the distributions
quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99)
quantile_results <- matrix(nrow = length(quantiles), ncol = length(sample_sizes))
rownames(quantile_results) <- paste0(quantiles*100, "%")
colnames(quantile_results) <- paste0("n=", sample_sizes)

for (i in 1:length(sample_sizes)) {
  n <- sample_sizes[i]
  samples <- gamma_samples[[as.character(n)]]
  quantile_results[, i] <- quantile(samples, probs = quantiles)
}

cat("\nEmpirical quantiles:\n")
print(quantile_results)

cat("\nTheoretical quantiles (from gamma distribution):\n")
theoretical_quantiles <- matrix(nrow = length(quantiles), ncol = length(sample_sizes))
rownames(theoretical_quantiles) <- paste0(quantiles*100, "%")
colnames(theoretical_quantiles) <- paste0("n=", sample_sizes)


for (i in 1:length(sample_sizes)) {
  n <- sample_sizes[i]
  theoretical_quantiles[, i] <- qgamma(quantiles, shape = n, rate = n*rate)
}

print(theoretical_quantiles)


