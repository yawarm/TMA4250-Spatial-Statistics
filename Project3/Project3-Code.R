# Yawar Mahmood
# TMA4250
# Project 3: Gaussian Markov Random Fields

# Importing packages


library(spdep)
library(Matrix)
library(lattice)
library(ggplot2)
library(sf)
library(MASS)
library(sf)
library(viridis)

# This funciton is from functions.R, and is used for plotting the maps

## plotAreaCol
# This functions plots the values and saves the figure to the specified file.
# It is advisable to not plot directly in R since it will take a long time.
# The arguments are:
#   fNamme: file name for saving figure
#   width: width of figure in inches
#   height: height of figure in inches
#   estVal: the k values to be plotted on geoMap
#   geoMap: the map containing k regions
#   leg: name to use on top of legend
#   colLim: control lower limit and upper limit of color scale (colLim = c(lowVal, highVal))
plotAreaCol = function(fName, width, height, estVal, geoMap, leg, colLim = NULL){
  if(is.null(colLim)){
    colLim = range(estVal)
  }
  
  # Set up data object for plotting
  nigeriaMapTmp = geoMap
  nigeriaMapTmp$MCV1 = estVal
  
  # Plot
  map = ggplot() +
    geom_sf(data = nigeriaMapTmp,
            aes(fill = MCV1),
            color = 'gray', size = .2)+
    scale_fill_viridis_c(direction = 1,
                         begin = 1,
                         end = 0,
                         limit = colLim,
                         name = leg) + 
    theme(text = element_text(size=40),
          legend.key.height = unit(4, 'cm'),
          legend.key.width  = unit(1.75, 'cm'))
  ggsave(filename = fName,
         plot = map,
         width = width, 
         height = height)
}

## Problem 1

# a)

# Function to:
# Calculate precision matrix, and its Rank and Dimensions
# Calculate proportions of zeros
# Plot sparsity patterns
process_and_plot <- function(file_path, tau) {
  # Read adjacency matrix
  adj_matrix <- as.matrix(read.table(file_path, header = TRUE, row.names = 1))
  
  # Convert adjacency matrix a sparse matrix
  R <- as(adj_matrix, "CsparseMatrix")
  
  # Adjust the diagonal entries to ensure the sum of each row equals 0
  diag(R) <- -rowSums(R)
  
  # Calculate Precision matrix
  Q <- tau * R
  
  # Compute proportion of non-zero elements
  prop_non_zero <- sum(Q != 0) / length(Q)
  
  cat("Proportion of non-zero elements:", prop_non_zero, "\n")
  
  # Prepare Q for plotting
  Q_dense <- as.matrix(Q)
  
  cat("Dimensions of Q:", dim(Q_dense), "\n")
  cat("Rank of Q:", qr(Q_dense)$rank, "\n")
  
  # Plotting sparsity pattern
  image(Q_dense != 0, main = paste("Sparsity Pattern for", basename(file_path)), col = c("white", "black"), xlab = "Column", ylab = "Row")
}

process_and_plot("/Users/yawarmahmood/Downloads/Admin1Graph.txt", 1)
process_and_plot("/Users/yawarmahmood/Downloads/AdminGraph2.txt", 1)


# b)

# Function to simulate from the Besag model
simulate_besag_cholesky <- function(Q, n, sum_zero_constraint = TRUE) {
  # Add small value for numerical stability
  Q <- Q + 1e-10
  
  # Calculate the Cholesky decomposition of Q
  L <- chol(Q)
  
  # Simulate from N(0, I)
  z <- matrix(rnorm(n * nrow(Q)), ncol = n)
  
  # Backsolving L^T * x = z to get draws from N(0, Q^-1)
  simulations <- backsolve(L, z)
  
  # Sum-zero constraint
  if (sum_zero_constraint) {
    simulations <- apply(simulations, 2, function(x) x - mean(x))
  }
  
  return(simulations)
}

# Function to simulate from N(0, I_37)
simulate_normal <- function(n, dim) {
  simulations <- MASS::mvrnorm(n = n, mu = rep(0, dim), Sigma = diag(dim))
  return(simulations)
}

# Load your geoMap for admin1
geoMap <- load("/Users/yawarmahmood/Downloads/Admin1Geography.RData")

# Generate simulations
adj_matrix <- as.matrix(read.table("/Users/yawarmahmood/Downloads/Admin1Graph.txt", header = TRUE, row.names = 1))
tau_1 <- 1
R <- -adj_matrix
diag(R) <- rowSums(adj_matrix)
# Precision matrix
Q <- tau_1 * R

besag_simulations <- simulate_besag_cholesky(Q, n = 2)

normal_simulations <- simulate_normal(n = 2, dim=37)

# Plot simulations
plotAreaCol("/Users/yawarmahmood/Downloads/Besag_1_task_1b.png", 20, 20, besag_simulations[,1], nigeriaAdm1, "Besag Model, sim1", colLim = c(-3,3))
plotAreaCol("/Users/yawarmahmood/Downloads/Besag_2_task_1b.png", 20, 20, besag_simulations[,2], nigeriaAdm1, "Besag Model, sim2", colLim = c(-3,3))
plotAreaCol("/Users/yawarmahmood/Downloads/Normal_1_task_1b.png", 20, 20, normal_simulations[1,], nigeriaAdm1, "Normal Model, sim1", colLim = c(-3,3))
plotAreaCol("/Users/yawarmahmood/Downloads/Normal_2_task_1b.png", 20, 20, normal_simulations[2,], nigeriaAdm1, "Normal Model, sim2", colLim = c(-3,3))


# c)

# Load your geoMap for admin2
geoMap2 <- load("/Users/yawarmahmood/Downloads/Admin2Geography(2).RData")

# Read adjacency matrix for admin2
adj_matrix_admin2 <- as.matrix(read.table("/Users/yawarmahmood/Downloads/AdminGraph2.txt", header = TRUE, row.names = 1))

# Generate simulations
tau_2 <- 1
R_admin2 <- -adj_matrix_admin2
diag(R_admin2) <- rowSums(adj_matrix_admin2)
# Precision matrix
Q_admin2 <- tau_2 * R_admin2

besag_simulations_admin2 <- simulate_besag_cholesky(Q_admin2, n = 2)
normal_simulations_admin2 <- simulate_normal(n=2, dim=775)

# Plot simulations
plotAreaCol("/Users/yawarmahmood/Downloads/Besag_1_task_1c.png", 20, 20, besag_simulations_admin2[,1], nigeriaAdm2, "Besag Model, sim1", colLim = c(-3,3))
plotAreaCol("/Users/yawarmahmood/Downloads/Besag_2_task_1c.png", 20, 20, besag_simulations_admin2[,2], nigeriaAdm2, "Besag Model, sim2", colLim = c(-3,3))
plotAreaCol("/Users/yawarmahmood/Downloads/Normal_1_task_1c.png", 20, 20, normal_simulations_admin2[1,], nigeriaAdm2, "Normal Model, sim1", colLim = c(-3,3))
plotAreaCol("/Users/yawarmahmood/Downloads/Normal_2_task_1c.png", 20, 20, normal_simulations_admin2[2,], nigeriaAdm2, "Normal Model, sim2", colLim = c(-3,3))


# d)

# Simulate 100 realizations from Besag Model
sim_besag_admin2 <- simulate_besag_cholesky(Q_admin2, n = 100)

# Compute empirical marginal variances
empirical_variances <- apply(sim_besag_admin2, 1, var)

# Plot empirical marginal variances
plotAreaCol("/Users/yawarmahmood/Downloads/Besag_marginal_variances_1d.png", 20, 20, empirical_variances, nigeriaAdm2, "Besag Model, Empirical variances", colLim = c(-3,3))

# Extract Gubio simulations
gubio_simulations <- sim_besag_admin2[150,]

# Compute correlations between Gubio and every area
correlations_with_gubio <- apply(sim_besag_admin2, 1, function(x) cor(x, gubio_simulations))

# Plot correlations between Gubio and every area
plotAreaCol("/Users/yawarmahmood/Downloads/Besag_gubio_correlations_1d.png", 20, 20, correlations_with_gubio, nigeriaAdm2, "Besag Model, Gubio, Correlations", colLim = c(-3,3))





## Problem 2

# a)

# Extract observed proportions
DirectEstimates <- read.csv("/Users/yawarmahmood/Downloads/DirectEstimates.txt", sep="")

# inverse logit function
invlogit <- function(x) {
  1/(1+ exp(-x))
}

# logit function
logit <- function(x) {
  log(x / (1 - x))
}

# inverse logit the data to get actual observed proportions
DirectEstimates_hat <- invlogit(DirectEstimates$Observation)

# Plot actual observed proportions
plotAreaCol("/Users/yawarmahmood/Downloads/2a.png",20, 20, DirectEstimates_hat ,nigeriaAdm1, leg="p-hat", colLim = c(0,1))


# b)

# Extract observations
y <- DirectEstimates$Observation

# Variance of the vague prior
sigma_squared <- 100^2

# Extract Variance
variance_vec <- DirectEstimates$StdDev^2

# Covariance matrix for X given Y
Sigma_X_given_Y <- diag(sigma_squared * variance_vec / (sigma_squared + variance_vec))

# Mean matrix for X given Y
Mu_X_given_Y <- diag(sigma_squared / (sigma_squared + variance_vec)) %*% y

# Simulate samples from X given Y
X_samples <- mvrnorm(100, Mu_X_given_Y, Sigma_X_given_Y)

# Inverse logit the data to get actual observed proportions
P_samples <- invlogit(X_samples)

# Calculate the median and coefficient of variation
P_medians <- apply(P_samples, 2, median)
P_coeff_variance <- apply(P_samples, 2, function(x) sd(x) / mean(x))

# Plot the median and coefficient of variation
plotAreaCol("/Users/yawarmahmood/Downloads/2b_medians.png", 20, 20, P_medians, nigeriaAdm1, "Median", colLim = c(0,1))
plotAreaCol("/Users/yawarmahmood/Downloads/2b_variances.png", 20, 20, P_coeff_variance, nigeriaAdm1, "Coefficient of variation", colLim = c(0,1))


# c)

# Generate simulations
adj_matrix <- as.matrix(read.table("/Users/yawarmahmood/Downloads/Admin1Graph.txt", header = TRUE, row.names = 1))
tau_1 <- 1
D <- diag(rowSums(adj_matrix))
variance_vec <- DirectEstimates$StdDev^2
Sigma_inv <- diag(1/variance_vec)
# Precision matrix
Q <- tau_1 * (D - adj_matrix)

# Posterior Precision
posterior_precision <- Sigma_inv + Q

# Posterior Mean
posterior_mean <- solve(posterior_precision, Sigma_inv %*% y)

# Simulate samples from X given Y
set.seed(123)
post_samples <- mvrnorm(n = 100, mu = posterior_mean, Sigma = solve(posterior_precision))

# Inverse logit the data to get actual observed proportions
P_samples <- invlogit(post_samples)

# Computing the Median and the coefficient of variation for P_a given Y
P_medians <- apply(P_samples, 2, median)
P_coeff_variance <- apply(P_samples, 2, function(x) sd(x) / mean(x))

# Plot the median and coefficient of variation
plotAreaCol("/Users/yawarmahmood/Downloads/2c_medians.png", 20, 20, P_medians, nigeriaAdm1, "Median", colLim = c(0,1))
plotAreaCol("/Users/yawarmahmood/Downloads/2c_variances.png", 20, 20, P_coeff_variance, nigeriaAdm1, "Coefficient of variation", colLim = c(0,1))


# d)

# Finding row number for Kaduna and defining it
DirectEstimates
kaduna_index <- 19

# Update precision matrix 
updated_precision <- posterior_precision
updated_precision[kaduna_index, kaduna_index] <- updated_precision[kaduna_index, kaduna_index] + 1 / 0.01

# Update mean 
updated_mean <- posterior_mean
updated_mean[kaduna_index] <- updated_mean[kaduna_index] + (1 / 0.01) * logit(0.5)
updated_mean <- solve(updated_precision, updated_mean)

# Using updated mean and precision matrix to run 100 realizations
updated_samples <- mvrnorm(100, updated_mean, solve(updated_precision))

# Inverse logit the data to get actual observed proportions
P_samples <- invlogit(updated_samples)

# Computing the Median and the coefficient of variation for P_a given Y
P_medians <- apply(P_samples, 2, median)
P_coeff_variance <- apply(P_samples, 2, function(x) sd(x) / mean(x))

# Plot the median and coefficient of variation
plotAreaCol("/Users/yawarmahmood/Downloads/2d_medians.png", 20, 20, P_medians, nigeriaAdm1, "Median", colLim = c(0,1))
plotAreaCol("/Users/yawarmahmood/Downloads/2d_variances.png", 20, 20, P_coeff_variance, nigeriaAdm1, "Coefficient of variation", colLim = c(0,1))


# e)

# Generate simulations
adj_matrix <- as.matrix(read.table("/Users/yawarmahmood/Downloads/Admin1Graph.txt", header = TRUE, row.names = 1))
D <- diag(rowSums(adj_matrix))
Sigma_inv <- diag(1/variance_vec)

# Define values for precision parameters
tau_01 <- 0.1
tau_1 <- 1
tau_10 <- 10

# Precesion matrices
Q_01 <- tau_01 * (D - adj_matrix)
Q_1 <- tau_1 * (D - adj_matrix)
Q_10 <- tau_10 * (D - adj_matrix)

# Posterior Precision matrices
posterior_precision_01 <- Sigma_inv + Q_01
posterior_precision_1 <- Sigma_inv + Q_1
posterior_precision_10 <- Sigma_inv + Q_10

# Posterior mean
posterior_mean_01 <- solve(posterior_precision_01, Sigma_inv %*% y)
posterior_mean_1 <- solve(posterior_precision_1, Sigma_inv %*% y)
posterior_mean_10 <- solve(posterior_precision_10, Sigma_inv %*% y)

# Simulate samples 
set.seed(123)
post_samples_01 <- mvrnorm(n = 100, mu = posterior_mean_01, Sigma = solve(posterior_precision_01))
post_samples_1 <- mvrnorm(n = 100, mu = posterior_mean_1, Sigma = solve(posterior_precision_1))
post_samples_10 <- mvrnorm(n = 100, mu = posterior_mean_10, Sigma = solve(posterior_precision_10))

# Inverse logit the data to get actual observed proportions
P_samples_01 <- invlogit(post_samples_01)
P_samples_1 <- invlogit(post_samples_1)
P_samples_10 <- invlogit(post_samples_10)

# Computing the Median and the coefficient of variation for P_a given Y
P_medians <- apply(P_samples_01, 2, median)
P_coeff_variance <- apply(P_samples_01, 2, function(x) sd(x) / mean(x))

# Plot the median and coefficient of variation
plotAreaCol("/Users/yawarmahmood/Downloads/2e01_medians.png", 20, 20, P_medians, nigeriaAdm1, "Median", colLim = c(0,1))
plotAreaCol("/Users/yawarmahmood/Downloads/2e01_variances.png", 20, 20, P_coeff_variance, nigeriaAdm1, "Coefficient of variation", colLim = c(0,1))

# Computing the Median and the coefficient of variation for P_a given Y
P_medians <- apply(P_samples_1, 2, median)
P_coeff_variance <- apply(P_samples_1, 2, function(x) sd(x) / mean(x))

# Plot the median and coefficient of variation
plotAreaCol("/Users/yawarmahmood/Downloads/2e1_medians.png", 20, 20, P_medians, nigeriaAdm1, "Median", colLim = c(0,1))
plotAreaCol("/Users/yawarmahmood/Downloads/2e1_variances.png", 20, 20, P_coeff_variance, nigeriaAdm1, "Coefficient of variation", colLim = c(0,1))

# Computing the Median and the coefficient of variation for P_a given Y
P_medians <- apply(P_samples_10, 2, median)
P_coeff_variance <- apply(P_samples_10, 2, function(x) sd(x) / mean(x))

# Plot the median and coefficient of variation
plotAreaCol("/Users/yawarmahmood/Downloads/2e10_medians.png", 20, 20, P_medians, nigeriaAdm1, "Median", colLim = c(0,1))
plotAreaCol("/Users/yawarmahmood/Downloads/2e10_variances.png", 20, 20, P_coeff_variance, nigeriaAdm1, "Coefficient of variation", colLim = c(0,1))


# f)

# Generate simulations
adj_matrix <- as.matrix(read.table("/Users/yawarmahmood/Downloads/Admin1Graph.txt", header = TRUE, row.names = 1))
D <- diag(rowSums(adj_matrix))
Sigma_inv <- diag(1/variance_vec)
R <- (D - adj_matrix)
y <- DirectEstimates$Observation

# Define the log-likelihood function
log_likelihood <- function(tau) {
  Q <- tau*R
  
  posterior_precision <- Sigma_inv + Q
  
  posterior_mean <- solve(posterior_precision, Sigma_inv %*% y)
  
  x <- mvrnorm(n = 1, mu = posterior_mean, Sigma = solve(posterior_precision))
  
  log_likelihood_value <- 36/2*log(tau) - tau * 0.5 * t(x)%*%R%*%x - 0.5*t(y - x)%*%solve(D)%*%(y-x) - 0.5*log(det(posterior_precision)) + 0.5*t(x - posterior_mean)%*%posterior_precision%*%(x - posterior_mean)
  return(-log_likelihood_value)
}

# Storage for multiple tau_hat values
tau_hats <- numeric(10)

# Simulate multiple tau_hat values and store them
set.seed(123)
for (i in 1:10) {
  optimize_result <- optimize(log_likelihood, interval = c(0.1, 10), maximum = TRUE)
  tau_hats[i] <- optimize_result$maximum
}

# Use the average of stored tau_hat values as the estimator for tau_hat
tau_hat <- mean(tau_hats)
tau_hat

# Generate simulations
adj_matrix <- as.matrix(read.table("/Users/yawarmahmood/Downloads/Admin1Graph.txt", header = TRUE, row.names = 1))
D <- diag(rowSums(adj_matrix))
Sigma_inv <- diag(1/variance_vec)
# Precision matrix
Q_hat <- tau_hat * (D - adj_matrix)

# Posterior precision matrix
posterior_precision_hat <- Sigma_inv + Q_hat

# Posterior mean
posterior_mean_hat <- solve(posterior_precision_hat, Sigma_inv %*% y)

# Simulate 100 samples for P_a given Y
set.seed(123)
post_samples_hat <- mvrnorm(n = 100, mu = posterior_mean_hat, Sigma = solve(posterior_precision_hat))

# Inverse logit the data to get actual observed proportions
P_samples_hat <- invlogit(post_samples_hat)

# Computing the Median and the coefficient of variation for P_a given Y
P_medians_hat <- apply(P_samples_hat, 2, median)
P_coeff_variance_hat <- apply(P_samples_hat, 2, function(x) sd(x) / mean(x))

# Plot the median and coefficient of variation
plotAreaCol("/Users/yawarmahmood/Downloads/2f_medians.png", 20, 20, P_medians_hat, nigeriaAdm1, "Median", colLim = c(0,1))
plotAreaCol("/Users/yawarmahmood/Downloads/2f_variances.png", 20, 20, P_coeff_variance_hat, nigeriaAdm1, "Coefficient of variation", colLim = c(0,1))
