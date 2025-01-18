#Load libraries
library(ggplot2)
library(spatial)
library(reshape2)
library(gridExtra)
library(dplyr)
library(spatstat)


# PROBLEM 1

cells = ppinit('/Users/yawarmahmood/Downloads/pp_cells.dat')
redwood = ppinit('/Users/yawarmahmood/Downloads/pp_redwood.dat')
pines = ppinit('/Users/yawarmahmood/Downloads/pp_pines.dat')




#a) ########

#Define data frames
cells_df = data.frame(x = cells$x, y = cells$y)
redwood_df = data.frame(x = redwood$x, y = redwood$y)
pines_df = data.frame(x = pines$x, y = pines$y)

# Plot each dataframe
ggplot(data = cells_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Cells data") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data = redwood_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Redwood data") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data = pines_df, aes(x=x, y=y)) + 
  geom_point() + 
  ggtitle("Pines data") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))




#b) ########

# Convert the data frames to ppp objects
window <- owin(xrange=c(0, 1), yrange=c(0, 1))
cells_ppp <- ppp(cells_df$x, cells_df$y, window=window)
redwood_ppp <- ppp(redwood_df$x, redwood_df$y, window=window)
pines_ppp <- ppp(pines_df$x, pines_df$y, window=window)

# Compute the empirical L-function for each point pattern
L_cells <- envelope(cells_ppp, fun="Lest", nsim = 99)
L_redwood <- envelope(redwood_ppp, fun="Lest", nsim = 99)
L_pines <- envelope(pines_ppp, fun="Lest", nsim = 99)

# Plot the empirical L-functions with theoretical curve
plot(L_cells, main="L-function for Cells Data")
plot(L_redwood, main="L-function for Redwood Data")
plot(L_pines, main="L-function for Pines Data")




#c) #####

simulate_and_plot_L <- function(df, title) {
  window <- owin(xrange=c(0, 1), yrange=c(0, 1))
  observed_ppp <- ppp(df$x, df$y, window=window)
  
  # Simulate Poisson processes
  n_sims <- 100
  L_vals <- vector("list", n_sims)
  for (i in 1:n_sims) {
    sim_ppp <- rpoispp(lambda = observed_ppp$n/area(observed_ppp$window), win=observed_ppp$window)
    L_vals[[i]] <- Lest(sim_ppp)$iso
  }
  
  # Calculate prediction intervals
  L_combined <- do.call(cbind, L_vals)
  lower <- apply(L_combined, 1, quantile, probs=0.05)
  upper <- apply(L_combined, 1, quantile, probs=0.95)
  r <- Lest(observed_ppp)$r
  observed_L <- Lest(observed_ppp)$iso
  
  # Plot
  ggplot() + 
    geom_ribbon(aes(x=r, ymin=lower, ymax=upper), fill="grey80") +
    geom_line(aes(x=r, y=observed_L), color="blue") +
    ggtitle(title) +
    xlab("Distance r") + 
    ylab(expression(L(r) - r)) +
    theme_minimal()
}

# Plot data
plot_cells <- simulate_and_plot_L(cells_df, "Cells Data")
plot_redwood <- simulate_and_plot_L(redwood_df, "Redwood Data")
plot_pines <- simulate_and_plot_L(pines_df, "Pines Data")

# Display Plots
grid.arrange(plot_cells, plot_redwood, plot_pines, ncol = 1)



# PROBLEM 2

#a) #####

# Load detection probabilities
obs_prob <- read.table("/Users/yawarmahmood/Downloads/obsprob.txt", header=TRUE)

# Load observed counts
obs_counts <- read.table("/Users/yawarmahmood/Downloads/obspines.txt", header=TRUE)

par(mfrow = c(1, 2))

# Plot probabilities and counts
ggplot(obs_prob, aes(x, y, fill = alpha)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Detection Probabilities (αij)", x = "X coordinate", y = "Y coordinate", fill = "Prob")

ggplot(obs_counts, aes(x, y, fill = N_obs)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Counts of Detected Pine Trees (Mij)", x = "X coordinate", y = "Y coordinate", fill = "Count")



#b) #####

# No Code

#c) ####

# Calculate the total observed counts
total_M <- sum(obs_counts$N_obs)

# Calculate the average detection probability
avg_alpha <- mean(obs_prob$alpha)

# Calculate the correction factor C
C <- 1 / avg_alpha

# Calculate the unbiased estimator for lambda
A <- 90000
lambda_hat <- (C * total_M) / A

# Generate realizations and simulate point patterns
set.seed(123)
simulations <- list()

for (i in 1:3) {
  # Generate realizations of the true counts N using the estimated lambda
  lambda_realization <- rpois(n = nrow(obs_counts), lambda = lambda_hat * 100)
  
  # Simulate point patterns within each grid cell
  points <- lapply(1:length(lambda_realization), function(idx) {
    if (lambda_realization[idx] > 0) {
      x_coords <- runif(lambda_realization[idx], min = obs_counts$x[idx] - 5, max = obs_counts$x[idx] + 5)
      y_coords <- runif(lambda_realization[idx], min = obs_counts$y[idx] - 5, max = obs_counts$y[idx] + 5)
      return(data.frame(x = x_coords, y = y_coords))
    } else {
      return(data.frame(x = numeric(0), y = numeric(0)))
    }
  })
  
  simulation_points <- do.call(rbind, points)
  simulations[[i]] <- simulation_points
}

# Function to plot each simulation
plot_simulation <- function(simulated_points, title) {
  ggplot(simulated_points, aes(x = x, y = y)) +
    geom_point(alpha = 0.6, color = "blue") +
    labs(title = title, x = "X coordinate", y = "Y coordinate") +
    theme_minimal() +
    coord_fixed()
}

# Plot
for (i in 1:length(simulations)) {
  print(plot_simulation(simulations[[i]], paste("Simulation", i)))
}




#d) ####


# Function to simulate point patterns within grid cells based on posterior assumptions
simulate_posterior_patterns <- function(obs_counts, obs_prob, lambda_hat, num_simulations = 3) {
  set.seed(123) 
  
  # Calculate adjusted lambda for each cell
  adjusted_lambda <- lambda_hat / obs_prob$alpha
  
  simulations <- list()
  
  for (i in 1:num_simulations) {
    points <- lapply(1:nrow(obs_counts), function(idx) {
      # Generate realizations of the true counts N using the adjusted lambda
      n_points <- rpois(1, lambda = adjusted_lambda[idx] * 100)
      
      # Simulate point locations within the cell
      if (n_points > 0) {
        x_coords <- runif(n_points, min = obs_counts$x[idx] - 5, max = obs_counts$x[idx] + 5)
        y_coords <- runif(n_points, min = obs_counts$y[idx] - 5, max = obs_counts$y[idx] + 5)
        return(data.frame(x = x_coords, y = y_coords))
      } else {
        return(data.frame(x = numeric(0), y = numeric(0)))
      }
    })
    
    simulation_points <- do.call(rbind, points)
    simulations[[i]] <- simulation_points
  }
  
  simulations
}

# Plot
plot_simulated_patterns <- function(simulations) {
  plots <- lapply(1:length(simulations), function(i) {
    ggplot(simulations[[i]], aes(x = x, y = y)) +
      geom_point(alpha = 0.6, color = "blue") +
      labs(title = paste("Posterior Simulation", i), x = "X coordinate", y = "Y coordinate") +
      theme_minimal() +
      coord_fixed()
  })
  
  plots
}

# Generate and plot the simulations
simulations <- simulate_posterior_patterns(obs_counts, obs_prob, lambda_hat)
plot_list <- plot_simulated_patterns(simulations)

for (plot in plot_list) {
  print(plot)
}



#e) ####

num_realizations <- 500
grid_size <- 30 
cell_area <- 100 

# Initialize matrices
prior_sum <- matrix(0, nrow = grid_size, ncol = grid_size)
prior_sum_sq <- matrix(0, nrow = grid_size, ncol = grid_size)
posterior_sum <- matrix(0, nrow = grid_size, ncol = grid_size)
posterior_sum_sq <- matrix(0, nrow = grid_size, ncol = grid_size)

# Averaged detection probability
avg_alpha <- mean(obs_prob$alpha)

# Loop to generate realizations
for (i in 1:num_realizations) {

  prior_N <- rpois(n = grid_size^2, lambda = lambda_hat * cell_area)
  prior_sum <- prior_sum + matrix(prior_N, nrow = grid_size, byrow = TRUE)
  prior_sum_sq <- prior_sum_sq + (matrix(prior_N, nrow = grid_size, byrow = TRUE))^2
  
  posterior_N <- rpois(n = grid_size^2, lambda = lambda_hat * cell_area / avg_alpha)
  posterior_sum <- posterior_sum + matrix(posterior_N, nrow = grid_size, byrow = TRUE)
  posterior_sum_sq <- posterior_sum_sq + (matrix(posterior_N, nrow = grid_size, byrow = TRUE))^2
}

# Compute averages and standard deviations
prior_mean <- prior_sum / num_realizations
posterior_mean <- posterior_sum / num_realizations
prior_sd <- sqrt(prior_sum_sq / num_realizations - prior_mean^2)
posterior_sd <- sqrt(posterior_sum_sq / num_realizations - posterior_mean^2)

# Function to plot matrices
plot_matrix <- function(mat, title, low_color = "white", high_color = "red") {
  melted_mat <- melt(mat)
  ggplot(data = melted_mat, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = low_color, high = high_color) +
    labs(title = title, x = "Row", y = "Column") +
    theme_minimal() +
    coord_fixed()
}

# Plotting
plot1 <- plot_matrix(prior_mean, "A Priori Expected Value (E[N])")
plot2 <- plot_matrix(posterior_mean, "A Posteriori Expected Value (E[N|M = m])")
plot3 <- plot_matrix(prior_sd, "A Priori Standard Deviation")
plot4 <- plot_matrix(posterior_sd, "A Posteriori Standard Deviation")

gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)

par(mfrow = c(1, 1))




# PROBLEM 3

#a) #####


# Simulation function for Neyman-Scott process
neumanScottProcess <- function(parent_intensity, offspring_dispersion, mean_offspring_per_parent, area, original_data_count) {
  offspring_points <- data.frame(x = numeric(), y = numeric())
  
  # Continue generating until we reach the desired number of points
  while(nrow(offspring_points) < original_data_count) {
    # Generate parent point
    parent_x <- runif(1, 0, 1)
    parent_y <- runif(1, 0, 1)
    
    # Determine number of offspring for this parent
    num_offspring <- rpois(1, mean_offspring_per_parent)
    
    # Generate offspring points around parent
    offspring_x <- rnorm(num_offspring, mean = parent_x, sd = offspring_dispersion)
    offspring_y <- rnorm(num_offspring, mean = parent_y, sd = offspring_dispersion)
    
    valid_indices <- which(offspring_x >= 0 & offspring_x <= 1 & offspring_y >= 0 & offspring_y <= 1)
    valid_offspring <- data.frame(x = offspring_x[valid_indices], y = offspring_y[valid_indices])
    
    # Add valid offspring points to offspring_points
    offspring_points <- rbind(offspring_points, valid_offspring)
    
    if (nrow(offspring_points) > original_data_count) {
      offspring_points <- offspring_points[1:original_data_count, ]
    }
  }
  
  return(ppp(offspring_points$x, offspring_points$y, window = owin(c(0,1), c(0,1))))
}

nsim <- 100
area <- 1

# Define parameters based on empirical fit
parent_intensity <- 0.05  # lambda
mean_offspring_per_parent <- 6  # mu
offspring_dispersion <- 0.1  # sigma

realizations <- lapply(1:3, function(x) neumanScottProcess(parent_intensity, offspring_dispersion, mean_offspring_per_parent, 1, 62))

# Plotting function for point patterns
plot_ppp <- function(ppp_data, title) {
  plot(ppp_data, main=title, xlab="x", ylab="y", pch=19, cex=0.5)
}

empty_plot <- function() {
  plot(1, type="n", axes=FALSE, xlab="", ylab="", main="")
}

# Display the redwood tree dataset and the three model realizations
par(mfrow=c(2,2), mar=c(4,4,2,1))

plot_ppp(redwood_ppp, "Redwood Data")
for (i in 1:length(realizations)) {
  plot_ppp(realizations[[i]], paste("Realization", i))
}

par(mfrow=c(1,1), mar=c(5,4,4,2) + 0.1)

window <- owin(xrange=c(0, 1), yrange=c(0, 1))

# Compute the empirical L-function for each point pattern
sim1 <- envelope(realizations[[1]], fun="Lest", nsim = 99)
sim2 <- envelope(realizations[[2]], fun="Lest", nsim = 99)
sim3 <- envelope(realizations[[2]], fun="Lest", nsim = 99)

par(mfrow=c(2,2))
# Plot the empirical L-functions with theoretical curve
plot(L_redwood, main="L-function for Redwood Data")
plot(sim1, main="L-function for sim 1")
plot(sim2, main="L-function for sim 2")
plot(sim3, main="L-function for sim 3")






# PROBLEM 4
par(mfrow=c(1,1))

#a) #####

# Parameters based on empirical guesses
lambda_guess <- 40  # number of points per unit area
r0_guess <- 0.1 # interaction radius
beta_guess <- 2  # interaction strength
gamma_guess <- exp(-beta_guess)

# Window of observation
win <- owin(c(0,1), c(0,1))

# Empirical L-function of the dataset
L_empirical <- Lest(cells_ppp)

# Simulate 100 realizations
nsim <- 100
L_sims <- vector("list", nsim)
for (i in 1:nsim) {
  sim_pattern <- rStrauss(W = win, lambda_guess, gamma_guess, r0_guess)
  L_sims[[i]] <- Lest(sim_pattern)
}

# Compute the 90% prediction intervals
r <- L_sims[[1]]$r
L_values <- sapply(L_sims, `[[`, "iso")
L_mean <- apply(L_values, 1, mean)
L_lower <- apply(L_values, 1, quantile, probs = 0.05)
L_upper <- apply(L_values, 1, quantile, probs = 0.95)

# Plotting
plot(r, L_empirical$iso, type = "l", col = "red", xlab = "Distance r", ylab = "L-function", main = "Empirical L-Function vs. Simulated")
lines(r, L_mean, col = "blue")
matlines(r, t(apply(L_values, 1, quantile, probs = c(0.05, 0.95))), lty = c(2,2), col = "blue")
legend("topleft", legend = c("Empirical", "Mean Simulated", "90% Prediction Interval"), col = c("red", "blue", "blue"), lty = c(1, 1, 2))

