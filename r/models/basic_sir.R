# Basic SIR Model for Malaria Transmission
# ========================================
#
# This script implements a basic Susceptible-Infected-Recovered (SIR) model
# specifically adapted for malaria transmission dynamics.
#
# Key malaria-specific considerations:
# - Recovery rate accounts for natural immunity and treatment
# - Transmission rate includes vector-borne transmission
# - Basic reproduction number (R0) calculation
# - Population dynamics with birth and death rates

# Load required libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Define the MalariaSIRModel class (using R's S3 system)
MalariaSIRModel <- function(beta = 0.3,      # Transmission rate (per day)
                           gamma = 0.1,     # Recovery rate (per day)
                           mu = 0.0001,     # Birth/death rate (per day)
                           N = 10000,       # Total population
                           I0 = 10,         # Initial infected
                           R0 = 0) {        # Initial recovered
  
  # Calculate initial susceptible population
  S0 <- N - I0 - R0
  
  # Calculate basic reproduction number
  R0_basic <- beta / (gamma + mu)
  
  # Create model object
  model <- list(
    beta = beta,
    gamma = gamma,
    mu = mu,
    N = N,
    S0 = S0,
    I0 = I0,
    R0 = R0,
    R0_basic = R0_basic
  )
  
  class(model) <- "MalariaSIRModel"
  return(model)
}

# SIR differential equations
sir_equations <- function(t, y, params) {
  with(as.list(c(y, params)), {
    # SIR equations with birth and death
    dS <- mu * N - beta * S * I / N - mu * S
    dI <- beta * S * I / N - (gamma + mu) * I
    dR <- gamma * I - mu * R
    
    return(list(c(dS, dI, dR)))
  })
}

# Simulate the SIR model
simulate.MalariaSIRModel <- function(model, t_span = c(0, 365), num_points = 1000) {
  # Create time vector
  times <- seq(from = t_span[1], to = t_span[2], length.out = num_points)
  
  # Initial conditions
  y0 <- c(S = model$S0, I = model$I0, R = model$R0)
  
  # Parameters for the differential equation
  params <- list(
    beta = model$beta,
    gamma = model$gamma,
    mu = model$mu,
    N = model$N
  )
  
  # Solve the differential equations
  solution <- ode(y = y0, times = times, func = sir_equations, parms = params)
  
  return(solution)
}

# Calculate epidemic metrics
get_epidemic_metrics <- function(model, solution) {
  # Extract the solution data
  S <- solution[, "S"]
  I <- solution[, "I"]
  R <- solution[, "R"]
  
  # Calculate metrics
  peak_infection <- max(I)
  peak_time <- which.max(I)
  final_recovered <- R[length(R)]
  attack_rate <- final_recovered / model$N
  
  metrics <- list(
    R0_basic = model$R0_basic,
    peak_infection = peak_infection,
    peak_infection_rate = peak_infection / model$N,
    final_attack_rate = attack_rate,
    time_to_peak = peak_time
  )
  
  return(metrics)
}

# Plot the SIR model results
plot_results <- function(model, solution, save_path = NULL) {
  # Convert solution to data frame for easier plotting
  df <- as.data.frame(solution)
  df$time <- df$time
  
  # Create plots
  p1 <- ggplot(df, aes(x = time)) +
    geom_line(aes(y = S, color = "Susceptible"), size = 1) +
    geom_line(aes(y = I, color = "Infected"), size = 1) +
    geom_line(aes(y = R, color = "Recovered"), size = 1) +
    labs(title = "Malaria SIR Model Dynamics",
         x = "Time (days)",
         y = "Number of individuals",
         color = "Compartment") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Infection rate plot
  df$infection_rate <- df$I / model$N * 100
  p2 <- ggplot(df, aes(x = time, y = infection_rate)) +
    geom_line(color = "red", size = 1) +
    labs(title = "Infection Rate Over Time",
         x = "Time (days)",
         y = "Infection Rate (%)") +
    theme_minimal()
  
  # Phase plane plot
  p3 <- ggplot(df, aes(x = S, y = I)) +
    geom_path(color = "purple", size = 1) +
    labs(title = "Phase Plane (S vs I)",
         x = "Susceptible",
         y = "Infected") +
    theme_minimal()
  
  # Cumulative cases plot
  df$cumulative_cases <- df$R + df$I
  p4 <- ggplot(df, aes(x = time, y = cumulative_cases)) +
    geom_line(color = "orange", size = 1) +
    labs(title = "Cumulative Malaria Cases",
         x = "Time (days)",
         y = "Cumulative Cases") +
    theme_minimal()
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = 12, height = 8, dpi = 300)
  }
  
  return(combined_plot)
}

# Parameter sensitivity analysis
parameter_sensitivity <- function(model, param_name, param_values, t_span = c(0, 365)) {
  results <- data.frame()
  
  for (value in param_values) {
    # Create temporary model with modified parameter
    temp_model <- MalariaSIRModel(
      beta = ifelse(param_name == "beta", value, model$beta),
      gamma = ifelse(param_name == "gamma", value, model$gamma),
      mu = ifelse(param_name == "mu", value, model$mu),
      N = model$N,
      I0 = model$I0,
      R0 = model$R0
    )
    
    # Simulate
    solution <- simulate(temp_model, t_span)
    metrics <- get_epidemic_metrics(temp_model, solution)
    
    # Add to results
    results <- rbind(results, data.frame(
      parameter = param_name,
      value = value,
      R0_basic = metrics$R0_basic,
      peak_infection = metrics$peak_infection,
      peak_infection_rate = metrics$peak_infection_rate,
      final_attack_rate = metrics$final_attack_rate,
      time_to_peak = metrics$time_to_peak
    ))
  }
  
  return(results)
}

# Main function to demonstrate the malaria SIR model
main <- function() {
  cat("Malaria SIR Model Simulation\n")
  cat("============================\n\n")
  
  # Create model with malaria-specific parameters
  model <- MalariaSIRModel(
    beta = 0.4,      # Higher transmission rate for malaria
    gamma = 0.05,    # Slower recovery (20 days average)
    mu = 0.0001,     # Low birth/death rate
    N = 10000,       # Population of 10,000
    I0 = 50,         # Start with 50 infected
    R0 = 0           # No initially recovered
  )
  
  cat("Model Parameters:\n")
  cat("  Transmission rate (β):", model$beta, "\n")
  cat("  Recovery rate (γ):", model$gamma, "\n")
  cat("  Birth/death rate (μ):", model$mu, "\n")
  cat("  Population size (N):", model$N, "\n")
  cat("  Basic reproduction number (R₀):", round(model$R0_basic, 2), "\n\n")
  
  # Simulate the model
  solution <- simulate(model, t_span = c(0, 1000), num_points = 1000)
  
  # Calculate and display metrics
  metrics <- get_epidemic_metrics(model, solution)
  cat("Epidemic Metrics:\n")
  cat("  R0_basic:", round(metrics$R0_basic, 3), "\n")
  cat("  Peak infection:", round(metrics$peak_infection, 0), "\n")
  cat("  Peak infection rate:", round(metrics$peak_infection_rate, 3), "\n")
  cat("  Final attack rate:", round(metrics$final_attack_rate, 3), "\n")
  cat("  Time to peak:", round(metrics$time_to_peak, 0), "days\n\n")
  
  # Plot results
  plot_results(model, solution, save_path = "../results/plots/malaria_sir_basic.png")
  
  # Parameter sensitivity analysis
  cat("Parameter Sensitivity Analysis:\n")
  cat("------------------------------\n")
  
  # Test different transmission rates
  beta_values <- c(0.2, 0.3, 0.4, 0.5, 0.6)
  beta_sensitivity <- parameter_sensitivity(model, "beta", beta_values)
  cat("Transmission Rate (β) Sensitivity:\n")
  print(round(beta_sensitivity[, c("value", "R0_basic", "peak_infection_rate", "final_attack_rate")], 3))
  cat("\n")
  
  # Test different recovery rates
  gamma_values <- c(0.02, 0.05, 0.1, 0.15, 0.2)
  gamma_sensitivity <- parameter_sensitivity(model, "gamma", gamma_values)
  cat("Recovery Rate (γ) Sensitivity:\n")
  print(round(gamma_sensitivity[, c("value", "R0_basic", "peak_infection_rate", "final_attack_rate")], 3))
  
  return(list(model = model, solution = solution, metrics = metrics))
}

# Run the main function if script is executed directly
if (interactive()) {
  results <- main()
}


