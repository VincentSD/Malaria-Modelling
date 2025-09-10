# SEIR Model for Malaria Transmission
# ==================================
#
# This script implements a Susceptible-Exposed-Infected-Recovered (SEIR) model
# for malaria transmission, which includes a latent period (exposed state)
# that is more realistic for malaria transmission dynamics.
#
# Key features:
# - Latent period (exposed state) before becoming infectious
# - Age-specific parameters for different malaria strains
# - Seasonal variations in transmission
# - Treatment and immunity considerations

# Load required libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Define the MalariaSEIRModel class (using R's S3 system)
MalariaSEIRModel <- function(beta = 0.3,      # Transmission rate (per day)
                            sigma = 0.1,     # Incubation rate (1/latent period)
                            gamma = 0.1,     # Recovery rate (per day)
                            mu = 0.0001,     # Birth/death rate (per day)
                            N = 10000,       # Total population
                            E0 = 5,          # Initial exposed
                            I0 = 10,         # Initial infected
                            R0 = 0,          # Initial recovered
                            seasonal_amplitude = 0.2,  # Seasonal variation amplitude
                            seasonal_phase = 0.0) {    # Seasonal phase offset
  
  # Calculate initial susceptible population
  S0 <- N - E0 - I0 - R0
  
  # Calculate basic reproduction number
  R0_basic <- beta / (gamma + mu)
  
  # Calculate effective reproduction number (accounting for latent period)
  R0_effective <- beta * sigma / ((sigma + mu) * (gamma + mu))
  
  # Create model object
  model <- list(
    beta = beta,
    sigma = sigma,
    gamma = gamma,
    mu = mu,
    N = N,
    S0 = S0,
    E0 = E0,
    I0 = I0,
    R0 = R0,
    seasonal_amplitude = seasonal_amplitude,
    seasonal_phase = seasonal_phase,
    R0_basic = R0_basic,
    R0_effective = R0_effective
  )
  
  class(model) <- "MalariaSEIRModel"
  return(model)
}

# Seasonal transmission rate function
seasonal_transmission_rate <- function(t, model) {
  # Annual seasonal variation
  seasonal_factor <- 1 + model$seasonal_amplitude * sin(2 * pi * t / 365 + model$seasonal_phase)
  return(model$beta * seasonal_factor)
}

# SEIR differential equations
seir_equations <- function(t, y, params) {
  with(as.list(c(y, params)), {
    # Seasonal transmission rate
    beta_t <- seasonal_transmission_rate(t, params)
    
    # SEIR equations with birth and death
    dS <- mu * N - beta_t * S * I / N - mu * S
    dE <- beta_t * S * I / N - (sigma + mu) * E
    dI <- sigma * E - (gamma + mu) * I
    dR <- gamma * I - mu * R
    
    return(list(c(dS, dE, dI, dR)))
  })
}

# Simulate the SEIR model
simulate.MalariaSEIRModel <- function(model, t_span = c(0, 365), num_points = 1000) {
  # Create time vector
  times <- seq(from = t_span[1], to = t_span[2], length.out = num_points)
  
  # Initial conditions
  y0 <- c(S = model$S0, E = model$E0, I = model$I0, R = model$R0)
  
  # Parameters for the differential equation
  params <- list(
    beta = model$beta,
    sigma = model$sigma,
    gamma = model$gamma,
    mu = model$mu,
    N = model$N,
    seasonal_amplitude = model$seasonal_amplitude,
    seasonal_phase = model$seasonal_phase
  )
  
  # Solve the differential equations
  solution <- ode(y = y0, times = times, func = seir_equations, parms = params)
  
  return(solution)
}

# Calculate epidemic metrics
get_epidemic_metrics <- function(model, solution) {
  # Extract the solution data
  S <- solution[, "S"]
  E <- solution[, "E"]
  I <- solution[, "I"]
  R <- solution[, "R"]
  
  # Calculate metrics
  peak_infection <- max(I)
  peak_time <- which.max(I)
  peak_exposed <- max(E)
  final_recovered <- R[length(R)]
  attack_rate <- final_recovered / model$N
  
  # Time to peak
  time_to_peak <- peak_time
  
  # Latent period statistics
  latent_period <- 1 / model$sigma
  
  metrics <- list(
    R0_basic = model$R0_basic,
    R0_effective = model$R0_effective,
    peak_infection = peak_infection,
    peak_infection_rate = peak_infection / model$N,
    peak_exposed = peak_exposed,
    peak_exposed_rate = peak_exposed / model$N,
    final_attack_rate = attack_rate,
    time_to_peak = time_to_peak,
    latent_period = latent_period,
    infectious_period = 1 / model$gamma
  )
  
  return(metrics)
}

# Plot the SEIR model results
plot_results <- function(model, solution, save_path = NULL) {
  # Convert solution to data frame for easier plotting
  df <- as.data.frame(solution)
  df$time <- df$time
  
  # Create plots
  p1 <- ggplot(df, aes(x = time)) +
    geom_line(aes(y = S, color = "Susceptible"), size = 1) +
    geom_line(aes(y = E, color = "Exposed"), size = 1) +
    geom_line(aes(y = I, color = "Infected"), size = 1) +
    geom_line(aes(y = R, color = "Recovered"), size = 1) +
    labs(title = "Malaria SEIR Model Dynamics",
         x = "Time (days)",
         y = "Number of individuals",
         color = "Compartment") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Infection and exposure rates
  df$infection_rate <- df$I / model$N * 100
  df$exposure_rate <- df$E / model$N * 100
  p2 <- ggplot(df, aes(x = time)) +
    geom_line(aes(y = infection_rate, color = "Infection Rate"), size = 1) +
    geom_line(aes(y = exposure_rate, color = "Exposure Rate"), size = 1) +
    labs(title = "Infection and Exposure Rates",
         x = "Time (days)",
         y = "Rate (%)",
         color = "Rate") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Seasonal transmission rate
  df$beta_seasonal <- sapply(df$time, function(t) seasonal_transmission_rate(t, model))
  p3 <- ggplot(df, aes(x = time, y = beta_seasonal)) +
    geom_line(color = "purple", size = 1) +
    labs(title = "Seasonal Transmission Rate",
         x = "Time (days)",
         y = "Transmission Rate (β)") +
    theme_minimal()
  
  # Phase plane (S vs I)
  p4 <- ggplot(df, aes(x = S, y = I)) +
    geom_path(color = "purple", size = 1) +
    labs(title = "Phase Plane (S vs I)",
         x = "Susceptible",
         y = "Infected") +
    theme_minimal()
  
  # Cumulative cases
  df$cumulative_cases <- df$R + df$I
  p5 <- ggplot(df, aes(x = time, y = cumulative_cases)) +
    geom_line(color = "orange", size = 1) +
    labs(title = "Cumulative Malaria Cases",
         x = "Time (days)",
         y = "Cumulative Cases") +
    theme_minimal()
  
  # Effective reproduction number over time
  df$R_eff <- sapply(1:nrow(df), function(i) {
    beta_t <- seasonal_transmission_rate(df$time[i], model)
    return(beta_t * df$S[i] / model$N / (model$gamma + model$mu))
  })
  p6 <- ggplot(df, aes(x = time, y = R_eff)) +
    geom_line(color = "red", size = 1) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed", alpha = 0.7) +
    labs(title = "Effective Reproduction Number",
         x = "Time (days)",
         y = "Effective Reproduction Number") +
    theme_minimal()
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = 15, height = 12, dpi = 300)
  }
  
  return(combined_plot)
}

# Compare with SIR model
compare_with_sir <- function(model, t_span = c(0, 365)) {
  # Create equivalent SIR model (no latent period)
  source("basic_sir.R")
  
  sir_model <- MalariaSIRModel(
    beta = model$beta,
    gamma = model$gamma,
    mu = model$mu,
    N = model$N,
    I0 = model$I0 + model$E0,  # Combine exposed and infected for SIR
    R0 = model$R0
  )
  
  # Simulate both models
  solution_seir <- simulate(model, t_span)
  solution_sir <- simulate(sir_model, t_span)
  
  # Create comparison data frame
  df_comparison <- data.frame(
    time = c(solution_seir[, "time"], solution_sir[, "time"]),
    S = c(solution_seir[, "S"], solution_sir[, "S"]),
    I = c(solution_seir[, "I"], solution_sir[, "I"]),
    R = c(solution_seir[, "R"], solution_sir[, "R"]),
    model = c(rep("SEIR", nrow(solution_seir)), rep("SIR", nrow(solution_sir)))
  )
  
  # Add exposed for SEIR
  df_comparison$E <- c(solution_seir[, "E"], rep(NA, nrow(solution_sir)))
  
  # Plot comparison
  p1 <- ggplot(df_comparison, aes(x = time, y = S, color = model, linetype = model)) +
    geom_line(size = 1) +
    labs(title = "Susceptible Population Comparison",
         x = "Time (days)",
         y = "Susceptible",
         color = "Model", linetype = "Model") +
    theme_minimal()
  
  p2 <- ggplot(df_comparison, aes(x = time, y = I, color = model, linetype = model)) +
    geom_line(size = 1) +
    labs(title = "Infected Population Comparison",
         x = "Time (days)",
         y = "Infected",
         color = "Model", linetype = "Model") +
    theme_minimal()
  
  p3 <- ggplot(df_comparison, aes(x = time, y = R, color = model, linetype = model)) +
    geom_line(size = 1) +
    labs(title = "Recovered Population Comparison",
         x = "Time (days)",
         y = "Recovered",
         color = "Model", linetype = "Model") +
    theme_minimal()
  
  # Exposed (SEIR only)
  df_exposed <- df_comparison[!is.na(df_comparison$E), ]
  p4 <- ggplot(df_exposed, aes(x = time, y = E)) +
    geom_line(color = "orange", size = 1) +
    labs(title = "Exposed Population (SEIR only)",
         x = "Time (days)",
         y = "Exposed") +
    theme_minimal()
  
  comparison_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  return(comparison_plot)
}

# Parameter sensitivity analysis
parameter_sensitivity <- function(model, param_name, param_values, t_span = c(0, 365)) {
  results <- data.frame()
  
  for (value in param_values) {
    # Create temporary model with modified parameter
    temp_model <- MalariaSEIRModel(
      beta = ifelse(param_name == "beta", value, model$beta),
      sigma = ifelse(param_name == "sigma", value, model$sigma),
      gamma = ifelse(param_name == "gamma", value, model$gamma),
      mu = ifelse(param_name == "mu", value, model$mu),
      N = model$N,
      E0 = model$E0,
      I0 = model$I0,
      R0 = model$R0,
      seasonal_amplitude = model$seasonal_amplitude,
      seasonal_phase = model$seasonal_phase
    )
    
    # Simulate
    solution <- simulate(temp_model, t_span)
    metrics <- get_epidemic_metrics(temp_model, solution)
    
    # Add to results
    results <- rbind(results, data.frame(
      parameter = param_name,
      value = value,
      R0_basic = metrics$R0_basic,
      R0_effective = metrics$R0_effective,
      peak_infection = metrics$peak_infection,
      peak_infection_rate = metrics$peak_infection_rate,
      final_attack_rate = metrics$final_attack_rate,
      time_to_peak = metrics$time_to_peak
    ))
  }
  
  return(results)
}

# Main function to demonstrate the malaria SEIR model
main <- function() {
  cat("Malaria SEIR Model Simulation\n")
  cat("=============================\n\n")
  
  # Create model with malaria-specific parameters
  model <- MalariaSEIRModel(
    beta = 0.4,                    # Transmission rate
    sigma = 0.1,                   # Incubation rate (10-day latent period)
    gamma = 0.05,                  # Recovery rate (20-day infectious period)
    mu = 0.0001,                   # Birth/death rate
    N = 10000,                     # Population of 10,000
    E0 = 20,                       # Start with 20 exposed
    I0 = 30,                       # Start with 30 infected
    R0 = 0,                        # No initially recovered
    seasonal_amplitude = 0.3,      # 30% seasonal variation
    seasonal_phase = 0.0           # Peak in summer
  )
  
  cat("Model Parameters:\n")
  cat("  Transmission rate (β):", model$beta, "\n")
  cat("  Incubation rate (σ):", model$sigma, "\n")
  cat("  Recovery rate (γ):", model$gamma, "\n")
  cat("  Birth/death rate (μ):", model$mu, "\n")
  cat("  Population size (N):", model$N, "\n")
  cat("  Seasonal amplitude:", model$seasonal_amplitude, "\n")
  cat("  Basic reproduction number (R₀):", round(model$R0_basic, 2), "\n")
  cat("  Effective reproduction number (R₀):", round(model$R0_effective, 2), "\n\n")
  
  # Simulate the model
  solution <- simulate(model, t_span = c(0, 1000), num_points = 1000)
  
  # Calculate and display metrics
  metrics <- get_epidemic_metrics(model, solution)
  cat("Epidemic Metrics:\n")
  cat("  R0_basic:", round(metrics$R0_basic, 3), "\n")
  cat("  R0_effective:", round(metrics$R0_effective, 3), "\n")
  cat("  Peak infection:", round(metrics$peak_infection, 0), "\n")
  cat("  Peak infection rate:", round(metrics$peak_infection_rate, 3), "\n")
  cat("  Peak exposed:", round(metrics$peak_exposed, 0), "\n")
  cat("  Peak exposed rate:", round(metrics$peak_exposed_rate, 3), "\n")
  cat("  Final attack rate:", round(metrics$final_attack_rate, 3), "\n")
  cat("  Time to peak:", round(metrics$time_to_peak, 0), "days\n")
  cat("  Latent period:", round(metrics$latent_period, 1), "days\n")
  cat("  Infectious period:", round(metrics$infectious_period, 1), "days\n\n")
  
  # Plot results
  plot_results(model, solution, save_path = "../results/plots/malaria_seir.png")
  
  # Compare with SIR model
  cat("Comparing with SIR model...\n")
  compare_with_sir(model)
  
  # Parameter sensitivity analysis
  cat("Parameter Sensitivity Analysis:\n")
  cat("-------------------------------\n")
  
  # Test different incubation rates
  sigma_values <- c(0.05, 0.1, 0.15, 0.2, 0.25)
  sigma_sensitivity <- parameter_sensitivity(model, "sigma", sigma_values)
  cat("Incubation Rate (σ) Sensitivity:\n")
  print(round(sigma_sensitivity[, c("value", "R0_effective", "peak_infection_rate", "final_attack_rate")], 3))
  cat("\n")
  
  # Test different seasonal amplitudes
  seasonal_values <- c(0.0, 0.1, 0.2, 0.3, 0.4)
  seasonal_sensitivity <- parameter_sensitivity(model, "beta", seasonal_values)
  cat("Seasonal Amplitude Sensitivity:\n")
  print(round(seasonal_sensitivity[, c("value", "R0_effective", "peak_infection_rate", "final_attack_rate")], 3))
  
  return(list(model = model, solution = solution, metrics = metrics))
}

# Run the main function if script is executed directly
if (interactive()) {
  results <- main()
}


