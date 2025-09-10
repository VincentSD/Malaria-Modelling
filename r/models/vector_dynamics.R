# Vector Dynamics Model for Malaria Transmission
# =============================================
#
# This script implements a comprehensive malaria transmission model that includes
# both human and mosquito (vector) population dynamics. This is essential for
# malaria modeling as malaria is a vector-borne disease.
#
# Key features:
# - Human population: Susceptible-Exposed-Infected-Recovered (SEIR)
# - Mosquito population: Susceptible-Exposed-Infected (SEI)
# - Human-mosquito transmission cycles
# - Seasonal variations in mosquito populations
# - Vector control interventions

# Load required libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Define the MalariaVectorModel class (using R's S3 system)
MalariaVectorModel <- function(beta_h = 0.3,      # Human transmission rate
                              sigma_h = 0.1,     # Human incubation rate
                              gamma_h = 0.05,    # Human recovery rate
                              mu_h = 0.0001,     # Human birth/death rate
                              N_h = 10000,       # Human population
                              beta_v = 0.2,      # Vector transmission rate
                              sigma_v = 0.2,     # Vector incubation rate
                              mu_v = 0.1,        # Vector death rate
                              N_v = 50000,       # Vector population
                              a = 0.3,           # Biting rate
                              b = 0.5,           # Human to vector transmission probability
                              c = 0.5,           # Vector to human transmission probability
                              seasonal_amplitude = 0.5,  # Seasonal variation amplitude
                              seasonal_phase = 0.0,      # Seasonal phase offset
                              E_h0 = 10,         # Initial exposed humans
                              I_h0 = 20,         # Initial infected humans
                              R_h0 = 0,          # Initial recovered humans
                              E_v0 = 100,        # Initial exposed vectors
                              I_v0 = 200) {      # Initial infected vectors
  
  # Calculate initial susceptible populations
  S_h0 <- N_h - E_h0 - I_h0 - R_h0
  S_v0 <- N_v - E_v0 - I_v0
  
  # Calculate basic reproduction number
  R0 <- (a^2 * b * c * N_v) / (mu_v * (gamma_h + mu_h) * N_h)
  
  # Create model object
  model <- list(
    beta_h = beta_h,
    sigma_h = sigma_h,
    gamma_h = gamma_h,
    mu_h = mu_h,
    N_h = N_h,
    beta_v = beta_v,
    sigma_v = sigma_v,
    mu_v = mu_v,
    N_v = N_v,
    a = a,
    b = b,
    c = c,
    seasonal_amplitude = seasonal_amplitude,
    seasonal_phase = seasonal_phase,
    S_h0 = S_h0,
    E_h0 = E_h0,
    I_h0 = I_h0,
    R_h0 = R_h0,
    S_v0 = S_v0,
    E_v0 = E_v0,
    I_v0 = I_v0,
    R0 = R0
  )
  
  class(model) <- "MalariaVectorModel"
  return(model)
}

# Seasonal vector population function
seasonal_vector_population <- function(t, model) {
  # Annual seasonal variation
  seasonal_factor <- 1 + model$seasonal_amplitude * sin(2 * pi * t / 365 + model$seasonal_phase)
  return(seasonal_factor)
}

# Vector dynamics differential equations
vector_dynamics_equations <- function(t, y, params) {
  with(as.list(c(y, params)), {
    # Extract compartments
    S_h <- y[1]
    E_h <- y[2]
    I_h <- y[3]
    R_h <- y[4]
    S_v <- y[5]
    E_v <- y[6]
    I_v <- y[7]
    
    # Seasonal vector population
    N_v_seasonal <- N_v * seasonal_vector_population(t, params)
    
    # Force of infection from vectors to humans
    lambda_h <- a * c * I_v / N_v_seasonal
    
    # Force of infection from humans to vectors
    lambda_v <- a * b * I_h / N_h
    
    # Human dynamics
    dS_h <- mu_h * N_h - lambda_h * S_h - mu_h * S_h
    dE_h <- lambda_h * S_h - (sigma_h + mu_h) * E_h
    dI_h <- sigma_h * E_h - (gamma_h + mu_h) * I_h
    dR_h <- gamma_h * I_h - mu_h * R_h
    
    # Vector dynamics
    dS_v <- mu_v * N_v_seasonal - lambda_v * S_v - mu_v * S_v
    dE_v <- lambda_v * S_v - (sigma_v + mu_v) * E_v
    dI_v <- sigma_v * E_v - mu_v * I_v
    
    return(list(c(dS_h, dE_h, dI_h, dR_h, dS_v, dE_v, dI_v)))
  })
}

# Simulate the vector dynamics model
simulate.MalariaVectorModel <- function(model, t_span = c(0, 365), num_points = 1000) {
  # Create time vector
  times <- seq(from = t_span[1], to = t_span[2], length.out = num_points)
  
  # Initial conditions
  y0 <- c(S_h = model$S_h0, E_h = model$E_h0, I_h = model$I_h0, R_h = model$R_h0,
          S_v = model$S_v0, E_v = model$E_v0, I_v = model$I_v0)
  
  # Parameters for the differential equation
  params <- list(
    beta_h = model$beta_h,
    sigma_h = model$sigma_h,
    gamma_h = model$gamma_h,
    mu_h = model$mu_h,
    N_h = model$N_h,
    beta_v = model$beta_v,
    sigma_v = model$sigma_v,
    mu_v = model$mu_v,
    N_v = model$N_v,
    a = model$a,
    b = model$b,
    c = model$c,
    seasonal_amplitude = model$seasonal_amplitude,
    seasonal_phase = model$seasonal_phase
  )
  
  # Solve the differential equations
  solution <- ode(y = y0, times = times, func = vector_dynamics_equations, parms = params)
  
  return(solution)
}

# Calculate epidemic metrics
get_epidemic_metrics <- function(model, solution) {
  # Extract the solution data
  S_h <- solution[, "S_h"]
  E_h <- solution[, "E_h"]
  I_h <- solution[, "I_h"]
  R_h <- solution[, "R_h"]
  S_v <- solution[, "S_v"]
  E_v <- solution[, "E_v"]
  I_v <- solution[, "I_v"]
  
  # Human metrics
  peak_human_infection <- max(I_h)
  peak_human_infection_rate <- peak_human_infection / model$N_h
  final_human_attack_rate <- R_h[length(R_h)] / model$N_h
  
  # Vector metrics
  peak_vector_infection <- max(I_v)
  peak_vector_infection_rate <- peak_vector_infection / model$N_v
  
  # Time to peak
  time_to_peak_human <- which.max(I_h)
  time_to_peak_vector <- which.max(I_v)
  
  metrics <- list(
    R0_basic = model$R0,
    peak_human_infection = peak_human_infection,
    peak_human_infection_rate = peak_human_infection_rate,
    final_human_attack_rate = final_human_attack_rate,
    peak_vector_infection = peak_vector_infection,
    peak_vector_infection_rate = peak_vector_infection_rate,
    time_to_peak_human = time_to_peak_human,
    time_to_peak_vector = time_to_peak_vector,
    vector_human_ratio = model$N_v / model$N_h
  )
  
  return(metrics)
}

# Plot the vector dynamics model results
plot_results <- function(model, solution, save_path = NULL) {
  # Convert solution to data frame for easier plotting
  df <- as.data.frame(solution)
  df$time <- df$time
  
  # Create plots
  p1 <- ggplot(df, aes(x = time)) +
    geom_line(aes(y = S_h, color = "Susceptible"), size = 1) +
    geom_line(aes(y = E_h, color = "Exposed"), size = 1) +
    geom_line(aes(y = I_h, color = "Infected"), size = 1) +
    geom_line(aes(y = R_h, color = "Recovered"), size = 1) +
    labs(title = "Human Population Dynamics",
         x = "Time (days)",
         y = "Number of humans",
         color = "Compartment") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p2 <- ggplot(df, aes(x = time)) +
    geom_line(aes(y = S_v, color = "Susceptible"), size = 1) +
    geom_line(aes(y = E_v, color = "Exposed"), size = 1) +
    geom_line(aes(y = I_v, color = "Infected"), size = 1) +
    labs(title = "Vector Population Dynamics",
         x = "Time (days)",
         y = "Number of vectors",
         color = "Compartment") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Infection rates comparison
  df$human_infection_rate <- df$I_h / model$N_h * 100
  df$vector_infection_rate <- df$I_v / model$N_v * 100
  p3 <- ggplot(df, aes(x = time)) +
    geom_line(aes(y = human_infection_rate, color = "Human"), size = 1) +
    geom_line(aes(y = vector_infection_rate, color = "Vector"), size = 1) +
    labs(title = "Infection Rates Comparison",
         x = "Time (days)",
         y = "Infection Rate (%)",
         color = "Population") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Seasonal vector population
  df$N_v_seasonal <- sapply(df$time, function(t) model$N_v * seasonal_vector_population(t, model))
  p4 <- ggplot(df, aes(x = time, y = N_v_seasonal)) +
    geom_line(color = "green", size = 1) +
    labs(title = "Seasonal Vector Population",
         x = "Time (days)",
         y = "Vector Population") +
    theme_minimal()
  
  # Cumulative human cases
  df$cumulative_cases <- df$R_h + df$I_h
  p5 <- ggplot(df, aes(x = time, y = cumulative_cases)) +
    geom_line(color = "orange", size = 1) +
    labs(title = "Cumulative Human Malaria Cases",
         x = "Time (days)",
         y = "Cumulative Cases") +
    theme_minimal()
  
  # Vector-human transmission cycle
  p6 <- ggplot(df, aes(x = I_h, y = I_v)) +
    geom_path(color = "purple", size = 1) +
    labs(title = "Vector-Human Transmission Cycle",
         x = "Infected Humans",
         y = "Infected Vectors") +
    theme_minimal()
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = 15, height = 12, dpi = 300)
  }
  
  return(combined_plot)
}

# Intervention analysis
intervention_analysis <- function(model, intervention_type, intervention_strength, 
                                intervention_start, intervention_duration, 
                                t_span = c(0, 1000)) {
  # Create modified model with intervention
  if (intervention_type == "insecticide") {
    # Reduce vector population
    modified_model <- MalariaVectorModel(
      beta_h = model$beta_h, sigma_h = model$sigma_h, gamma_h = model$gamma_h, 
      mu_h = model$mu_h, N_h = model$N_h,
      beta_v = model$beta_v, sigma_v = model$sigma_v, 
      mu_v = model$mu_v * (1 + intervention_strength), N_v = model$N_v,
      a = model$a, b = model$b, c = model$c,
      seasonal_amplitude = model$seasonal_amplitude, seasonal_phase = model$seasonal_phase,
      E_h0 = model$E_h0, I_h0 = model$I_h0, R_h0 = model$R_h0, 
      E_v0 = model$E_v0, I_v0 = model$I_v0
    )
  } else if (intervention_type == "bed_nets") {
    # Reduce biting rate
    modified_model <- MalariaVectorModel(
      beta_h = model$beta_h, sigma_h = model$sigma_h, gamma_h = model$gamma_h, 
      mu_h = model$mu_h, N_h = model$N_h,
      beta_v = model$beta_v, sigma_v = model$sigma_v, mu_v = model$mu_v, N_v = model$N_v,
      a = model$a * (1 - intervention_strength), b = model$b, c = model$c,
      seasonal_amplitude = model$seasonal_amplitude, seasonal_phase = model$seasonal_phase,
      E_h0 = model$E_h0, I_h0 = model$I_h0, R_h0 = model$R_h0, 
      E_v0 = model$E_v0, I_v0 = model$I_v0
    )
  } else if (intervention_type == "larvicide") {
    # Reduce vector birth rate (affects seasonal population)
    modified_model <- MalariaVectorModel(
      beta_h = model$beta_h, sigma_h = model$sigma_h, gamma_h = model$gamma_h, 
      mu_h = model$mu_h, N_h = model$N_h,
      beta_v = model$beta_v, sigma_v = model$sigma_v, mu_v = model$mu_v, N_v = model$N_v,
      a = model$a, b = model$b, c = model$c,
      seasonal_amplitude = model$seasonal_amplitude * (1 - intervention_strength), 
      seasonal_phase = model$seasonal_phase,
      E_h0 = model$E_h0, I_h0 = model$I_h0, R_h0 = model$R_h0, 
      E_v0 = model$E_v0, I_v0 = model$I_v0
    )
  } else {
    stop("Invalid intervention type")
  }
  
  # Simulate both models
  solution_baseline <- simulate(model, t_span)
  solution_intervention <- simulate(modified_model, t_span)
  
  # Calculate metrics
  metrics_baseline <- get_epidemic_metrics(model, solution_baseline)
  metrics_intervention <- get_epidemic_metrics(modified_model, solution_intervention)
  
  # Calculate intervention effectiveness
  effectiveness <- list(
    intervention_type = intervention_type,
    intervention_strength = intervention_strength,
    baseline_metrics = metrics_baseline,
    intervention_metrics = metrics_intervention,
    reduction_in_peak_infection = (metrics_baseline$peak_human_infection_rate - 
                                  metrics_intervention$peak_human_infection_rate) / 
                                 metrics_baseline$peak_human_infection_rate,
    reduction_in_attack_rate = (metrics_baseline$final_human_attack_rate - 
                               metrics_intervention$final_human_attack_rate) / 
                              metrics_baseline$final_human_attack_rate
  )
  
  return(effectiveness)
}

# Main function to demonstrate the malaria vector dynamics model
main <- function() {
  cat("Malaria Vector Dynamics Model Simulation\n")
  cat("=======================================\n\n")
  
  # Create model with realistic parameters
  model <- MalariaVectorModel(
    # Human parameters
    beta_h = 0.3, sigma_h = 0.1, gamma_h = 0.05, mu_h = 0.0001, N_h = 10000,
    # Vector parameters
    beta_v = 0.2, sigma_v = 0.2, mu_v = 0.1, N_v = 50000,
    # Biting parameters
    a = 0.3, b = 0.5, c = 0.5,
    # Seasonal parameters
    seasonal_amplitude = 0.5, seasonal_phase = 0.0,
    # Initial conditions
    E_h0 = 20, I_h0 = 30, R_h0 = 0, E_v0 = 500, I_v0 = 1000
  )
  
  cat("Model Parameters:\n")
  cat("  Human population:", format(model$N_h, big.mark = ","), "\n")
  cat("  Vector population:", format(model$N_v, big.mark = ","), "\n")
  cat("  Vector-human ratio:", round(model$N_v/model$N_h, 1), "\n")
  cat("  Biting rate (a):", model$a, "\n")
  cat("  Human to vector transmission (b):", model$b, "\n")
  cat("  Vector to human transmission (c):", model$c, "\n")
  cat("  Basic reproduction number (R₀):", round(model$R0, 2), "\n\n")
  
  # Simulate the model
  solution <- simulate(model, t_span = c(0, 1000), num_points = 1000)
  
  # Calculate and display metrics
  metrics <- get_epidemic_metrics(model, solution)
  cat("Epidemic Metrics:\n")
  cat("  R0_basic:", round(metrics$R0_basic, 3), "\n")
  cat("  Peak human infection:", round(metrics$peak_human_infection, 0), "\n")
  cat("  Peak human infection rate:", round(metrics$peak_human_infection_rate, 3), "\n")
  cat("  Final human attack rate:", round(metrics$final_human_attack_rate, 3), "\n")
  cat("  Peak vector infection:", round(metrics$peak_vector_infection, 0), "\n")
  cat("  Peak vector infection rate:", round(metrics$peak_vector_infection_rate, 3), "\n")
  cat("  Time to peak (human):", round(metrics$time_to_peak_human, 0), "days\n")
  cat("  Time to peak (vector):", round(metrics$time_to_peak_vector, 0), "days\n")
  cat("  Vector-human ratio:", round(metrics$vector_human_ratio, 1), "\n\n")
  
  # Plot results
  plot_results(model, solution, save_path = "../results/plots/malaria_vector_dynamics.png")
  
  # Intervention analysis
  cat("Intervention Analysis:\n")
  cat("---------------------\n")
  
  interventions <- list(
    list("insecticide", 0.5, 200, 100),
    list("bed_nets", 0.7, 200, 100),
    list("larvicide", 0.6, 200, 100)
  )
  
  for (intervention in interventions) {
    intervention_type <- intervention[[1]]
    strength <- intervention[[2]]
    start <- intervention[[3]]
    duration <- intervention[[4]]
    
    effectiveness <- intervention_analysis(
      model, intervention_type, strength, start, duration
    )
    
    cat(paste0(toupper(intervention_type), " (strength: ", strength, "):\n"))
    cat("  Reduction in peak infection:", round(effectiveness$reduction_in_peak_infection * 100, 1), "%\n")
    cat("  Reduction in attack rate:", round(effectiveness$reduction_in_attack_rate * 100, 1), "%\n\n")
  }
  
  return(list(model = model, solution = solution, metrics = metrics))
}

# Run the main function if script is executed directly
if (interactive()) {
  results <- main()
}


