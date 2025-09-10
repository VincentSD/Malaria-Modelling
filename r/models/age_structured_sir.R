# Age-Structured SIR Model for Malaria Transmission
# ================================================
#
# This script implements an age-structured SIR model for malaria transmission,
# where the population is divided into different age groups with age-specific
# transmission rates, recovery rates, and contact patterns.
#
# Key features:
# - Multiple age groups with different epidemiological parameters
# - Age-specific contact matrices
# - Age-dependent susceptibility and infectiousness
# - Realistic population age distribution

# Load required libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

# Define the AgeStructuredMalariaSIR class (using R's S3 system)
AgeStructuredMalariaSIR <- function(age_groups = NULL,
                                  age_distribution = NULL,
                                  contact_matrix = NULL,
                                  beta_matrix = NULL,
                                  gamma = NULL,
                                  mu = NULL,
                                  N_total = 100000,
                                  I0 = NULL) {
  
  # Default age groups (WHO standard)
  if (is.null(age_groups)) {
    age_groups <- c('0-4', '5-14', '15-29', '30-44', '45-59', '60+')
  }
  
  n_age_groups <- length(age_groups)
  
  # Default age distribution (typical for malaria-endemic regions)
  if (is.null(age_distribution)) {
    age_distribution <- c(0.15, 0.25, 0.25, 0.15, 0.10, 0.10)
  }
  
  # Ensure age distribution sums to 1
  age_distribution <- age_distribution / sum(age_distribution)
  
  # Population in each age group
  N <- age_distribution * N_total
  
  # Default contact matrix (symmetric, higher within-group contact)
  if (is.null(contact_matrix)) {
    contact_matrix <- create_default_contact_matrix(n_age_groups)
  }
  
  # Default transmission rates (age-specific)
  if (is.null(beta_matrix)) {
    beta_matrix <- create_default_beta_matrix(n_age_groups)
  }
  
  # Default recovery rates (age-specific)
  if (is.null(gamma)) {
    gamma <- c(0.08, 0.06, 0.05, 0.05, 0.04, 0.03)  # Slower recovery in older ages
  }
  
  # Default birth/death rates (age-specific)
  if (is.null(mu)) {
    mu <- rep(0.0001, n_age_groups)
  }
  
  # Initial conditions
  if (is.null(I0)) {
    I0 <- c(5, 10, 15, 10, 5, 5)  # Initial infected by age
  }
  
  S0 <- N - I0  # Initial susceptible
  R0 <- rep(0, n_age_groups)  # Initial recovered
  
  # Calculate age-specific R0
  R0_age <- calculate_age_specific_r0(beta_matrix, contact_matrix, gamma, mu)
  
  # Create model object
  model <- list(
    age_groups = age_groups,
    n_age_groups = n_age_groups,
    age_distribution = age_distribution,
    N = N,
    N_total = N_total,
    contact_matrix = contact_matrix,
    beta_matrix = beta_matrix,
    gamma = gamma,
    mu = mu,
    S0 = S0,
    I0 = I0,
    R0 = R0,
    R0_age = R0_age
  )
  
  class(model) <- "AgeStructuredMalariaSIR"
  return(model)
}

# Helper function to create default contact matrix
create_default_contact_matrix <- function(n_age_groups) {
  # Higher contact within age groups, moderate between adjacent groups
  contact_matrix <- matrix(1, nrow = n_age_groups, ncol = n_age_groups)
  
  # Within-group contact is highest
  diag(contact_matrix) <- 2.0
  
  # Adjacent age groups have moderate contact
  for (i in 1:(n_age_groups - 1)) {
    contact_matrix[i, i + 1] <- 1.5
    contact_matrix[i + 1, i] <- 1.5
  }
  
  # Non-adjacent groups have lower contact
  for (i in 1:n_age_groups) {
    for (j in 1:n_age_groups) {
      if (abs(i - j) > 1) {
        contact_matrix[i, j] <- 0.5
      }
    }
  }
  
  return(contact_matrix)
}

# Helper function to create default beta matrix
create_default_beta_matrix <- function(n_age_groups) {
  # Base transmission rate
  beta_base <- 0.3
  
  # Age-specific susceptibility and infectiousness
  susceptibility <- c(1.2, 1.0, 0.8, 0.7, 0.6, 0.5)  # Higher in children
  infectiousness <- c(1.1, 1.0, 0.9, 0.8, 0.7, 0.6)  # Higher in children
  
  # Create transmission matrix
  beta_matrix <- outer(susceptibility, infectiousness) * beta_base
  
  return(beta_matrix)
}

# Helper function to calculate age-specific R0
calculate_age_specific_r0 <- function(beta_matrix, contact_matrix, gamma, mu) {
  n_age_groups <- length(gamma)
  R0_age <- numeric(n_age_groups)
  
  for (i in 1:n_age_groups) {
    # R0 for age group i
    R0_age[i] <- (beta_matrix[i, i] * contact_matrix[i, i]) / (gamma[i] + mu[i])
  }
  
  return(R0_age)
}

# Age-structured SIR differential equations
age_structured_sir_equations <- function(t, y, params) {
  with(as.list(c(y, params)), {
    # Extract compartments
    S <- y[1:n_age_groups]
    I <- y[(n_age_groups + 1):(2 * n_age_groups)]
    R <- y[(2 * n_age_groups + 1):(3 * n_age_groups)]
    
    # Initialize derivatives
    dSdt <- numeric(n_age_groups)
    dIdt <- numeric(n_age_groups)
    dRdt <- numeric(n_age_groups)
    
    # Calculate force of infection for each age group
    force_of_infection <- numeric(n_age_groups)
    
    for (i in 1:n_age_groups) {
      for (j in 1:n_age_groups) {
        # Force of infection from age group j to age group i
        force_of_infection[i] <- force_of_infection[i] + 
          (beta_matrix[i, j] * contact_matrix[i, j] * I[j] / N[j])
      }
    }
    
    # SIR equations for each age group
    for (i in 1:n_age_groups) {
      # Births enter susceptible compartment
      births <- mu[i] * N[i]
      
      # Susceptible dynamics
      dSdt[i] <- births - force_of_infection[i] * S[i] - mu[i] * S[i]
      
      # Infected dynamics
      dIdt[i] <- force_of_infection[i] * S[i] - (gamma[i] + mu[i]) * I[i]
      
      # Recovered dynamics
      dRdt[i] <- gamma[i] * I[i] - mu[i] * R[i]
    }
    
    return(list(c(dSdt, dIdt, dRdt)))
  })
}

# Simulate the age-structured SIR model
simulate.AgeStructuredMalariaSIR <- function(model, t_span = c(0, 365), num_points = 1000) {
  # Create time vector
  times <- seq(from = t_span[1], to = t_span[2], length.out = num_points)
  
  # Initial conditions
  y0 <- c(model$S0, model$I0, model$R0)
  
  # Parameters for the differential equation
  params <- list(
    n_age_groups = model$n_age_groups,
    beta_matrix = model$beta_matrix,
    contact_matrix = model$contact_matrix,
    gamma = model$gamma,
    mu = model$mu,
    N = model$N
  )
  
  # Solve the differential equations
  solution <- ode(y = y0, times = times, func = age_structured_sir_equations, parms = params)
  
  return(solution)
}

# Calculate age-specific epidemic metrics
get_age_specific_metrics <- function(model, solution) {
  # Extract the solution data
  S <- solution[, 2:(model$n_age_groups + 1)]
  I <- solution[, (model$n_age_groups + 2):(2 * model$n_age_groups + 1)]
  R <- solution[, (2 * model$n_age_groups + 2):(3 * model$n_age_groups + 1)]
  
  # Calculate metrics
  peak_infection <- apply(I, 2, max)
  peak_infection_rate <- peak_infection / model$N
  final_attack_rate <- R[nrow(R), ] / model$N
  time_to_peak <- apply(I, 2, which.max)
  
  metrics <- list(
    peak_infection = peak_infection,
    peak_infection_rate = peak_infection_rate,
    final_attack_rate = final_attack_rate,
    time_to_peak = time_to_peak,
    R0_age = model$R0_age
  )
  
  return(metrics)
}

# Plot the age-structured SIR model results
plot_age_structured_results <- function(model, solution, save_path = NULL) {
  # Extract compartments
  S <- solution[, 2:(model$n_age_groups + 1)]
  I <- solution[, (model$n_age_groups + 2):(2 * model$n_age_groups + 1)]
  R <- solution[, (2 * model$n_age_groups + 2):(3 * model$n_age_groups + 1)]
  
  # Create data frames for plotting
  time <- solution[, 1]
  
  # Total population dynamics
  total_S <- rowSums(S)
  total_I <- rowSums(I)
  total_R <- rowSums(R)
  
  df_total <- data.frame(
    time = time,
    S = total_S,
    I = total_I,
    R = total_R
  )
  
  # Age-specific infection rates
  infection_rates <- I / matrix(rep(model$N, each = nrow(I)), nrow = nrow(I)) * 100
  df_infection <- data.frame(
    time = rep(time, model$n_age_groups),
    infection_rate = as.vector(infection_rates),
    age_group = rep(model$age_groups, each = length(time))
  )
  
  # Final attack rates
  final_attack_rates <- R[nrow(R), ] / model$N * 100
  df_attack <- data.frame(
    age_group = model$age_groups,
    attack_rate = final_attack_rates
  )
  
  # Contact matrix for heatmap
  df_contact <- melt(model$contact_matrix)
  names(df_contact) <- c("Contactor", "Contactee", "Contact_Rate")
  df_contact$Contactor <- model$age_groups[df_contact$Contactor]
  df_contact$Contactee <- model$age_groups[df_contact$Contactee]
  
  # Transmission matrix for heatmap
  df_transmission <- melt(model$beta_matrix)
  names(df_transmission) <- c("Susceptible", "Infectious", "Transmission_Rate")
  df_transmission$Susceptible <- model$age_groups[df_transmission$Susceptible]
  df_transmission$Infectious <- model$age_groups[df_transmission$Infectious]
  
  # Age-specific R0
  df_r0 <- data.frame(
    age_group = model$age_groups,
    R0 = model$R0_age
  )
  
  # Create plots
  p1 <- ggplot(df_total, aes(x = time)) +
    geom_line(aes(y = S, color = "Susceptible"), size = 1) +
    geom_line(aes(y = I, color = "Infected"), size = 1) +
    geom_line(aes(y = R, color = "Recovered"), size = 1) +
    labs(title = "Total Population Dynamics",
         x = "Time (days)",
         y = "Number of individuals",
         color = "Compartment") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p2 <- ggplot(df_infection, aes(x = time, y = infection_rate, color = age_group)) +
    geom_line(size = 1) +
    labs(title = "Age-Specific Infection Rates",
         x = "Time (days)",
         y = "Infection Rate (%)",
         color = "Age Group") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p3 <- ggplot(df_attack, aes(x = age_group, y = attack_rate)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    geom_text(aes(label = sprintf("%.1f%%", attack_rate)), vjust = -0.5) +
    labs(title = "Age-Specific Final Attack Rates",
         x = "Age Group",
         y = "Final Attack Rate (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p4 <- ggplot(df_contact, aes(x = Contactee, y = Contactor, fill = Contact_Rate)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(title = "Contact Matrix",
         x = "Age Group (Contactee)",
         y = "Age Group (Contactor)",
         fill = "Contact Rate") +
    theme_minimal()
  
  p5 <- ggplot(df_transmission, aes(x = Infectious, y = Susceptible, fill = Transmission_Rate)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(title = "Transmission Rate Matrix",
         x = "Age Group (Infectious)",
         y = "Age Group (Susceptible)",
         fill = "Transmission Rate") +
    theme_minimal()
  
  p6 <- ggplot(df_r0, aes(x = age_group, y = R0)) +
    geom_bar(stat = "identity", fill = "red") +
    geom_text(aes(label = sprintf("%.2f", R0)), vjust = -0.5) +
    labs(title = "Age-Specific Basic Reproduction Numbers",
         x = "Age Group",
         y = "Basic Reproduction Number (R₀)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = 18, height = 12, dpi = 300)
  }
  
  return(combined_plot)
}

# Compare with basic SIR model
compare_with_basic_sir <- function(model, t_span = c(0, 365)) {
  # Create equivalent basic SIR model
  source("basic_sir.R")
  
  # Calculate population-weighted average parameters
  avg_beta <- mean(diag(model$beta_matrix))
  avg_gamma <- mean(model$gamma)
  avg_mu <- mean(model$mu)
  
  basic_model <- MalariaSIRModel(
    beta = avg_beta,
    gamma = avg_gamma,
    mu = avg_mu,
    N = model$N_total,
    I0 = sum(model$I0),
    R0 = 0
  )
  
  # Simulate both models
  solution_age <- simulate(model, t_span)
  solution_basic <- simulate(basic_model, t_span)
  
  # Extract total populations
  S_age <- rowSums(solution_age[, 2:(model$n_age_groups + 1)])
  I_age <- rowSums(solution_age[, (model$n_age_groups + 2):(2 * model$n_age_groups + 1)])
  R_age <- rowSums(solution_age[, (2 * model$n_age_groups + 2):(3 * model$n_age_groups + 1)])
  
  S_basic <- solution_basic[, "S"]
  I_basic <- solution_basic[, "I"]
  R_basic <- solution_basic[, "R"]
  
  time_age <- solution_age[, 1]
  time_basic <- solution_basic[, "time"]
  
  # Create comparison data frame
  df_comparison <- data.frame(
    time = c(time_age, time_basic),
    S = c(S_age, S_basic),
    I = c(I_age, I_basic),
    R = c(R_age, R_basic),
    model = c(rep("Age-structured", length(time_age)), rep("Basic SIR", length(time_basic)))
  )
  
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
  
  comparison_plot <- grid.arrange(p1, p2, p3, ncol = 1)
  
  return(comparison_plot)
}

# Main function to demonstrate the age-structured malaria SIR model
main <- function() {
  cat("Age-Structured Malaria SIR Model Simulation\n")
  cat("==========================================\n\n")
  
  # Create age-structured model
  model <- AgeStructuredMalariaSIR(
    N_total = 100000,
    I0 = c(10, 20, 30, 20, 10, 10)  # Initial infected by age
  )
  
  cat("Model Parameters:\n")
  cat("  Number of age groups:", model$n_age_groups, "\n")
  cat("  Age groups:", paste(model$age_groups, collapse = ", "), "\n")
  cat("  Total population:", format(model$N_total, big.mark = ","), "\n")
  cat("  Age distribution:", paste(round(model$age_distribution, 3), collapse = ", "), "\n\n")
  
  cat("Age-specific R₀ values:\n")
  for (i in 1:model$n_age_groups) {
    cat("  Age", model$age_groups[i], ":", round(model$R0_age[i], 2), "\n")
  }
  cat("\n")
  
  # Simulate the model
  solution <- simulate(model, t_span = c(0, 1000), num_points = 1000)
  
  # Calculate and display metrics
  metrics <- get_age_specific_metrics(model, solution)
  
  cat("Age-specific epidemic metrics:\n")
  cat("------------------------------\n")
  for (i in 1:model$n_age_groups) {
    cat("Age", model$age_groups[i], ":\n")
    cat("  Peak infection rate:", round(metrics$peak_infection_rate[i], 3), "\n")
    cat("  Final attack rate:", round(metrics$final_attack_rate[i], 3), "\n")
    cat("  Time to peak:", metrics$time_to_peak[i], "days\n\n")
  }
  
  # Plot results
  plot_age_structured_results(model, solution, save_path = "../results/plots/malaria_age_structured.png")
  
  # Compare with basic SIR
  cat("Comparing with basic SIR model...\n")
  compare_with_basic_sir(model)
  
  return(list(model = model, solution = solution, metrics = metrics))
}

# Run the main function if script is executed directly
if (interactive()) {
  results <- main()
}


