# Malaria Model Simulation Example (R Version)
# ============================================
#
# This script demonstrates the complete malaria model simulation project
# using R implementations of all models.

# Load required libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Source all model files
source("models/basic_sir.R")
source("models/seir_model.R")
source("models/age_structured_sir.R")
source("models/vector_dynamics.R")
source("analysis/visualization.R")
source("analysis/parameter_analysis.R")

# Main function to run complete example
main <- function() {
  cat("Malaria Model Simulation Project (R Version)\n")
  cat("============================================\n\n")
  
  # Create results directory
  if (!dir.exists("results")) {
    dir.create("results", recursive = TRUE)
  }
  if (!dir.exists("results/plots")) {
    dir.create("results/plots", recursive = TRUE)
  }
  if (!dir.exists("results/analysis")) {
    dir.create("results/analysis", recursive = TRUE)
  }
  
  # 1. Basic SIR Model
  cat("1. Running Basic SIR Model...\n")
  cat("-----------------------------\n")
  
  sir_model <- MalariaSIRModel(
    beta = 0.4,      # Transmission rate
    gamma = 0.05,    # Recovery rate
    mu = 0.0001,     # Birth/death rate
    N = 10000,       # Population size
    I0 = 50,         # Initial infected
    R0 = 0           # Initial recovered
  )
  
  sir_solution <- simulate(sir_model, t_span = c(0, 365), num_points = 1000)
  sir_metrics <- get_epidemic_metrics(sir_model, sir_solution)
  
  cat("Basic SIR Model Results:\n")
  cat("  R0:", round(sir_metrics$R0, 3), "\n")
  cat("  Peak infection:", round(sir_metrics$peak_infection, 0), "\n")
  cat("  Peak infection rate:", round(sir_metrics$peak_infection_rate, 3), "\n")
  cat("  Final attack rate:", round(sir_metrics$final_attack_rate, 3), "\n")
  cat("  Time to peak:", round(sir_metrics$time_to_peak, 0), "days\n\n")
  
  # Plot basic SIR results
  plot_results(sir_model, sir_solution, save_path = "results/plots/basic_sir.png")
  
  # 2. SEIR Model
  cat("2. Running SEIR Model...\n")
  cat("------------------------\n")
  
  seir_model <- MalariaSEIRModel(
    beta = 0.4,                    # Transmission rate
    sigma = 0.1,                   # Incubation rate (10-day latent period)
    gamma = 0.05,                  # Recovery rate (20-day infectious period)
    mu = 0.0001,                   # Birth/death rate
    N = 10000,                     # Population size
    E0 = 20,                       # Initial exposed
    I0 = 30,                       # Initial infected
    R0 = 0,                        # Initial recovered
    seasonal_amplitude = 0.3,      # 30% seasonal variation
    seasonal_phase = 0.0           # Peak in summer
  )
  
  seir_solution <- simulate(seir_model, t_span = c(0, 365), num_points = 1000)
  seir_metrics <- get_epidemic_metrics(seir_model, seir_solution)
  
  cat("SEIR Model Results:\n")
  cat("  R0_basic:", round(seir_metrics$R0_basic, 3), "\n")
  cat("  R0_effective:", round(seir_metrics$R0_effective, 3), "\n")
  cat("  Peak infection:", round(seir_metrics$peak_infection, 0), "\n")
  cat("  Peak infection rate:", round(seir_metrics$peak_infection_rate, 3), "\n")
  cat("  Peak exposed:", round(seir_metrics$peak_exposed, 0), "\n")
  cat("  Final attack rate:", round(seir_metrics$final_attack_rate, 3), "\n")
  cat("  Time to peak:", round(seir_metrics$time_to_peak, 0), "days\n")
  cat("  Latent period:", round(seir_metrics$latent_period, 1), "days\n\n")
  
  # Plot SEIR results
  plot_results(seir_model, seir_solution, save_path = "results/plots/seir_model.png")
  
  # 3. Age-Structured SIR Model
  cat("3. Running Age-Structured SIR Model...\n")
  cat("--------------------------------------\n")
  
  age_sir_model <- AgeStructuredMalariaSIR(
    N_total = 10000,               # Total population
    age_groups = c("0-4", "5-14", "15-24", "25-44", "45-64", "65+"),
    N = c(1500, 2000, 2000, 2500, 1500, 500),  # Age group populations
    I0 = c(5, 10, 15, 10, 5, 5),  # Initial infected by age group
    R0 = rep(0, 6),                # Initial recovered by age group
    beta = 0.4,                    # Base transmission rate
    gamma = 0.05,                  # Recovery rate
    mu = 0.0001,                   # Birth/death rate
    contact_matrix = NULL          # Will be generated automatically
  )
  
  age_sir_solution <- simulate(age_sir_model, t_span = c(0, 365), num_points = 1000)
  age_sir_metrics <- get_epidemic_metrics(age_sir_model, age_sir_solution)
  
  cat("Age-Structured SIR Model Results:\n")
  cat("  R0_basic:", round(age_sir_metrics$R0_basic, 3), "\n")
  cat("  Peak infection:", round(age_sir_metrics$peak_infection, 0), "\n")
  cat("  Peak infection rate:", round(age_sir_metrics$peak_infection_rate, 3), "\n")
  cat("  Final attack rate:", round(age_sir_metrics$final_attack_rate, 3), "\n")
  cat("  Time to peak:", round(age_sir_metrics$time_to_peak, 0), "days\n")
  cat("  Age groups:", length(age_sir_model$age_groups), "\n\n")
  
  # Plot age-structured results
  plot_age_structured_results(age_sir_model, age_sir_solution, 
                            save_path = "results/plots/age_structured_sir.png")
  
  # 4. Vector Dynamics Model
  cat("4. Running Vector Dynamics Model...\n")
  cat("-----------------------------------\n")
  
  vector_model <- MalariaVectorModel(
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
  
  vector_solution <- simulate(vector_model, t_span = c(0, 365), num_points = 1000)
  vector_metrics <- get_epidemic_metrics(vector_model, vector_solution)
  
  cat("Vector Dynamics Model Results:\n")
  cat("  R0_basic:", round(vector_metrics$R0_basic, 3), "\n")
  cat("  Peak human infection:", round(vector_metrics$peak_human_infection, 0), "\n")
  cat("  Peak human infection rate:", round(vector_metrics$peak_human_infection_rate, 3), "\n")
  cat("  Final human attack rate:", round(vector_metrics$final_human_attack_rate, 3), "\n")
  cat("  Peak vector infection:", round(vector_metrics$peak_vector_infection, 0), "\n")
  cat("  Peak vector infection rate:", round(vector_metrics$peak_vector_infection_rate, 3), "\n")
  cat("  Vector-human ratio:", round(vector_metrics$vector_human_ratio, 1), "\n\n")
  
  # Plot vector dynamics results
  plot_vector_dynamics_results(vector_model, vector_solution, 
                             save_path = "results/plots/vector_dynamics.png")
  
  # 5. Model Comparison
  cat("5. Creating Model Comparison...\n")
  cat("-------------------------------\n")
  
  models <- list(
    "Basic SIR" = sir_model,
    "SEIR" = seir_model,
    "Age-Structured SIR" = age_sir_model,
    "Vector Dynamics" = vector_model
  )
  
  plot_model_comparison(models, t_span = c(0, 365), 
                       save_path = "results/plots/model_comparison.png")
  
  # 6. Parameter Sensitivity Analysis
  cat("6. Running Parameter Sensitivity Analysis...\n")
  cat("---------------------------------------------\n")
  
  # Define parameter ranges
  param_ranges <- list(
    beta = c(0.1, 0.8),
    gamma = c(0.01, 0.2),
    mu = c(0.00005, 0.0005)
  )
  
  # Run sensitivity analysis on basic SIR model
  sensitivity_results <- parameter_sensitivity_analysis(sir_model, param_ranges)
  
  # Plot sensitivity results
  plot_sensitivity_results(sensitivity_results, 
                          save_path = "results/analysis/sensitivity_analysis.png")
  
  # Save sensitivity results
  write.csv(sensitivity_results, "results/analysis/sensitivity_results.csv", row.names = FALSE)
  
  # 7. Intervention Analysis
  cat("7. Running Intervention Analysis...\n")
  cat("-----------------------------------\n")
  
  # Test different interventions on vector model
  interventions <- list(
    list("insecticide", 0.5, 200, 100),
    list("bed_nets", 0.7, 200, 100),
    list("larvicide", 0.6, 200, 100)
  )
  
  intervention_results <- list()
  
  for (intervention in interventions) {
    intervention_type <- intervention[[1]]
    strength <- intervention[[2]]
    start <- intervention[[3]]
    duration <- intervention[[4]]
    
    effectiveness <- intervention_analysis(
      vector_model, intervention_type, strength, start, duration
    )
    
    intervention_results[[intervention_type]] <- effectiveness
    
    cat(paste0(toupper(intervention_type), " (strength: ", strength, "):\n"))
    cat("  Reduction in peak infection:", round(effectiveness$reduction_in_peak_infection * 100, 1), "%\n")
    cat("  Reduction in attack rate:", round(effectiveness$reduction_in_attack_rate * 100, 1), "%\n\n")
  }
  
  # Plot intervention effectiveness
  plot_intervention_effectiveness(intervention_results, 
                                 save_path = "results/analysis/intervention_effectiveness.png")
  
  # 8. Summary Report
  cat("8. Generating Summary Report...\n")
  cat("-------------------------------\n")
  
  # Create summary data frame
  summary_data <- data.frame(
    Model = c("Basic SIR", "SEIR", "Age-Structured SIR", "Vector Dynamics"),
    R0 = c(sir_metrics$R0, seir_metrics$R0_effective, age_sir_metrics$R0_basic, vector_metrics$R0_basic),
    Peak_Infection_Rate = c(sir_metrics$peak_infection_rate, seir_metrics$peak_infection_rate, 
                           age_sir_metrics$peak_infection_rate, vector_metrics$peak_human_infection_rate),
    Final_Attack_Rate = c(sir_metrics$final_attack_rate, seir_metrics$final_attack_rate, 
                         age_sir_metrics$final_attack_rate, vector_metrics$final_human_attack_rate),
    Time_to_Peak = c(sir_metrics$time_to_peak, seir_metrics$time_to_peak, 
                    age_sir_metrics$time_to_peak, vector_metrics$time_to_peak_human)
  )
  
  # Save summary
  write.csv(summary_data, "results/analysis/model_summary.csv", row.names = FALSE)
  
  # Print summary
  cat("Model Summary:\n")
  print(round(summary_data, 3))
  
  cat("\nAnalysis complete! All results saved to the results/ directory.\n")
  cat("Generated files:\n")
  cat("  - results/plots/: Model visualizations\n")
  cat("  - results/analysis/: Analysis results and summaries\n")
  
  return(list(
    models = models,
    solutions = list(sir_solution, seir_solution, age_sir_solution, vector_solution),
    metrics = list(sir_metrics, seir_metrics, age_sir_metrics, vector_metrics),
    sensitivity_results = sensitivity_results,
    intervention_results = intervention_results,
    summary_data = summary_data
  ))
}

# Run the main function if script is executed directly
if (interactive()) {
  results <- main()
}


