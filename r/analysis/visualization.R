# Advanced Visualization Tools for Malaria Models (R Version)
# ==========================================================
#
# This script provides comprehensive visualization tools for malaria transmission models,
# including comparison plots, parameter sensitivity analysis, and interactive visualizations.

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(viridis)

# Model comparison visualization
plot_model_comparison <- function(models, t_span = c(0, 365), save_path = NULL) {
  # Prepare data for plotting
  plot_data <- data.frame()
  
  for (model_name in names(models)) {
    model <- models[[model_name]]
    
    # Simulate model
    if (inherits(model, "MalariaSIRModel")) {
      solution <- simulate(model, t_span)
      S <- solution[, "S"]
      I <- solution[, "I"]
      R <- solution[, "R"]
      time <- solution[, "time"]
      N <- model$N
    } else if (inherits(model, "MalariaSEIRModel")) {
      solution <- simulate(model, t_span)
      S <- solution[, "S"]
      I <- solution[, "I"]
      R <- solution[, "R"]
      time <- solution[, "time"]
      N <- model$N
    } else if (inherits(model, "AgeStructuredMalariaSIR")) {
      solution <- simulate(model, t_span)
      S <- rowSums(solution[, 2:(model$n_age_groups + 1)])
      I <- rowSums(solution[, (model$n_age_groups + 2):(2 * model$n_age_groups + 1)])
      R <- rowSums(solution[, (2 * model$n_age_groups + 2):(3 * model$n_age_groups + 1)])
      time <- solution[, 1]
      N <- model$N_total
    } else if (inherits(model, "MalariaVectorModel")) {
      solution <- simulate(model, t_span)
      S <- solution[, "S_h"]
      I <- solution[, "I_h"]
      R <- solution[, "R_h"]
      time <- solution[, "time"]
      N <- model$N_h
    }
    
    # Calculate infection rate
    infection_rate <- I / N * 100
    
    # Add to plot data
    model_data <- data.frame(
      time = time,
      S = S,
      I = I,
      R = R,
      infection_rate = infection_rate,
      model = model_name
    )
    
    plot_data <- rbind(plot_data, model_data)
  }
  
  # Create comparison plots
  p1 <- ggplot(plot_data, aes(x = time, y = S, color = model)) +
    geom_line(size = 1) +
    labs(title = "Susceptible Population Comparison",
         x = "Time (days)",
         y = "Number of individuals",
         color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p2 <- ggplot(plot_data, aes(x = time, y = I, color = model)) +
    geom_line(size = 1) +
    labs(title = "Infected Population Comparison",
         x = "Time (days)",
         y = "Number of individuals",
         color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p3 <- ggplot(plot_data, aes(x = time, y = R, color = model)) +
    geom_line(size = 1) +
    labs(title = "Recovered Population Comparison",
         x = "Time (days)",
         y = "Number of individuals",
         color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p4 <- ggplot(plot_data, aes(x = time, y = infection_rate, color = model)) +
    geom_line(size = 1) +
    labs(title = "Infection Rate Comparison",
         x = "Time (days)",
         y = "Infection Rate (%)",
         color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = 15, height = 10, dpi = 300)
  }
  
  return(combined_plot)
}

# Parameter sensitivity visualization
plot_parameter_sensitivity <- function(sensitivity_results, save_path = NULL) {
  # Get unique parameters and metrics
  parameters <- unique(sensitivity_results$parameter)
  metrics <- c("R0_basic", "peak_infection_rate", "final_attack_rate", "time_to_peak")
  metric_titles <- c("Basic R₀", "Peak Infection Rate", "Final Attack Rate", "Time to Peak")
  
  # Create subplots
  plots <- list()
  
  for (i in seq_along(metrics)) {
    metric <- metrics[i]
    title <- metric_titles[i]
    
    # Filter data for this metric
    metric_data <- sensitivity_results[sensitivity_results[[metric]] != 0, ]
    
    if (nrow(metric_data) > 0) {
      p <- ggplot(metric_data, aes_string(x = "value", y = metric, color = "parameter")) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        labs(title = paste(title, "Sensitivity"),
             x = "Parameter Value",
             y = title,
             color = "Parameter") +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      plots[[i]] <- p
    }
  }
  
  # Remove NULL plots
  plots <- plots[!sapply(plots, is.null)]
  
  # Combine plots
  if (length(plots) > 0) {
    combined_plot <- do.call(grid.arrange, c(plots, ncol = 2))
    
    # Save plot if path provided
    if (!is.null(save_path)) {
      ggsave(save_path, combined_plot, width = 15, height = 10, dpi = 300)
    }
    
    return(combined_plot)
  }
  
  return(NULL)
}

# Intervention effectiveness visualization
plot_intervention_effectiveness <- function(intervention_results, save_path = NULL) {
  # Extract data for plotting
  intervention_types <- sapply(intervention_results, function(x) x$intervention_type)
  peak_reductions <- sapply(intervention_results, function(x) x$reduction_in_peak_infection * 100)
  attack_reductions <- sapply(intervention_results, function(x) x$reduction_in_attack_rate * 100)
  
  # Create data frame
  df <- data.frame(
    intervention = intervention_types,
    peak_reduction = peak_reductions,
    attack_reduction = attack_reductions
  )
  
  # Create bar plots
  p1 <- ggplot(df, aes(x = intervention, y = peak_reduction)) +
    geom_bar(stat = "identity", fill = "skyblue", alpha = 0.8) +
    geom_text(aes(label = paste0(round(peak_reduction, 1), "%")), 
              vjust = -0.5, size = 3) +
    labs(title = "Intervention Effectiveness: Peak Infection Reduction",
         x = "Intervention Type",
         y = "Reduction (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(df, aes(x = intervention, y = attack_reduction)) +
    geom_bar(stat = "identity", fill = "lightcoral", alpha = 0.8) +
    geom_text(aes(label = paste0(round(attack_reduction, 1), "%")), 
              vjust = -0.5, size = 3) +
    labs(title = "Intervention Effectiveness: Attack Rate Reduction",
         x = "Intervention Type",
         y = "Reduction (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, ncol = 2)
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = 12, height = 6, dpi = 300)
  }
  
  return(combined_plot)
}

# Age-structured model visualization
plot_age_structured_results <- function(model, solution, save_path = NULL) {
  # Extract compartments
  S <- solution[, 2:(model$n_age_groups + 1)]
  I <- solution[, (model$n_age_groups + 2):(2 * model$n_age_groups + 1)]
  R <- solution[, (2 * model$n_age_groups + 2):(3 * model$n_age_groups + 1)]
  time <- solution[, 1]
  
  # Create data frames for plotting
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
  
  p5 <- ggplot(df_r0, aes(x = age_group, y = R0)) +
    geom_bar(stat = "identity", fill = "red") +
    geom_text(aes(label = sprintf("%.2f", R0)), vjust = -0.5) +
    labs(title = "Age-Specific Basic Reproduction Numbers",
         x = "Age Group",
         y = "Basic Reproduction Number (R₀)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, p4, p5, ncol = 2)
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = 15, height = 12, dpi = 300)
  }
  
  return(combined_plot)
}

# Vector dynamics model visualization
plot_vector_dynamics_results <- function(model, solution, save_path = NULL) {
  # Extract compartments
  S_h <- solution[, "S_h"]
  E_h <- solution[, "E_h"]
  I_h <- solution[, "I_h"]
  R_h <- solution[, "R_h"]
  S_v <- solution[, "S_v"]
  E_v <- solution[, "E_v"]
  I_v <- solution[, "I_v"]
  time <- solution[, "time"]
  
  # Create data frame
  df <- data.frame(
    time = time,
    S_h = S_h, E_h = E_h, I_h = I_h, R_h = R_h,
    S_v = S_v, E_v = E_v, I_v = I_v
  )
  
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
  df$N_v_seasonal <- sapply(df$time, function(t) model$N_v * (1 + model$seasonal_amplitude * sin(2 * pi * t / 365 + model$seasonal_phase)))
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

# Comprehensive analysis function
create_comprehensive_analysis <- function(models, t_span = c(0, 1000), save_dir = "results/plots/") {
  cat("Creating comprehensive analysis...\n")
  
  # Model comparison
  cat("1. Creating model comparison plots...\n")
  plot_model_comparison(models, t_span, save_path = paste0(save_dir, "model_comparison.png"))
  
  # Individual model plots
  for (model_name in names(models)) {
    cat(paste("2. Creating plots for", model_name, "...\n"))
    
    model <- models[[model_name]]
    solution <- simulate(model, t_span)
    
    if (inherits(model, "AgeStructuredMalariaSIR")) {
      plot_age_structured_results(model, solution, 
                                save_path = paste0(save_dir, model_name, "_age_structured.png"))
    } else if (inherits(model, "MalariaVectorModel")) {
      plot_vector_dynamics_results(model, solution, 
                                 save_path = paste0(save_dir, model_name, "_vector_dynamics.png"))
    } else {
      # Basic SIR or SEIR
      plot_results(model, solution, 
                  save_path = paste0(save_dir, model_name, "_results.png"))
    }
  }
  
  cat("Analysis complete! Plots saved to", save_dir, "\n")
}

# Example usage function
demo_visualization <- function() {
  cat("Malaria Model Visualization Demo (R Version)\n")
  cat("===========================================\n\n")
  
  # Load model files
  source("models/basic_sir.R")
  source("models/seir_model.R")
  source("models/age_structured_sir.R")
  source("models/vector_dynamics.R")
  
  # Create sample models
  models <- list(
    "Basic SIR" = MalariaSIRModel(beta = 0.4, gamma = 0.05, N = 10000, I0 = 50),
    "SEIR" = MalariaSEIRModel(beta = 0.4, sigma = 0.1, gamma = 0.05, N = 10000, I0 = 30, E0 = 20),
    "Age-Structured SIR" = AgeStructuredMalariaSIR(N_total = 10000, I0 = c(5, 10, 15, 10, 5, 5)),
    "Vector Dynamics" = MalariaVectorModel(N_h = 10000, N_v = 50000, I_h0 = 30, E_h0 = 20)
  )
  
  # Run comprehensive analysis
  create_comprehensive_analysis(models)
  
  cat("Demo completed! Check the results/plots/ directory for generated visualizations.\n")
}

# Run demo if script is executed directly
if (interactive()) {
  demo_visualization()
}


