# Parameter Analysis Tools for Malaria Models (R Version)
# ======================================================
#
# This script provides comprehensive parameter analysis tools for malaria transmission models,
# including sensitivity analysis, parameter estimation, and uncertainty quantification.

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

# Parameter sensitivity analysis
parameter_sensitivity_analysis <- function(model, param_ranges, num_points = 10, 
                                         t_span = c(0, 365), metrics = NULL) {
  if (is.null(metrics)) {
    metrics <- c("peak_infection_rate", "final_attack_rate", "time_to_peak")
  }
  
  results <- data.frame()
  
  for (param_name in names(param_ranges)) {
    param_range <- param_ranges[[param_name]]
    param_values <- seq(param_range[1], param_range[2], length.out = num_points)
    
    cat(paste("Analyzing parameter:", param_name, "\n"))
    
    for (value in param_values) {
      # Create model with modified parameter
      modified_params <- list()
      modified_params[[param_name]] <- value
      
      # Create temporary model
      temp_model <- do.call(class(model)[1], c(as.list(model), modified_params))
      
      tryCatch({
        # Simulate model
        solution <- simulate(temp_model, t_span)
        
        # Calculate metrics
        model_metrics <- get_epidemic_metrics(temp_model, solution)
        
        # Store results
        result_row <- data.frame(
          parameter = param_name,
          value = value,
          normalized_value = (value - param_range[1]) / (param_range[2] - param_range[1])
        )
        
        for (metric in metrics) {
          if (metric %in% names(model_metrics)) {
            result_row[[metric]] <- model_metrics[[metric]]
          } else {
            result_row[[metric]] <- NA
          }
        }
        
        results <- rbind(results, result_row)
        
      }, error = function(e) {
        cat(paste("Error with", param_name, "=", value, ":", e$message, "\n"))
      })
    }
  }
  
  return(results)
}

# Plot sensitivity results
plot_sensitivity_results <- function(sensitivity_df, save_path = NULL) {
  # Get unique parameters and metrics
  parameters <- unique(sensitivity_df$parameter)
  metrics <- c("peak_infection_rate", "final_attack_rate", "time_to_peak")
  metric_titles <- c("Peak Infection Rate", "Final Attack Rate", "Time to Peak")
  
  # Create subplots
  plots <- list()
  
  for (i in seq_along(metrics)) {
    metric <- metrics[i]
    title <- metric_titles[i]
    
    # Filter data for this metric
    metric_data <- sensitivity_df[!is.na(sensitivity_df[[metric]]), ]
    
    if (nrow(metric_data) > 0) {
      p <- ggplot(metric_data, aes_string(x = "normalized_value", y = metric, color = "parameter")) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        labs(title = paste(title, "Sensitivity"),
             x = "Normalized Parameter Value",
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

# Parameter estimation using optimization
parameter_estimation <- function(observed_data, param_bounds, model_class, 
                               t_span = c(0, 365), method = "L-BFGS-B") {
  # Objective function for parameter estimation
  objective_function <- function(params) {
    # Create model with current parameters
    model <- do.call(model_class, params)
    
    tryCatch({
      # Simulate model
      solution <- simulate(model, t_span)
      
      # Extract infected data
      if (inherits(model, "MalariaSIRModel")) {
        I_pred <- solution[, "I"]
      } else if (inherits(model, "MalariaSEIRModel")) {
        I_pred <- solution[, "I"]
      } else if (inherits(model, "AgeStructuredMalariaSIR")) {
        I_pred <- rowSums(solution[, (model$n_age_groups + 2):(2 * model$n_age_groups + 1)])
      } else if (inherits(model, "MalariaVectorModel")) {
        I_pred <- solution[, "I_h"]
      }
      
      # Interpolate to match observed data time points
      I_pred_interp <- approx(solution[, "time"], I_pred, observed_data$time)$y
      
      # Calculate sum of squared errors
      sse <- sum((observed_data$infected - I_pred_interp)^2, na.rm = TRUE)
      
      return(sse)
      
    }, error = function(e) {
      return(1e10)  # Large penalty for invalid parameters
    })
  }
  
  # Prepare bounds for optimization
  param_names <- names(param_bounds)
  lower_bounds <- sapply(param_bounds, function(x) x[1])
  upper_bounds <- sapply(param_bounds, function(x) x[2])
  
  # Initial guess (middle of bounds)
  x0 <- sapply(param_bounds, function(x) (x[1] + x[2]) / 2)
  names(x0) <- param_names
  
  # Run optimization
  result <- optim(x0, objective_function, method = method,
                  lower = lower_bounds, upper = upper_bounds,
                  control = list(maxit = 1000))
  
  # Create results
  estimated_params <- result$par
  names(estimated_params) <- param_names
  
  return(list(
    estimated_parameters = estimated_params,
    optimization_success = result$convergence == 0,
    objective_value = result$value,
    optimization_result = result
  ))
}

# Uncertainty quantification using Monte Carlo sampling
uncertainty_quantification <- function(model, param_distributions, num_samples = 1000, 
                                     t_span = c(0, 365)) {
  results <- data.frame()
  
  for (i in 1:num_samples) {
    if (i %% 100 == 0) {
      cat(paste("Monte Carlo sample", i, "/", num_samples, "\n"))
    }
    
    # Sample parameters
    sampled_params <- list()
    
    for (param_name in names(param_distributions)) {
      dist_info <- param_distributions[[param_name]]
      dist_type <- dist_info$type
      param1 <- dist_info$param1
      param2 <- dist_info$param2
      
      if (dist_type == "uniform") {
        value <- runif(1, param1, param2)
      } else if (dist_type == "normal") {
        value <- rnorm(1, param1, param2)
      } else {
        stop(paste("Unknown distribution type:", dist_type))
      }
      
      sampled_params[[param_name]] <- value
    }
    
    tryCatch({
      # Create and simulate model
      temp_model <- do.call(class(model)[1], c(as.list(model), sampled_params))
      solution <- simulate(temp_model, t_span)
      
      # Calculate metrics
      metrics <- get_epidemic_metrics(temp_model, solution)
      
      # Store results
      result_row <- data.frame(sample = i)
      result_row <- cbind(result_row, as.data.frame(sampled_params))
      result_row <- cbind(result_row, as.data.frame(metrics))
      
      results <- rbind(results, result_row)
      
    }, error = function(e) {
      cat(paste("Error in sample", i, ":", e$message, "\n"))
    })
  }
  
  return(results)
}

# Plot uncertainty results
plot_uncertainty_results <- function(uncertainty_df, model_params, save_path = NULL) {
  # Get parameter and metric columns
  param_cols <- names(model_params)
  metric_cols <- setdiff(names(uncertainty_df), c("sample", param_cols))
  
  # Create subplots
  n_metrics <- length(metric_cols)
  n_cols <- min(3, n_metrics)
  n_rows <- ceiling(n_metrics / n_cols)
  
  plots <- list()
  
  for (i in seq_along(metric_cols)) {
    metric <- metric_cols[i]
    
    p <- ggplot(uncertainty_df, aes_string(x = metric)) +
      geom_histogram(bins = 50, alpha = 0.7, fill = "skyblue") +
      geom_vline(aes_string(xintercept = paste0("mean(", metric, ")")), 
                color = "red", linestyle = "dashed", size = 1) +
      geom_vline(aes_string(xintercept = paste0("median(", metric, ")")), 
                color = "green", linestyle = "dashed", size = 1) +
      labs(title = paste("Distribution of", metric),
           x = metric,
           y = "Frequency") +
      theme_minimal()
    
    plots[[i]] <- p
  }
  
  # Combine plots
  if (length(plots) > 0) {
    combined_plot <- do.call(grid.arrange, c(plots, ncol = n_cols))
    
    # Save plot if path provided
    if (!is.null(save_path)) {
      ggsave(save_path, combined_plot, width = 5 * n_cols, height = 4 * n_rows, dpi = 300)
    }
    
    return(combined_plot)
  }
  
  return(NULL)
}

# Correlation analysis
correlation_analysis <- function(uncertainty_df, model_params) {
  # Get parameter and metric columns
  param_cols <- names(model_params)
  metric_cols <- setdiff(names(uncertainty_df), c("sample", param_cols))
  
  # Calculate correlations
  correlation_data <- data.frame()
  
  for (param in param_cols) {
    for (metric in metric_cols) {
      corr <- cor(uncertainty_df[[param]], uncertainty_df[[metric]], use = "complete.obs")
      correlation_data <- rbind(correlation_data, data.frame(
        parameter = param,
        metric = metric,
        correlation = corr,
        abs_correlation = abs(corr)
      ))
    }
  }
  
  # Create correlation heatmap
  pivot_table <- reshape2::dcast(correlation_data, parameter ~ metric, value.var = "correlation")
  rownames(pivot_table) <- pivot_table$parameter
  pivot_table$parameter <- NULL
  
  # Convert to matrix for heatmap
  corr_matrix <- as.matrix(pivot_table)
  
  # Create heatmap
  p <- ggplot(reshape2::melt(corr_matrix), aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0, limit = c(-1, 1), space = "Lab", 
                        name = "Correlation") +
    labs(title = "Parameter-Metric Correlation Matrix",
         x = "Parameter",
         y = "Metric") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(correlation_data = correlation_data, heatmap = p))
}

# Comprehensive parameter analysis
comprehensive_parameter_analysis <- function(model, param_ranges, save_dir = "results/analysis/") {
  cat("Starting comprehensive parameter analysis...\n")
  
  # Create save directory if it doesn't exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # 1. Sensitivity analysis
  cat("1. Performing sensitivity analysis...\n")
  sensitivity_results <- parameter_sensitivity_analysis(model, param_ranges)
  
  # Plot sensitivity results
  plot_sensitivity_results(sensitivity_results, 
                          save_path = paste0(save_dir, "sensitivity_analysis.png"))
  
  # Save sensitivity results
  write.csv(sensitivity_results, paste0(save_dir, "sensitivity_results.csv"), row.names = FALSE)
  
  # 2. Uncertainty quantification
  cat("2. Performing uncertainty quantification...\n")
  
  # Create parameter distributions
  param_distributions <- list()
  for (param_name in names(param_ranges)) {
    param_range <- param_ranges[[param_name]]
    param_distributions[[param_name]] <- list(
      type = "uniform",
      param1 = param_range[1],
      param2 = param_range[2]
    )
  }
  
  uncertainty_results <- uncertainty_quantification(model, param_distributions, num_samples = 500)
  
  # Plot uncertainty results
  plot_uncertainty_results(uncertainty_results, as.list(model), 
                          save_path = paste0(save_dir, "uncertainty_analysis.png"))
  
  # Save uncertainty results
  write.csv(uncertainty_results, paste0(save_dir, "uncertainty_results.csv"), row.names = FALSE)
  
  # 3. Correlation analysis
  cat("3. Performing correlation analysis...\n")
  correlation_results <- correlation_analysis(uncertainty_results, as.list(model))
  
  # Save correlation results
  write.csv(correlation_results$correlation_data, 
            paste0(save_dir, "correlation_analysis.csv"), row.names = FALSE)
  
  # Save correlation heatmap
  ggsave(paste0(save_dir, "correlation_heatmap.png"), 
         correlation_results$heatmap, width = 10, height = 8, dpi = 300)
  
  cat("Analysis complete! Results saved to", save_dir, "\n")
  
  return(list(
    sensitivity_results = sensitivity_results,
    uncertainty_results = uncertainty_results,
    correlation_results = correlation_results
  ))
}

# Example usage function
demo_parameter_analysis <- function() {
  cat("Malaria Model Parameter Analysis Demo (R Version)\n")
  cat("================================================\n\n")
  
  # Load model files
  source("models/basic_sir.R")
  
  # Create sample model
  model <- MalariaSIRModel(beta = 0.4, gamma = 0.05, mu = 0.0001, N = 10000, I0 = 50, R0 = 0)
  
  # Define parameter ranges for analysis
  param_ranges <- list(
    beta = c(0.1, 0.8),
    gamma = c(0.01, 0.2),
    mu = c(0.00005, 0.0005)
  )
  
  # Run comprehensive analysis
  results <- comprehensive_parameter_analysis(model, param_ranges)
  
  cat("Demo completed! Check the results/analysis/ directory for generated files.\n")
}

# Run demo if script is executed directly
if (interactive()) {
  demo_parameter_analysis()
}


