# R Package Dependencies for Malaria Model Simulation Project
# ==========================================================

# Core packages for differential equations and simulation
if (!require(deSolve)) {
  install.packages("deSolve")
  library(deSolve)
}

# Data manipulation and analysis
if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

if (!require(tidyr)) {
  install.packages("tidyr")
  library(tidyr)
}

# Visualization
if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(gridExtra)) {
  install.packages("gridExtra")
  library(gridExtra)
}

if (!require(plotly)) {
  install.packages("plotly")
  library(plotly)
}

# Statistical analysis
if (!require(stats)) {
  install.packages("stats")
  library(stats)
}

# Matrix operations and linear algebra
if (!require(Matrix)) {
  install.packages("Matrix")
  library(Matrix)
}

# Optional: Advanced visualization
if (!require(ggthemes)) {
  install.packages("ggthemes")
  library(ggthemes)
}

if (!require(viridis)) {
  install.packages("viridis")
  library(viridis)
}

# Optional: Interactive plots
if (!require(plotly)) {
  install.packages("plotly")
  library(plotly)
}

# Optional: Parallel processing
if (!require(parallel)) {
  install.packages("parallel")
  library(parallel)
}

# Optional: Progress bars
if (!require(pbapply)) {
  install.packages("pbapply")
  library(pbapply)
}

# Optional: Data export
if (!require(openxlsx)) {
  install.packages("openxlsx")
  library(openxlsx)
}

# Print package versions for reproducibility
cat("R Package Versions:\n")
cat("==================\n")
cat("R version:", R.version.string, "\n")
cat("deSolve:", packageVersion("deSolve"), "\n")
cat("dplyr:", packageVersion("dplyr"), "\n")
cat("ggplot2:", packageVersion("ggplot2"), "\n")
cat("gridExtra:", packageVersion("gridExtra"), "\n")


