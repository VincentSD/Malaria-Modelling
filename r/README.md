# Malaria Model Simulation Project - R Implementation
# ==================================================

This directory contains the R implementation of the malaria model simulation project, providing comprehensive epidemiological modeling tools for malaria transmission dynamics.

## 📁 Project Structure

```
r/
├── models/                    # Core malaria transmission models
│   ├── basic_sir.R           # Basic SIR model
│   ├── seir_model.R          # SEIR model with latent period
│   ├── age_structured_sir.R  # Age-structured SIR model
│   └── vector_dynamics.R     # Vector dynamics model
├── analysis/                 # Analysis and visualization tools
│   ├── visualization.R       # Advanced visualization functions
│   └── parameter_analysis.R  # Parameter sensitivity analysis
├── run_example.R            # Complete example script
├── packages.R               # Required R packages
└── README.md               # This file
```

## 🚀 Quick Start

### Prerequisites

- R (version 3.6 or higher)
- Required R packages (see `packages.R`)

### Installation

1. **Install required packages:**
   ```r
   source("packages.R")
   ```

2. **Run the complete example:**
   ```r
   source("run_example.R")
   ```

3. **Or run individual models:**
   ```r
   # Basic SIR model
   source("models/basic_sir.R")
   model <- MalariaSIRModel(beta = 0.4, gamma = 0.05, N = 10000, I0 = 50)
   solution <- simulate(model, t_span = c(0, 365))
   plot_results(model, solution)
   ```

## 🧬 Implemented Models

### 1. Basic SIR Model (`basic_sir.R`)
- **Purpose**: Foundational malaria transmission model
- **Compartments**: Susceptible (S), Infected (I), Recovered (R)
- **Key Features**:
  - Malaria-specific parameters
  - Epidemic metrics calculation
  - Visualization tools
  - Parameter sensitivity analysis

### 2. SEIR Model (`seir_model.R`)
- **Purpose**: Enhanced model with latent period
- **Compartments**: Susceptible (S), Exposed (E), Infected (I), Recovered (R)
- **Key Features**:
  - Incubation period modeling
  - Seasonal transmission variations
  - Comparison with SIR model
  - Latent period analysis

### 3. Age-Structured SIR Model (`age_structured_sir.R`)
- **Purpose**: Population heterogeneity modeling
- **Compartments**: Age-specific S, I, R compartments
- **Key Features**:
  - Multiple age groups
  - Age-specific contact matrices
  - Age-specific transmission rates
  - Demographic considerations

### 4. Vector Dynamics Model (`vector_dynamics.R`)
- **Purpose**: Complete human-vector transmission cycle
- **Compartments**: Human (S, E, I, R) + Vector (S, E, I)
- **Key Features**:
  - Mosquito population dynamics
  - Human-vector transmission
  - Seasonal vector populations
  - Intervention modeling

## 📊 Analysis Tools

### Visualization (`analysis/visualization.R`)
- Model comparison plots
- Parameter sensitivity visualization
- Intervention effectiveness plots
- Age-structured model visualization
- Vector dynamics visualization

### Parameter Analysis (`analysis/parameter_analysis.R`)
- Parameter sensitivity analysis
- Parameter estimation
- Uncertainty quantification
- Correlation analysis
- Monte Carlo sampling

## 🔬 Key Features

### Model Capabilities
- **Multiple Model Types**: SIR, SEIR, Age-Structured, Vector Dynamics
- **Realistic Parameters**: Malaria-specific transmission rates, recovery rates, etc.
- **Seasonal Variations**: Annual transmission patterns
- **Age Structure**: Population heterogeneity modeling
- **Vector Dynamics**: Complete human-mosquito transmission cycle

### Analysis Features
- **Epidemic Metrics**: R₀, peak infection, attack rate, time to peak
- **Sensitivity Analysis**: Parameter impact assessment
- **Intervention Modeling**: Insecticide, bed nets, larvicide effects
- **Uncertainty Quantification**: Monte Carlo analysis
- **Visualization**: Comprehensive plotting tools

### Scientific Accuracy
- **Malaria-Specific**: Parameters based on malaria epidemiology
- **Realistic Dynamics**: Incubation periods, recovery rates, transmission cycles
- **Population Structure**: Age-specific transmission patterns
- **Vector Biology**: Mosquito life cycle and transmission dynamics

## 📈 Example Usage

### Basic Model Simulation
```r
# Create and simulate a basic SIR model
model <- MalariaSIRModel(beta = 0.4, gamma = 0.05, N = 10000, I0 = 50)
solution <- simulate(model, t_span = c(0, 365))
metrics <- get_epidemic_metrics(model, solution)
plot_results(model, solution)
```

### Model Comparison
```r
# Compare different models
models <- list(
  "SIR" = MalariaSIRModel(beta = 0.4, gamma = 0.05, N = 10000, I0 = 50),
  "SEIR" = MalariaSEIRModel(beta = 0.4, sigma = 0.1, gamma = 0.05, N = 10000, I0 = 30, E0 = 20)
)
plot_model_comparison(models)
```

### Parameter Sensitivity Analysis
```r
# Analyze parameter sensitivity
param_ranges <- list(beta = c(0.1, 0.8), gamma = c(0.01, 0.2))
sensitivity_results <- parameter_sensitivity_analysis(model, param_ranges)
plot_sensitivity_results(sensitivity_results)
```

### Intervention Analysis
```r
# Test intervention effectiveness
effectiveness <- intervention_analysis(
  vector_model, "bed_nets", 0.7, 200, 100
)
```

## 🧪 Experimental Setups

The project includes comprehensive experimental scenarios (see `../EXPERIMENTAL_SETUPS.md`):

1. **R₀ Analysis**: Understanding basic reproduction number
2. **Recovery Rate Impact**: Effect of treatment on epidemic dynamics
3. **Latent Period Effects**: Incubation period influence
4. **Seasonal Transmission**: Annual variation patterns
5. **Age-Specific Transmission**: Population heterogeneity effects
6. **Vector-Human Ratio**: Mosquito population impact
7. **Intervention Effectiveness**: Control strategy comparison
8. **Parameter Uncertainty**: Sensitivity and robustness analysis
9. **Model Comparison**: SIR vs SEIR vs Age-Structured vs Vector
10. **Real-World Scenarios**: Country-specific parameter sets

## 📚 Scientific Background

### Malaria Transmission
- **Vector-Borne Disease**: Transmitted by Anopheles mosquitoes
- **Plasmodium Parasites**: Multiple species (P. falciparum, P. vivax, etc.)
- **Transmission Cycle**: Human-mosquito-human cycle
- **Incubation Period**: 7-30 days depending on species
- **Recovery**: Partial immunity, reinfection possible

### Model Assumptions
- **Homogeneous Mixing**: Within age groups and compartments
- **Constant Parameters**: Over simulation time (except seasonal)
- **No Migration**: Closed population systems
- **Perfect Immunity**: Recovered individuals immune to reinfection
- **Vector Dynamics**: Simplified mosquito life cycle

## 🔧 Dependencies

### Required R Packages
- `deSolve`: Differential equation solving
- `ggplot2`: Data visualization
- `dplyr`: Data manipulation
- `gridExtra`: Plot arrangement
- `reshape2`: Data reshaping
- `viridis`: Color scales

### Installation
```r
# Install packages
install.packages(c("deSolve", "ggplot2", "dplyr", "gridExtra", "reshape2", "viridis"))

# Or use the provided script
source("packages.R")
```

## 📖 Documentation

- **Model Documentation**: Each model file contains detailed documentation
- **Function Help**: Use `?function_name` for help on specific functions
- **Examples**: See `run_example.R` for comprehensive examples
- **Parameter Guide**: Check `../data/malaria_parameters.csv` for parameter ranges

## 🤝 Contributing

1. **Code Style**: Follow R coding conventions
2. **Documentation**: Document all functions and parameters
3. **Testing**: Test new features with example data
4. **Validation**: Compare results with Python implementation

## 📄 License

This project is licensed under the MIT License - see the main project README for details.

## ⚠️ Disclaimer

This software is for educational and research purposes only. It should not be used for clinical decision-making or public health policy without proper validation and expert review.

## 🆚 R vs Python Implementation

### R Advantages
- **Statistical Analysis**: Built-in statistical functions
- **Data Visualization**: ggplot2 for publication-quality plots
- **Epidemiological Focus**: Specialized packages for disease modeling
- **Academic Use**: Widely used in epidemiology research

### Python Advantages
- **Web Interface**: Streamlit dashboard
- **Machine Learning**: Advanced analysis capabilities
- **Integration**: Better integration with other tools
- **Performance**: Faster computation for large simulations

### Recommendation
- **Use R** for: Statistical analysis, research, academic work
- **Use Python** for: Web interfaces, machine learning, production systems
- **Use Both** for: Comprehensive analysis and validation

## 🚀 Getting Started

1. **Clone the repository**
2. **Install R and required packages**
3. **Run the example**: `source("run_example.R")`
4. **Explore individual models**
5. **Try the experimental setups**
6. **Compare with Python implementation**

For more information, see the main project README and experimental setups guide.


