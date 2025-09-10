# Malaria Model Simulation Project

This project implements comprehensive malaria transmission models using both R and Python, including basic SIR models and more complex epidemiological models with age structure and vector dynamics. The project provides a complete framework for understanding malaria transmission dynamics, parameter sensitivity analysis, and intervention effectiveness evaluation.

## 🦟 Project Overview

Malaria is a vector-borne disease caused by *Plasmodium* parasites transmitted through *Anopheles* mosquitoes. This project implements multiple mathematical models to simulate malaria transmission dynamics, from simple SIR models to complex age-structured models with vector dynamics.

## 📁 Project Structure

```
MalProj/
├── python/
│   ├── models/
│   │   ├── basic_sir.py              # Basic SIR model
│   │   ├── age_structured_sir.py     # Age-structured SIR model
│   │   ├── seir_model.py             # SEIR model with latent period
│   │   └── vector_dynamics.py        # Vector dynamics model
│   ├── analysis/
│   │   ├── visualization.py          # Advanced visualization tools
│   │   └── parameter_analysis.py     # Parameter sensitivity analysis
│   └── requirements.txt              # Python dependencies
├── r/
│   ├── models/
│   │   ├── basic_sir.R               # Basic SIR model (R version)
│   │   ├── age_structured_sir.R      # Age-structured SIR model (R version)
│   │   ├── seir_model.R              # SEIR model (R version)
│   │   └── vector_dynamics.R         # Vector dynamics model (R version)
│   ├── analysis/
│   │   ├── visualization.R           # Visualization tools (R version)
│   │   └── parameter_analysis.R      # Parameter analysis (R version)
│   └── packages.R                    # R package dependencies
├── data/
│   └── malaria_parameters.csv        # Model parameters database
├── results/
│   ├── plots/                        # Generated visualizations
│   └── analysis/                     # Analysis results
├── run_example.py                    # Complete example runner
└── README.md                         # This file
```

## 🧬 Models Implemented

### 1. Basic SIR Model
- **Purpose**: Foundation model for malaria transmission
- **Compartments**: Susceptible (S) → Infected (I) → Recovered (R)
- **Features**: 
  - Malaria-specific parameters
  - Basic reproduction number (R₀) calculation
  - Population dynamics with birth/death rates
- **Use Case**: Initial analysis and parameter estimation

### 2. Age-Structured SIR Model
- **Purpose**: Realistic population modeling with age-specific dynamics
- **Compartments**: S, I, R for each age group
- **Features**:
  - Multiple age groups (0-4, 5-14, 15-29, 30-44, 45-59, 60+)
  - Age-specific transmission and recovery rates
  - Contact matrices between age groups
  - Age-dependent susceptibility and infectiousness
- **Use Case**: Policy planning, age-targeted interventions

### 3. SEIR Model
- **Purpose**: Includes latent period for more realistic malaria modeling
- **Compartments**: Susceptible (S) → Exposed (E) → Infected (I) → Recovered (R)
- **Features**:
  - Latent period (exposed state)
  - Seasonal variations in transmission
  - Age-specific parameters for different malaria strains
  - Treatment and immunity considerations
- **Use Case**: Detailed transmission dynamics, seasonal analysis

### 4. Vector Dynamics Model
- **Purpose**: Complete human-vector transmission cycle
- **Compartments**: 
  - Human: S, E, I, R
  - Vector: S, E, I
- **Features**:
  - Mosquito population dynamics
  - Human-mosquito transmission cycles
  - Seasonal variations in vector populations
  - Vector control intervention modeling
- **Use Case**: Vector control strategies, comprehensive transmission modeling

## 🔬 Key Malaria Parameters

| Parameter | Symbol | Description | Typical Range |
|-----------|--------|-------------|---------------|
| **Transmission Rate** | β | Rate of infection from infected to susceptible | 0.1 - 0.8 |
| **Recovery Rate** | γ | Rate of recovery from infection | 0.01 - 0.2 |
| **Incubation Rate** | σ | Rate of progression from exposed to infected | 0.05 - 0.3 |
| **Birth/Death Rate** | μ | Population turnover rate | 0.00005 - 0.0005 |
| **Biting Rate** | a | Bites per mosquito per day | 0.1 - 0.5 |
| **Human→Vector Transmission** | b | Probability of transmission from human to mosquito | 0.3 - 0.8 |
| **Vector→Human Transmission** | c | Probability of transmission from mosquito to human | 0.3 - 0.8 |
| **Vector Death Rate** | μᵥ | Mosquito mortality rate | 0.05 - 0.2 |

## 🚀 Quick Start

### Python Installation and Usage

```bash
# Clone or download the project
cd MalProj

# Install Python dependencies
pip install -r python/requirements.txt

# Run the complete example
python run_example.py

# Or run individual models
python python/models/basic_sir.py
python python/models/seir_model.py
python python/models/age_structured_sir.py
python python/models/vector_dynamics.py
```

### R Installation and Usage

```bash
# Install R packages
Rscript r/packages.R

# Run individual models
Rscript r/models/basic_sir.R
Rscript r/models/seir_model.R
Rscript r/models/age_structured_sir.R
Rscript r/models/vector_dynamics.R
```

## 📊 Analysis Features

### Visualization Tools
- **Model Comparison**: Side-by-side comparison of different models
- **Parameter Sensitivity**: Heatmaps and plots showing parameter effects
- **Intervention Analysis**: Effectiveness of different control strategies
- **Interactive Plots**: Plotly-based interactive visualizations
- **Epidemic Curves**: Comprehensive epidemic progression analysis

### Parameter Analysis
- **Sensitivity Analysis**: Systematic parameter variation studies
- **Uncertainty Quantification**: Monte Carlo uncertainty analysis
- **Parameter Estimation**: Fitting models to observed data
- **Correlation Analysis**: Parameter-outcome relationships

### Intervention Modeling
- **Insecticide Spraying**: Vector population reduction
- **Bed Nets**: Biting rate reduction
- **Larvicide**: Vector birth rate reduction
- **Treatment**: Recovery rate enhancement

## 📈 Example Outputs

The project generates comprehensive visualizations including:

1. **Population Dynamics**: Time series of S, I, R compartments
2. **Infection Rates**: Percentage of population infected over time
3. **Phase Planes**: Trajectory plots in state space
4. **Age-Specific Analysis**: Age group comparisons
5. **Vector Dynamics**: Human-vector transmission cycles
6. **Intervention Effectiveness**: Before/after intervention comparisons

## 🔧 Customization

### Adding New Models
1. Create a new model class inheriting from base patterns
2. Implement `simulate()` and `get_epidemic_metrics()` methods
3. Add visualization methods
4. Update the example runner

### Parameter Customization
- Modify `data/malaria_parameters.csv` for parameter ranges
- Adjust model initialization parameters
- Use parameter analysis tools for optimization

### Adding New Interventions
1. Extend the intervention analysis framework
2. Add new intervention types to vector dynamics model
3. Implement effectiveness calculation methods

## 📚 Scientific Background

### Mathematical Framework
The models are based on systems of ordinary differential equations (ODEs) representing the flow between epidemiological compartments. The basic SIR model follows:

```
dS/dt = μN - βSI/N - μS
dI/dt = βSI/N - (γ + μ)I
dR/dt = γI - μR
```

### Malaria-Specific Considerations
- **Vector-borne transmission**: Requires mosquito intermediate host
- **Seasonal patterns**: Vector populations vary with climate
- **Age-dependent susceptibility**: Children more susceptible than adults
- **Immunity**: Partial immunity after recovery
- **Multiple strains**: Different *Plasmodium* species

## 📖 References

1. **Anderson, R.M. & May, R.M.** (1991). *Infectious Diseases of Humans: Dynamics and Control*. Oxford University Press.

2. **Smith, D.L. et al.** (2012). Ross, Macdonald, and a Theory for the Dynamics and Control of Mosquito-Transmitted Pathogens. *PLoS Pathogens*, 8(4), e1002588.

3. **World Health Organization** (2023). *World Malaria Report 2023*. WHO Press.

4. **Griffin, J.T. et al.** (2010). Reducing *Plasmodium falciparum* malaria transmission in Africa: a model-based evaluation of intervention strategies. *PLoS Medicine*, 7(8), e1000324.

5. **Eckhoff, P.A.** (2011). A malaria transmission-directed model of mosquito life cycle and ecology. *Malaria Journal*, 10, 303.

## 🤝 Contributing

Contributions are welcome! Please consider:
- Adding new model variants
- Improving visualization tools
- Adding real-world data integration
- Enhancing parameter estimation methods
- Adding new intervention types

## 📄 License

This project is open source and available under the MIT License.

## ⚠️ Disclaimer

This software is for educational and research purposes. It should not be used for clinical decision-making or public health policy without proper validation and expert consultation.
