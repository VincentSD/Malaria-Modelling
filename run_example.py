#!/usr/bin/env python3
"""
Malaria Model Simulation - Example Runner
========================================

This script demonstrates how to use the malaria transmission models
and provides examples of different analyses that can be performed.

Usage:
    python run_example.py
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add the python directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), 'python'))

from models.basic_sir import MalariaSIRModel
from models.seir_model import MalariaSEIRModel
from models.age_structured_sir import AgeStructuredMalariaSIR
from models.vector_dynamics import MalariaVectorModel
from analysis.visualization import MalariaVisualizer, create_comprehensive_analysis
from analysis.parameter_analysis import comprehensive_parameter_analysis


def run_basic_sir_example():
    """Run basic SIR model example."""
    print("=" * 60)
    print("BASIC SIR MODEL EXAMPLE")
    print("=" * 60)
    
    # Create basic SIR model
    model = MalariaSIRModel(
        beta=0.4,      # Transmission rate
        gamma=0.05,    # Recovery rate (20-day infectious period)
        mu=0.0001,     # Birth/death rate
        N=10000,       # Population of 10,000
        I0=50,         # Start with 50 infected
        R0=0           # No initially recovered
    )
    
    print(f"Model Parameters:")
    print(f"  Transmission rate (β): {model.beta}")
    print(f"  Recovery rate (γ): {model.gamma}")
    print(f"  Birth/death rate (μ): {model.mu}")
    print(f"  Population size (N): {model.N}")
    print(f"  Basic reproduction number (R₀): {model.R0_basic:.2f}")
    print()
    
    # Simulate the model
    t, solution = model.simulate(t_span=(0, 1000), num_points=1000)
    
    # Calculate and display metrics
    metrics = model.get_epidemic_metrics(solution)
    print("Epidemic Metrics:")
    for key, value in metrics.items():
        if 'rate' in key or 'R0' in key:
            print(f"  {key}: {value:.3f}")
        else:
            print(f"  {key}: {value:.0f}")
    print()
    
    # Plot results
    model.plot_results(t, solution, save_path="results/plots/basic_sir_example.png")
    
    return model, t, solution


def run_seir_example():
    """Run SEIR model example."""
    print("=" * 60)
    print("SEIR MODEL EXAMPLE")
    print("=" * 60)
    
    # Create SEIR model
    model = MalariaSEIRModel(
        beta=0.4,                    # Transmission rate
        sigma=0.1,                   # Incubation rate (10-day latent period)
        gamma=0.05,                  # Recovery rate (20-day infectious period)
        mu=0.0001,                   # Birth/death rate
        N=10000,                     # Population of 10,000
        E0=20,                       # Start with 20 exposed
        I0=30,                       # Start with 30 infected
        R0=0,                        # No initially recovered
        seasonal_amplitude=0.3,      # 30% seasonal variation
        seasonal_phase=0.0           # Peak in summer
    )
    
    print(f"Model Parameters:")
    print(f"  Transmission rate (β): {model.beta}")
    print(f"  Incubation rate (σ): {model.sigma}")
    print(f"  Recovery rate (γ): {model.gamma}")
    print(f"  Birth/death rate (μ): {model.mu}")
    print(f"  Population size (N): {model.N}")
    print(f"  Seasonal amplitude: {model.seasonal_amplitude}")
    print(f"  Basic reproduction number (R₀): {model.R0_basic:.2f}")
    print(f"  Effective reproduction number (R₀): {model.R0_effective:.2f}")
    print()
    
    # Simulate the model
    t, solution = model.simulate(t_span=(0, 1000), num_points=1000)
    
    # Calculate and display metrics
    metrics = model.get_epidemic_metrics(solution)
    print("Epidemic Metrics:")
    for key, value in metrics.items():
        if 'rate' in key or 'R0' in key or 'period' in key:
            print(f"  {key}: {value:.3f}")
        else:
            print(f"  {key}: {value:.0f}")
    print()
    
    # Plot results
    model.plot_results(t, solution, save_path="results/plots/seir_example.png")
    
    return model, t, solution


def run_age_structured_example():
    """Run age-structured SIR model example."""
    print("=" * 60)
    print("AGE-STRUCTURED SIR MODEL EXAMPLE")
    print("=" * 60)
    
    # Create age-structured model
    model = AgeStructuredMalariaSIR(
        N_total=100000,
        I0=np.array([10, 20, 30, 20, 10, 10])  # Initial infected by age
    )
    
    print(f"Model Parameters:")
    print(f"  Number of age groups: {model.n_age_groups}")
    print(f"  Age groups: {model.age_groups}")
    print(f"  Total population: {model.N_total:,}")
    print(f"  Age distribution: {model.age_distribution}")
    print()
    
    print("Age-specific R₀ values:")
    for i, age_group in enumerate(model.age_groups):
        print(f"  Age {age_group}: {model.R0_age[i]:.2f}")
    print()
    
    # Simulate the model
    t, solution = model.simulate(t_span=(0, 1000), num_points=1000)
    
    # Calculate and display metrics
    metrics = model.get_age_specific_metrics(solution)
    
    print("Age-specific epidemic metrics:")
    print("-" * 40)
    for i, age_group in enumerate(model.age_groups):
        print(f"Age {age_group}:")
        print(f"  Peak infection rate: {metrics['peak_infection_rate'][i]:.3f}")
        print(f"  Final attack rate: {metrics['final_attack_rate'][i]:.3f}")
        print(f"  Time to peak: {metrics['time_to_peak'][i]} days")
        print()
    
    # Plot results
    model.plot_age_structured_results(t, solution, 
                                    save_path="results/plots/age_structured_example.png")
    
    return model, t, solution


def run_vector_dynamics_example():
    """Run vector dynamics model example."""
    print("=" * 60)
    print("VECTOR DYNAMICS MODEL EXAMPLE")
    print("=" * 60)
    
    # Create vector dynamics model
    model = MalariaVectorModel(
        # Human parameters
        beta_h=0.3, sigma_h=0.1, gamma_h=0.05, mu_h=0.0001, N_h=10000,
        # Vector parameters
        beta_v=0.2, sigma_v=0.2, mu_v=0.1, N_v=50000,
        # Biting parameters
        a=0.3, b=0.5, c=0.5,
        # Seasonal parameters
        seasonal_amplitude=0.5, seasonal_phase=0.0,
        # Initial conditions
        E_h0=20, I_h0=30, R_h0=0, E_v0=500, I_v0=1000
    )
    
    print(f"Model Parameters:")
    print(f"  Human population: {model.N_h:,}")
    print(f"  Vector population: {model.N_v:,}")
    print(f"  Vector-human ratio: {model.N_v/model.N_h:.1f}")
    print(f"  Biting rate (a): {model.a}")
    print(f"  Human to vector transmission (b): {model.b}")
    print(f"  Vector to human transmission (c): {model.c}")
    print(f"  Basic reproduction number (R₀): {model.R0:.2f}")
    print()
    
    # Simulate the model
    t, solution = model.simulate(t_span=(0, 1000), num_points=1000)
    
    # Calculate and display metrics
    metrics = model.get_epidemic_metrics(solution)
    print("Epidemic Metrics:")
    for key, value in metrics.items():
        if 'rate' in key or 'R0' in key or 'ratio' in key:
            print(f"  {key}: {value:.3f}")
        else:
            print(f"  {key}: {value:.0f}")
    print()
    
    # Plot results
    model.plot_results(t, solution, save_path="results/plots/vector_dynamics_example.png")
    
    # Intervention analysis
    print("Intervention Analysis:")
    print("-" * 20)
    
    interventions = [
        ('insecticide', 0.5, 200, 100),
        ('bed_nets', 0.7, 200, 100),
        ('larvicide', 0.6, 200, 100)
    ]
    
    for intervention_type, strength, start, duration in interventions:
        effectiveness = model.intervention_analysis(
            intervention_type, strength, start, duration
        )
        
        print(f"{intervention_type.title()} (strength: {strength}):")
        print(f"  Reduction in peak infection: {effectiveness['reduction_in_peak_infection']:.1%}")
        print(f"  Reduction in attack rate: {effectiveness['reduction_in_attack_rate']:.1%}")
        print()
    
    return model, t, solution


def run_comprehensive_analysis():
    """Run comprehensive analysis comparing all models."""
    print("=" * 60)
    print("COMPREHENSIVE MODEL COMPARISON")
    print("=" * 60)
    
    # Create all models
    models = {
        'Basic SIR': MalariaSIRModel(beta=0.4, gamma=0.05, N=10000, I0=50),
        'SEIR': MalariaSEIRModel(beta=0.4, sigma=0.1, gamma=0.05, N=10000, I0=30, E0=20),
        'Age-Structured SIR': AgeStructuredMalariaSIR(N_total=10000, I0=np.array([5, 10, 15, 10, 5, 5])),
        'Vector Dynamics': MalariaVectorModel(N_h=10000, N_v=50000, I_h0=30, E_h0=20)
    }
    
    # Run comprehensive analysis
    create_comprehensive_analysis(models, t_span=(0, 1000), save_dir="results/plots/")
    
    return models


def run_parameter_analysis():
    """Run parameter sensitivity analysis."""
    print("=" * 60)
    print("PARAMETER SENSITIVITY ANALYSIS")
    print("=" * 60)
    
    # Define model parameters
    model_params = {
        'beta': 0.4,
        'gamma': 0.05,
        'mu': 0.0001,
        'N': 10000,
        'I0': 50,
        'R0': 0
    }
    
    # Define parameter ranges for analysis
    param_ranges = {
        'beta': (0.1, 0.8),
        'gamma': (0.01, 0.2),
        'mu': (0.00005, 0.0005)
    }
    
    # Run comprehensive analysis
    comprehensive_parameter_analysis(
        MalariaSIRModel, 
        model_params, 
        param_ranges,
        save_dir="results/analysis/"
    )


def main():
    """Main function to run all examples."""
    print("MALARIA MODEL SIMULATION PROJECT")
    print("================================")
    print("This script demonstrates various malaria transmission models")
    print("and analysis techniques.")
    print()
    
    # Create results directories
    os.makedirs("results/plots", exist_ok=True)
    os.makedirs("results/analysis", exist_ok=True)
    
    try:
        # Run individual model examples
        print("Running individual model examples...")
        print()
        
        basic_model, basic_t, basic_solution = run_basic_sir_example()
        print()
        
        seir_model, seir_t, seir_solution = run_seir_example()
        print()
        
        age_model, age_t, age_solution = run_age_structured_example()
        print()
        
        vector_model, vector_t, vector_solution = run_vector_dynamics_example()
        print()
        
        # Run comprehensive analysis
        print("Running comprehensive model comparison...")
        print()
        models = run_comprehensive_analysis()
        print()
        
        # Run parameter analysis
        print("Running parameter sensitivity analysis...")
        print()
        run_parameter_analysis()
        print()
        
        print("=" * 60)
        print("ALL EXAMPLES COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        print("Check the 'results/' directory for generated plots and analysis files.")
        print()
        print("Generated files:")
        print("- results/plots/ - All visualization plots")
        print("- results/analysis/ - Parameter analysis results")
        print()
        
    except Exception as e:
        print(f"Error running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()


