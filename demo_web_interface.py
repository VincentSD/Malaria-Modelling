#!/usr/bin/env python3
"""
Malaria Model Simulation - Web Interface Demo
===========================================

This script demonstrates the key features of the web interface
by running a quick simulation and showing the results.
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

def demo_basic_sir():
    """Demonstrate Basic SIR model features."""
    print("🦟 Basic SIR Model Demo")
    print("=" * 30)
    
    # Create model with different transmission rates
    scenarios = [
        {"name": "Low Transmission", "beta": 0.2, "color": "blue"},
        {"name": "Medium Transmission", "beta": 0.4, "color": "orange"},
        {"name": "High Transmission", "beta": 0.6, "color": "red"}
    ]
    
    plt.figure(figsize=(15, 10))
    
    for i, scenario in enumerate(scenarios):
        model = MalariaSIRModel(
            beta=scenario["beta"],
            gamma=0.05,
            mu=0.0001,
            N=10000,
            I0=50,
            R0=0
        )
        
        t, solution = model.simulate(t_span=(0, 1000), num_points=1000)
        S, I, R = solution.T
        
        # Plot infection rate
        plt.subplot(2, 2, 1)
        infection_rate = I / model.N * 100
        plt.plot(t, infection_rate, color=scenario["color"], 
                label=f"{scenario['name']} (R₀={model.R0_basic:.1f})", linewidth=2)
        
        # Plot cumulative cases
        plt.subplot(2, 2, 2)
        cumulative_cases = (R + I) / model.N * 100
        plt.plot(t, cumulative_cases, color=scenario["color"], 
                label=f"{scenario['name']}", linewidth=2)
        
        # Print metrics
        metrics = model.get_epidemic_metrics(solution)
        print(f"\n{scenario['name']}:")
        print(f"  R₀ = {metrics['R0_basic']:.2f}")
        print(f"  Peak infection rate = {metrics['peak_infection_rate']*100:.1f}%")
        print(f"  Final attack rate = {metrics['final_attack_rate']*100:.1f}%")
        print(f"  Time to peak = {metrics['time_to_peak']} days")
    
    # Format plots
    plt.subplot(2, 2, 1)
    plt.title('Infection Rate Over Time')
    plt.xlabel('Time (days)')
    plt.ylabel('Infection Rate (%)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(2, 2, 2)
    plt.title('Cumulative Cases Over Time')
    plt.xlabel('Time (days)')
    plt.ylabel('Cumulative Cases (%)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Parameter sensitivity
    plt.subplot(2, 2, 3)
    beta_values = np.linspace(0.1, 0.8, 20)
    r0_values = []
    peak_rates = []
    
    for beta in beta_values:
        temp_model = MalariaSIRModel(beta=beta, gamma=0.05, mu=0.0001, N=10000, I0=50, R0=0)
        t_temp, solution_temp = temp_model.simulate(t_span=(0, 1000), num_points=1000)
        metrics_temp = temp_model.get_epidemic_metrics(solution_temp)
        r0_values.append(metrics_temp['R0_basic'])
        peak_rates.append(metrics_temp['peak_infection_rate'] * 100)
    
    plt.plot(beta_values, r0_values, 'b-', linewidth=2, label='R₀')
    plt.xlabel('Transmission Rate (β)')
    plt.ylabel('Basic Reproduction Number (R₀)')
    plt.title('R₀ vs Transmission Rate')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(2, 2, 4)
    plt.plot(beta_values, peak_rates, 'r-', linewidth=2, label='Peak Infection Rate')
    plt.xlabel('Transmission Rate (β)')
    plt.ylabel('Peak Infection Rate (%)')
    plt.title('Peak Infection Rate vs Transmission Rate')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('results/plots/web_interface_demo.png', dpi=300, bbox_inches='tight')
    print(f"\n📊 Demo plot saved to: results/plots/web_interface_demo.png")
    plt.show()

def demo_model_comparison():
    """Demonstrate different model types."""
    print("\n🔄 Model Comparison Demo")
    print("=" * 30)
    
    models = {
        'Basic SIR': MalariaSIRModel(beta=0.4, gamma=0.05, N=10000, I0=50),
        'SEIR': MalariaSEIRModel(beta=0.4, sigma=0.1, gamma=0.05, N=10000, I0=30, E0=20),
        'Age-Structured': AgeStructuredMalariaSIR(N_total=10000, I0=np.array([5, 10, 15, 10, 5, 5])),
        'Vector Dynamics': MalariaVectorModel(N_h=10000, N_v=50000, I_h0=30, E_h0=20)
    }
    
    plt.figure(figsize=(15, 8))
    
    for i, (name, model) in enumerate(models.items()):
        t, solution = model.simulate(t_span=(0, 1000), num_points=1000)
        
        if hasattr(model, 'n_age_groups'):
            # Age-structured model
            I = np.sum(solution[:, model.n_age_groups:2*model.n_age_groups], axis=1)
            N = model.N_total
        elif hasattr(model, 'N_h'):
            # Vector dynamics model
            I = solution[:, 2]  # Infected humans
            N = model.N_h
        else:
            # Basic SIR or SEIR
            if solution.shape[1] == 4:  # SEIR
                I = solution[:, 2]  # Infected
            else:  # Basic SIR
                I = solution[:, 1]  # Infected
            N = model.N
        
        infection_rate = I / N * 100
        
        plt.subplot(2, 2, i+1)
        plt.plot(t, infection_rate, linewidth=2)
        plt.title(f'{name} Model')
        plt.xlabel('Time (days)')
        plt.ylabel('Infection Rate (%)')
        plt.grid(True, alpha=0.3)
        
        # Calculate and display key metrics
        peak_rate = np.max(infection_rate)
        time_to_peak = np.argmax(infection_rate)
        
        if hasattr(model, 'R0_basic'):
            r0 = model.R0_basic
        elif hasattr(model, 'R0'):
            r0 = model.R0
        else:
            r0 = "N/A"
        
        plt.text(0.05, 0.95, f'R₀ = {r0}\nPeak = {peak_rate:.1f}%\nTime = {time_to_peak} days', 
                transform=plt.gca().transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        print(f"\n{name}:")
        print(f"  R₀ = {r0}")
        print(f"  Peak infection rate = {peak_rate:.1f}%")
        print(f"  Time to peak = {time_to_peak} days")
    
    plt.tight_layout()
    plt.savefig('results/plots/model_comparison_demo.png', dpi=300, bbox_inches='tight')
    print(f"\n📊 Model comparison plot saved to: results/plots/model_comparison_demo.png")
    plt.show()

def main():
    """Run the web interface demo."""
    print("🌐 Malaria Model Simulation - Web Interface Demo")
    print("=" * 50)
    print()
    print("This demo shows the key features available in the web interface:")
    print("• Parameter adjustment and real-time results")
    print("• Multiple model types (SIR, SEIR, Age-Structured, Vector Dynamics)")
    print("• Interactive visualizations with Plotly")
    print("• Parameter sensitivity analysis")
    print("• Model comparison capabilities")
    print("• Intervention effectiveness analysis")
    print()
    print("🚀 To start the web interface, run:")
    print("   ./start_web_interface.sh")
    print()
    print("📖 For experimental setups, see:")
    print("   cat EXPERIMENTAL_SETUPS.md")
    print()
    
    # Run demos
    demo_basic_sir()
    demo_model_comparison()
    
    print("\n✅ Demo completed!")
    print("\n🌐 Next steps:")
    print("1. Start the web interface: ./start_web_interface.sh")
    print("2. Open your browser to: http://localhost:8501")
    print("3. Try the experimental setups in EXPERIMENTAL_SETUPS.md")
    print("4. Explore different model types and parameters")

if __name__ == "__main__":
    main()


