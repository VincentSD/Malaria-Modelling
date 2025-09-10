"""
Basic SIR Model for Malaria Transmission
========================================

This module implements a basic Susceptible-Infected-Recovered (SIR) model
specifically adapted for malaria transmission dynamics.

Key malaria-specific considerations:
- Recovery rate accounts for natural immunity and treatment
- Transmission rate includes vector-borne transmission
- Basic reproduction number (R0) calculation
- Population dynamics with birth and death rates
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
from typing import Tuple, Dict, List


class MalariaSIRModel:
    """
    Basic SIR model for malaria transmission with malaria-specific parameters.
    """
    
    def __init__(self, 
                 beta: float = 0.3,      # Transmission rate (per day)
                 gamma: float = 0.1,     # Recovery rate (per day)
                 mu: float = 0.0001,     # Birth/death rate (per day)
                 N: int = 10000,         # Total population
                 I0: int = 10,           # Initial infected
                 R0: int = 0):           # Initial recovered
        """
        Initialize the malaria SIR model.
        
        Parameters:
        -----------
        beta : float
            Transmission rate (contacts per unit time * probability of transmission)
        gamma : float
            Recovery rate (1/infectious period)
        mu : float
            Birth and death rate (assumed equal for simplicity)
        N : int
            Total population size
        I0 : int
            Initial number of infected individuals
        R0 : int
            Initial number of recovered individuals
        """
        self.beta = beta
        self.gamma = gamma
        self.mu = mu
        self.N = N
        self.S0 = N - I0 - R0  # Initial susceptible
        self.I0 = I0
        self.R0 = R0
        
        # Calculate basic reproduction number
        self.R0_basic = self.beta / (self.gamma + self.mu)
        
    def sir_equations(self, y: np.ndarray, t: np.ndarray) -> np.ndarray:
        """
        Define the SIR differential equations.
        
        Parameters:
        -----------
        y : np.ndarray
            Current state [S, I, R]
        t : np.ndarray
            Time points
            
        Returns:
        --------
        np.ndarray
            Derivatives [dS/dt, dI/dt, dR/dt]
        """
        S, I, R = y
        
        # SIR equations with birth and death
        dSdt = self.mu * self.N - self.beta * S * I / self.N - self.mu * S
        dIdt = self.beta * S * I / self.N - (self.gamma + self.mu) * I
        dRdt = self.gamma * I - self.mu * R
        
        return [dSdt, dIdt, dRdt]
    
    def simulate(self, t_span: Tuple[float, float] = (0, 365), 
                 num_points: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simulate the SIR model over time.
        
        Parameters:
        -----------
        t_span : tuple
            Time span (start, end) in days
        num_points : int
            Number of time points to simulate
            
        Returns:
        --------
        tuple
            (time_points, solution_array) where solution_array has shape (num_points, 3)
            with columns [S, I, R]
        """
        t = np.linspace(t_span[0], t_span[1], num_points)
        y0 = [self.S0, self.I0, self.R0]
        
        solution = odeint(self.sir_equations, y0, t)
        
        return t, solution
    
    def get_epidemic_metrics(self, solution: np.ndarray) -> Dict[str, float]:
        """
        Calculate key epidemic metrics from simulation results.
        
        Parameters:
        -----------
        solution : np.ndarray
            Solution array from simulation
            
        Returns:
        --------
        dict
            Dictionary containing epidemic metrics
        """
        S, I, R = solution.T
        
        # Peak infection
        peak_infection = np.max(I)
        peak_time = np.argmax(I)
        
        # Final epidemic size
        final_recovered = R[-1]
        attack_rate = final_recovered / self.N
        
        # Time to peak
        time_to_peak = peak_time
        
        return {
            'R0_basic': self.R0_basic,
            'peak_infection': peak_infection,
            'peak_infection_rate': peak_infection / self.N,
            'final_attack_rate': attack_rate,
            'time_to_peak': time_to_peak
        }
    
    def plot_results(self, t: np.ndarray, solution: np.ndarray, 
                    save_path: str = None) -> None:
        """
        Plot the SIR model results.
        
        Parameters:
        -----------
        t : np.ndarray
            Time points
        solution : np.ndarray
            Solution array [S, I, R]
        save_path : str, optional
            Path to save the plot
        """
        S, I, R = solution.T
        
        plt.figure(figsize=(12, 8))
        
        # Main SIR plot
        plt.subplot(2, 2, 1)
        plt.plot(t, S, 'b-', label='Susceptible', linewidth=2)
        plt.plot(t, I, 'r-', label='Infected', linewidth=2)
        plt.plot(t, R, 'g-', label='Recovered', linewidth=2)
        plt.xlabel('Time (days)')
        plt.ylabel('Number of individuals')
        plt.title('Malaria SIR Model Dynamics')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Infection rate
        plt.subplot(2, 2, 2)
        infection_rate = I / self.N * 100
        plt.plot(t, infection_rate, 'r-', linewidth=2)
        plt.xlabel('Time (days)')
        plt.ylabel('Infection Rate (%)')
        plt.title('Infection Rate Over Time')
        plt.grid(True, alpha=0.3)
        
        # Phase plane
        plt.subplot(2, 2, 3)
        plt.plot(S, I, 'purple', linewidth=2)
        plt.xlabel('Susceptible')
        plt.ylabel('Infected')
        plt.title('Phase Plane (S vs I)')
        plt.grid(True, alpha=0.3)
        
        # Cumulative cases
        plt.subplot(2, 2, 4)
        cumulative_cases = R + I
        plt.plot(t, cumulative_cases, 'orange', linewidth=2)
        plt.xlabel('Time (days)')
        plt.ylabel('Cumulative Cases')
        plt.title('Cumulative Malaria Cases')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def parameter_sensitivity(self, param_name: str, param_values: List[float],
                            t_span: Tuple[float, float] = (0, 365)) -> pd.DataFrame:
        """
        Perform parameter sensitivity analysis.
        
        Parameters:
        -----------
        param_name : str
            Name of parameter to vary ('beta', 'gamma', 'mu')
        param_values : list
            List of parameter values to test
        t_span : tuple
            Time span for simulation
            
        Returns:
        --------
        pd.DataFrame
            DataFrame with results for each parameter value
        """
        results = []
        
        for value in param_values:
            # Create temporary model with modified parameter
            temp_model = MalariaSIRModel(
                beta=self.beta if param_name != 'beta' else value,
                gamma=self.gamma if param_name != 'gamma' else value,
                mu=self.mu if param_name != 'mu' else value,
                N=self.N, I0=self.I0, R0=self.R0
            )
            
            t, solution = temp_model.simulate(t_span)
            metrics = temp_model.get_epidemic_metrics(solution)
            
            results.append({
                'parameter': param_name,
                'value': value,
                **metrics
            })
        
        return pd.DataFrame(results)


def main():
    """
    Main function to demonstrate the malaria SIR model.
    """
    print("Malaria SIR Model Simulation")
    print("=" * 40)
    
    # Create model with malaria-specific parameters
    model = MalariaSIRModel(
        beta=0.4,      # Higher transmission rate for malaria
        gamma=0.05,    # Slower recovery (20 days average)
        mu=0.0001,     # Low birth/death rate
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
    model.plot_results(t, solution, save_path="../results/plots/malaria_sir_basic.png")
    
    # Parameter sensitivity analysis
    print("Parameter Sensitivity Analysis:")
    print("-" * 30)
    
    # Test different transmission rates
    beta_values = [0.2, 0.3, 0.4, 0.5, 0.6]
    beta_sensitivity = model.parameter_sensitivity('beta', beta_values)
    print("Transmission Rate (β) Sensitivity:")
    print(beta_sensitivity[['value', 'R0_basic', 'peak_infection_rate', 'final_attack_rate']].round(3))
    print()
    
    # Test different recovery rates
    gamma_values = [0.02, 0.05, 0.1, 0.15, 0.2]
    gamma_sensitivity = model.parameter_sensitivity('gamma', gamma_values)
    print("Recovery Rate (γ) Sensitivity:")
    print(gamma_sensitivity[['value', 'R0_basic', 'peak_infection_rate', 'final_attack_rate']].round(3))


if __name__ == "__main__":
    main()


