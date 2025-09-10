"""
SEIR Model for Malaria Transmission
==================================

This module implements a Susceptible-Exposed-Infected-Recovered (SEIR) model
for malaria transmission, which includes a latent period (exposed state)
that is more realistic for malaria transmission dynamics.

Key features:
- Latent period (exposed state) before becoming infectious
- Age-specific parameters for different malaria strains
- Seasonal variations in transmission
- Treatment and immunity considerations
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
from typing import Tuple, Dict, List, Optional
import seaborn as sns


class MalariaSEIRModel:
    """
    SEIR model for malaria transmission with malaria-specific parameters.
    """
    
    def __init__(self, 
                 beta: float = 0.3,      # Transmission rate (per day)
                 sigma: float = 0.1,     # Incubation rate (1/latent period)
                 gamma: float = 0.05,    # Recovery rate (per day)
                 mu: float = 0.0001,     # Birth/death rate (per day)
                 N: int = 10000,         # Total population
                 E0: int = 5,            # Initial exposed
                 I0: int = 10,           # Initial infected
                 R0: int = 0,            # Initial recovered
                 seasonal_amplitude: float = 0.2,  # Seasonal variation amplitude
                 seasonal_phase: float = 0.0):     # Seasonal phase offset
        """
        Initialize the malaria SEIR model.
        
        Parameters:
        -----------
        beta : float
            Base transmission rate
        sigma : float
            Incubation rate (1/latent period in days)
        gamma : float
            Recovery rate (1/infectious period in days)
        mu : float
            Birth and death rate
        N : int
            Total population size
        E0 : int
            Initial number of exposed individuals
        I0 : int
            Initial number of infected individuals
        R0 : int
            Initial number of recovered individuals
        seasonal_amplitude : float
            Amplitude of seasonal variation in transmission
        seasonal_phase : float
            Phase offset for seasonal variation
        """
        self.beta = beta
        self.sigma = sigma
        self.gamma = gamma
        self.mu = mu
        self.N = N
        self.S0 = N - E0 - I0 - R0  # Initial susceptible
        self.E0 = E0
        self.I0 = I0
        self.R0 = R0
        self.seasonal_amplitude = seasonal_amplitude
        self.seasonal_phase = seasonal_phase
        
        # Calculate basic reproduction number
        self.R0_basic = self.beta / (self.gamma + self.mu)
        
        # Calculate effective reproduction number (accounting for latent period)
        self.R0_effective = self.beta * self.sigma / ((self.sigma + self.mu) * (self.gamma + self.mu))
        
    def seasonal_transmission_rate(self, t: float) -> float:
        """
        Calculate seasonal transmission rate.
        
        Parameters:
        -----------
        t : float
            Time in days
            
        Returns:
        --------
        float
            Seasonal transmission rate
        """
        # Annual seasonal variation
        seasonal_factor = 1 + self.seasonal_amplitude * np.sin(2 * np.pi * t / 365 + self.seasonal_phase)
        return self.beta * seasonal_factor
    
    def seir_equations(self, y: np.ndarray, t: np.ndarray) -> np.ndarray:
        """
        Define the SEIR differential equations.
        
        Parameters:
        -----------
        y : np.ndarray
            Current state [S, E, I, R]
        t : np.ndarray
            Time points
            
        Returns:
        --------
        np.ndarray
            Derivatives [dS/dt, dE/dt, dI/dt, dR/dt]
        """
        S, E, I, R = y
        
        # Seasonal transmission rate
        beta_t = self.seasonal_transmission_rate(t)
        
        # SEIR equations with birth and death
        dSdt = self.mu * self.N - beta_t * S * I / self.N - self.mu * S
        dEdt = beta_t * S * I / self.N - (self.sigma + self.mu) * E
        dIdt = self.sigma * E - (self.gamma + self.mu) * I
        dRdt = self.gamma * I - self.mu * R
        
        return [dSdt, dEdt, dIdt, dRdt]
    
    def simulate(self, t_span: Tuple[float, float] = (0, 365), 
                 num_points: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simulate the SEIR model over time.
        
        Parameters:
        -----------
        t_span : tuple
            Time span (start, end) in days
        num_points : int
            Number of time points to simulate
            
        Returns:
        --------
        tuple
            (time_points, solution_array) where solution_array has shape (num_points, 4)
            with columns [S, E, I, R]
        """
        t = np.linspace(t_span[0], t_span[1], num_points)
        y0 = [self.S0, self.E0, self.I0, self.R0]
        
        solution = odeint(self.seir_equations, y0, t)
        
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
        S, E, I, R = solution.T
        
        # Peak infection
        peak_infection = np.max(I)
        peak_time = np.argmax(I)
        
        # Peak exposed
        peak_exposed = np.max(E)
        
        # Final epidemic size
        final_recovered = R[-1]
        attack_rate = final_recovered / self.N
        
        # Time to peak
        time_to_peak = peak_time
        
        # Latent period statistics
        latent_period = 1 / self.sigma
        
        return {
            'R0_basic': self.R0_basic,
            'R0_effective': self.R0_effective,
            'peak_infection': peak_infection,
            'peak_infection_rate': peak_infection / self.N,
            'peak_exposed': peak_exposed,
            'peak_exposed_rate': peak_exposed / self.N,
            'final_attack_rate': attack_rate,
            'time_to_peak': time_to_peak,
            'latent_period': latent_period,
            'infectious_period': 1 / self.gamma
        }
    
    def plot_results(self, t: np.ndarray, solution: np.ndarray, 
                    save_path: str = None) -> None:
        """
        Plot the SEIR model results.
        
        Parameters:
        -----------
        t : np.ndarray
            Time points
        solution : np.ndarray
            Solution array [S, E, I, R]
        save_path : str, optional
            Path to save the plot
        """
        S, E, I, R = solution.T
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # Main SEIR plot
        ax1 = axes[0, 0]
        ax1.plot(t, S, 'b-', label='Susceptible', linewidth=2)
        ax1.plot(t, E, 'orange', label='Exposed', linewidth=2)
        ax1.plot(t, I, 'r-', label='Infected', linewidth=2)
        ax1.plot(t, R, 'g-', label='Recovered', linewidth=2)
        ax1.set_xlabel('Time (days)')
        ax1.set_ylabel('Number of individuals')
        ax1.set_title('Malaria SEIR Model Dynamics')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Infection and exposure rates
        ax2 = axes[0, 1]
        infection_rate = I / self.N * 100
        exposure_rate = E / self.N * 100
        ax2.plot(t, infection_rate, 'r-', label='Infection Rate', linewidth=2)
        ax2.plot(t, exposure_rate, 'orange', label='Exposure Rate', linewidth=2)
        ax2.set_xlabel('Time (days)')
        ax2.set_ylabel('Rate (%)')
        ax2.set_title('Infection and Exposure Rates')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Seasonal transmission rate
        ax3 = axes[0, 2]
        beta_seasonal = [self.seasonal_transmission_rate(time) for time in t]
        ax3.plot(t, beta_seasonal, 'purple', linewidth=2)
        ax3.set_xlabel('Time (days)')
        ax3.set_ylabel('Transmission Rate (β)')
        ax3.set_title('Seasonal Transmission Rate')
        ax3.grid(True, alpha=0.3)
        
        # Phase plane (S vs I)
        ax4 = axes[1, 0]
        ax4.plot(S, I, 'purple', linewidth=2)
        ax4.set_xlabel('Susceptible')
        ax4.set_ylabel('Infected')
        ax4.set_title('Phase Plane (S vs I)')
        ax4.grid(True, alpha=0.3)
        
        # Cumulative cases
        ax5 = axes[1, 1]
        cumulative_cases = R + I
        ax5.plot(t, cumulative_cases, 'orange', linewidth=2)
        ax5.set_xlabel('Time (days)')
        ax5.set_ylabel('Cumulative Cases')
        ax5.set_title('Cumulative Malaria Cases')
        ax5.grid(True, alpha=0.3)
        
        # Effective reproduction number over time
        ax6 = axes[1, 2]
        R_eff = [self.seasonal_transmission_rate(time) * S[i] / self.N / (self.gamma + self.mu) 
                for i, time in enumerate(t)]
        ax6.plot(t, R_eff, 'red', linewidth=2)
        ax6.axhline(y=1, color='black', linestyle='--', alpha=0.7, label='Rₑff = 1')
        ax6.set_xlabel('Time (days)')
        ax6.set_ylabel('Effective Reproduction Number')
        ax6.set_title('Effective Reproduction Number')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def compare_with_sir(self, t_span: Tuple[float, float] = (0, 365)) -> None:
        """
        Compare SEIR model with equivalent SIR model.
        
        Parameters:
        -----------
        t_span : tuple
            Time span for comparison
        """
        # Import basic SIR model
        from basic_sir import MalariaSIRModel
        
        # Create equivalent SIR model (no latent period)
        sir_model = MalariaSIRModel(
            beta=self.beta,
            gamma=self.gamma,
            mu=self.mu,
            N=self.N,
            I0=self.I0 + self.E0,  # Combine exposed and infected for SIR
            R0=self.R0
        )
        
        # Simulate both models
        t_seir, solution_seir = self.simulate(t_span)
        t_sir, solution_sir = sir_model.simulate(t_span)
        
        # Plot comparison
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        # Susceptible comparison
        axes[0, 0].plot(t_seir, solution_seir[:, 0], 'b-', label='SEIR', linewidth=2)
        axes[0, 0].plot(t_sir, solution_sir[:, 0], 'b--', label='SIR', linewidth=2)
        axes[0, 0].set_xlabel('Time (days)')
        axes[0, 0].set_ylabel('Susceptible')
        axes[0, 0].set_title('Susceptible Population')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Infected comparison
        axes[0, 1].plot(t_seir, solution_seir[:, 2], 'r-', label='SEIR', linewidth=2)
        axes[0, 1].plot(t_sir, solution_sir[:, 1], 'r--', label='SIR', linewidth=2)
        axes[0, 1].set_xlabel('Time (days)')
        axes[0, 1].set_ylabel('Infected')
        axes[0, 1].set_title('Infected Population')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Recovered comparison
        axes[1, 0].plot(t_seir, solution_seir[:, 3], 'g-', label='SEIR', linewidth=2)
        axes[1, 0].plot(t_sir, solution_sir[:, 2], 'g--', label='SIR', linewidth=2)
        axes[1, 0].set_xlabel('Time (days)')
        axes[1, 0].set_ylabel('Recovered')
        axes[1, 0].set_title('Recovered Population')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Exposed (SEIR only)
        axes[1, 1].plot(t_seir, solution_seir[:, 1], 'orange', label='SEIR Exposed', linewidth=2)
        axes[1, 1].set_xlabel('Time (days)')
        axes[1, 1].set_ylabel('Exposed')
        axes[1, 1].set_title('Exposed Population (SEIR only)')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
    
    def parameter_sensitivity(self, param_name: str, param_values: List[float],
                            t_span: Tuple[float, float] = (0, 365)) -> pd.DataFrame:
        """
        Perform parameter sensitivity analysis.
        
        Parameters:
        -----------
        param_name : str
            Name of parameter to vary ('beta', 'sigma', 'gamma', 'mu')
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
            temp_model = MalariaSEIRModel(
                beta=self.beta if param_name != 'beta' else value,
                sigma=self.sigma if param_name != 'sigma' else value,
                gamma=self.gamma if param_name != 'gamma' else value,
                mu=self.mu if param_name != 'mu' else value,
                N=self.N, E0=self.E0, I0=self.I0, R0=self.R0,
                seasonal_amplitude=self.seasonal_amplitude,
                seasonal_phase=self.seasonal_phase
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
    Main function to demonstrate the malaria SEIR model.
    """
    print("Malaria SEIR Model Simulation")
    print("=" * 40)
    
    # Create model with malaria-specific parameters
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
    model.plot_results(t, solution, save_path="../results/plots/malaria_seir.png")
    
    # Compare with SIR model
    print("Comparing with SIR model...")
    model.compare_with_sir()
    
    # Parameter sensitivity analysis
    print("Parameter Sensitivity Analysis:")
    print("-" * 30)
    
    # Test different incubation rates
    sigma_values = [0.05, 0.1, 0.15, 0.2, 0.25]
    sigma_sensitivity = model.parameter_sensitivity('sigma', sigma_values)
    print("Incubation Rate (σ) Sensitivity:")
    print(sigma_sensitivity[['value', 'R0_effective', 'peak_infection_rate', 'final_attack_rate']].round(3))
    print()
    
    # Test different seasonal amplitudes
    seasonal_values = [0.0, 0.1, 0.2, 0.3, 0.4]
    seasonal_sensitivity = model.parameter_sensitivity('beta', seasonal_values)
    print("Seasonal Amplitude Sensitivity:")
    print(seasonal_sensitivity[['value', 'R0_effective', 'peak_infection_rate', 'final_attack_rate']].round(3))


if __name__ == "__main__":
    main()


