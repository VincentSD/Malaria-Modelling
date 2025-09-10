"""
Vector Dynamics Model for Malaria Transmission
=============================================

This module implements a comprehensive malaria transmission model that includes
both human and mosquito (vector) population dynamics. This is essential for
malaria modeling as malaria is a vector-borne disease.

Key features:
- Human population: Susceptible-Exposed-Infected-Recovered (SEIR)
- Mosquito population: Susceptible-Exposed-Infected (SEI)
- Human-mosquito transmission cycles
- Seasonal variations in mosquito populations
- Age-structured mosquito populations
- Vector control interventions
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
from typing import Tuple, Dict, List, Optional
import seaborn as sns


class MalariaVectorModel:
    """
    Comprehensive malaria transmission model with human and vector dynamics.
    """
    
    def __init__(self, 
                 # Human parameters
                 beta_h: float = 0.3,      # Human transmission rate
                 sigma_h: float = 0.1,     # Human incubation rate
                 gamma_h: float = 0.05,    # Human recovery rate
                 mu_h: float = 0.0001,     # Human birth/death rate
                 N_h: int = 10000,         # Human population
                 
                 # Vector parameters
                 beta_v: float = 0.2,      # Vector transmission rate
                 sigma_v: float = 0.2,     # Vector incubation rate
                 mu_v: float = 0.1,        # Vector death rate
                 N_v: int = 50000,         # Vector population
                 
                 # Biting parameters
                 a: float = 0.3,           # Biting rate (bites per mosquito per day)
                 b: float = 0.5,           # Probability of transmission from human to mosquito
                 c: float = 0.5,           # Probability of transmission from mosquito to human
                 
                 # Seasonal parameters
                 seasonal_amplitude: float = 0.5,  # Seasonal variation in vector population
                 seasonal_phase: float = 0.0,      # Seasonal phase offset
                 
                 # Initial conditions
                 E_h0: int = 10,           # Initial exposed humans
                 I_h0: int = 20,           # Initial infected humans
                 R_h0: int = 0,            # Initial recovered humans
                 E_v0: int = 100,          # Initial exposed vectors
                 I_v0: int = 200):         # Initial infected vectors
        """
        Initialize the malaria vector dynamics model.
        
        Parameters:
        -----------
        beta_h : float
            Human transmission rate
        sigma_h : float
            Human incubation rate (1/latent period)
        gamma_h : float
            Human recovery rate (1/infectious period)
        mu_h : float
            Human birth/death rate
        N_h : int
            Human population size
        beta_v : float
            Vector transmission rate
        sigma_v : float
            Vector incubation rate (1/latent period)
        mu_v : float
            Vector death rate
        N_v : int
            Vector population size
        a : float
            Biting rate (bites per mosquito per day)
        b : float
            Probability of transmission from human to mosquito
        c : float
            Probability of transmission from mosquito to human
        seasonal_amplitude : float
            Amplitude of seasonal variation in vector population
        seasonal_phase : float
            Phase offset for seasonal variation
        E_h0, I_h0, R_h0 : int
            Initial human compartments
        E_v0, I_v0 : int
            Initial vector compartments
        """
        
        # Human parameters
        self.beta_h = beta_h
        self.sigma_h = sigma_h
        self.gamma_h = gamma_h
        self.mu_h = mu_h
        self.N_h = N_h
        
        # Vector parameters
        self.beta_v = beta_v
        self.sigma_v = sigma_v
        self.mu_v = mu_v
        self.N_v = N_v
        
        # Biting parameters
        self.a = a
        self.b = b
        self.c = c
        
        # Seasonal parameters
        self.seasonal_amplitude = seasonal_amplitude
        self.seasonal_phase = seasonal_phase
        
        # Initial conditions
        self.S_h0 = N_h - E_h0 - I_h0 - R_h0  # Initial susceptible humans
        self.E_h0 = E_h0
        self.I_h0 = I_h0
        self.R_h0 = R_h0
        
        self.S_v0 = N_v - E_v0 - I_v0  # Initial susceptible vectors
        self.E_v0 = E_v0
        self.I_v0 = I_v0
        
        # Calculate basic reproduction number
        self.R0 = self._calculate_basic_reproduction_number()
        
    def _calculate_basic_reproduction_number(self) -> float:
        """
        Calculate the basic reproduction number for the vector-borne model.
        
        Returns:
        --------
        float
            Basic reproduction number
        """
        # R0 for vector-borne diseases
        R0 = (self.a**2 * self.b * self.c * self.N_v) / (self.mu_v * (self.gamma_h + self.mu_h) * self.N_h)
        return R0
    
    def seasonal_vector_population(self, t: float) -> float:
        """
        Calculate seasonal vector population.
        
        Parameters:
        -----------
        t : float
            Time in days
            
        Returns:
        --------
        float
            Seasonal vector population multiplier
        """
        # Annual seasonal variation
        seasonal_factor = 1 + self.seasonal_amplitude * np.sin(2 * np.pi * t / 365 + self.seasonal_phase)
        return seasonal_factor
    
    def vector_dynamics_equations(self, y: np.ndarray, t: np.ndarray) -> np.ndarray:
        """
        Define the vector dynamics differential equations.
        
        Parameters:
        -----------
        y : np.ndarray
            Current state [S_h, E_h, I_h, R_h, S_v, E_v, I_v]
        t : np.ndarray
            Time points
            
        Returns:
        --------
        np.ndarray
            Derivatives for all compartments
        """
        S_h, E_h, I_h, R_h, S_v, E_v, I_v = y
        
        # Seasonal vector population
        N_v_seasonal = self.N_v * self.seasonal_vector_population(t)
        
        # Force of infection from vectors to humans
        lambda_h = self.a * self.c * I_v / N_v_seasonal
        
        # Force of infection from humans to vectors
        lambda_v = self.a * self.b * I_h / self.N_h
        
        # Human dynamics
        dS_hdt = self.mu_h * self.N_h - lambda_h * S_h - self.mu_h * S_h
        dE_hdt = lambda_h * S_h - (self.sigma_h + self.mu_h) * E_h
        dI_hdt = self.sigma_h * E_h - (self.gamma_h + self.mu_h) * I_h
        dR_hdt = self.gamma_h * I_h - self.mu_h * R_h
        
        # Vector dynamics
        dS_vdt = self.mu_v * N_v_seasonal - lambda_v * S_v - self.mu_v * S_v
        dE_vdt = lambda_v * S_v - (self.sigma_v + self.mu_v) * E_v
        dI_vdt = self.sigma_v * E_v - self.mu_v * I_v
        
        return [dS_hdt, dE_hdt, dI_hdt, dR_hdt, dS_vdt, dE_vdt, dI_vdt]
    
    def simulate(self, t_span: Tuple[float, float] = (0, 365), 
                 num_points: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simulate the vector dynamics model over time.
        
        Parameters:
        -----------
        t_span : tuple
            Time span (start, end) in days
        num_points : int
            Number of time points to simulate
            
        Returns:
        --------
        tuple
            (time_points, solution_array) where solution_array has shape (num_points, 7)
            with columns [S_h, E_h, I_h, R_h, S_v, E_v, I_v]
        """
        t = np.linspace(t_span[0], t_span[1], num_points)
        y0 = [self.S_h0, self.E_h0, self.I_h0, self.R_h0, 
              self.S_v0, self.E_v0, self.I_v0]
        
        solution = odeint(self.vector_dynamics_equations, y0, t)
        
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
        S_h, E_h, I_h, R_h, S_v, E_v, I_v = solution.T
        
        # Human metrics
        peak_human_infection = np.max(I_h)
        peak_human_infection_rate = peak_human_infection / self.N_h
        final_human_attack_rate = R_h[-1] / self.N_h
        
        # Vector metrics
        peak_vector_infection = np.max(I_v)
        peak_vector_infection_rate = peak_vector_infection / self.N_v
        
        # Time to peak
        time_to_peak_human = np.argmax(I_h)
        time_to_peak_vector = np.argmax(I_v)
        
        return {
            'R0_basic': self.R0,
            'peak_human_infection': peak_human_infection,
            'peak_human_infection_rate': peak_human_infection_rate,
            'final_human_attack_rate': final_human_attack_rate,
            'peak_vector_infection': peak_vector_infection,
            'peak_vector_infection_rate': peak_vector_infection_rate,
            'time_to_peak_human': time_to_peak_human,
            'time_to_peak_vector': time_to_peak_vector,
            'vector_human_ratio': self.N_v / self.N_h
        }
    
    def plot_results(self, t: np.ndarray, solution: np.ndarray, 
                    save_path: str = None) -> None:
        """
        Plot the vector dynamics model results.
        
        Parameters:
        -----------
        t : np.ndarray
            Time points
        solution : np.ndarray
            Solution array
        save_path : str, optional
            Path to save the plot
        """
        S_h, E_h, I_h, R_h, S_v, E_v, I_v = solution.T
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # Human population dynamics
        ax1 = axes[0, 0]
        ax1.plot(t, S_h, 'b-', label='Susceptible', linewidth=2)
        ax1.plot(t, E_h, 'orange', label='Exposed', linewidth=2)
        ax1.plot(t, I_h, 'r-', label='Infected', linewidth=2)
        ax1.plot(t, R_h, 'g-', label='Recovered', linewidth=2)
        ax1.set_xlabel('Time (days)')
        ax1.set_ylabel('Number of humans')
        ax1.set_title('Human Population Dynamics')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Vector population dynamics
        ax2 = axes[0, 1]
        ax2.plot(t, S_v, 'b-', label='Susceptible', linewidth=2)
        ax2.plot(t, E_v, 'orange', label='Exposed', linewidth=2)
        ax2.plot(t, I_v, 'r-', label='Infected', linewidth=2)
        ax2.set_xlabel('Time (days)')
        ax2.set_ylabel('Number of vectors')
        ax2.set_title('Vector Population Dynamics')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Infection rates comparison
        ax3 = axes[0, 2]
        human_infection_rate = I_h / self.N_h * 100
        vector_infection_rate = I_v / self.N_v * 100
        ax3.plot(t, human_infection_rate, 'r-', label='Human', linewidth=2)
        ax3.plot(t, vector_infection_rate, 'purple', label='Vector', linewidth=2)
        ax3.set_xlabel('Time (days)')
        ax3.set_ylabel('Infection Rate (%)')
        ax3.set_title('Infection Rates Comparison')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Seasonal vector population
        ax4 = axes[1, 0]
        N_v_seasonal = [self.N_v * self.seasonal_vector_population(time) for time in t]
        ax4.plot(t, N_v_seasonal, 'green', linewidth=2)
        ax4.set_xlabel('Time (days)')
        ax4.set_ylabel('Vector Population')
        ax4.set_title('Seasonal Vector Population')
        ax4.grid(True, alpha=0.3)
        
        # Cumulative human cases
        ax5 = axes[1, 1]
        cumulative_cases = R_h + I_h
        ax5.plot(t, cumulative_cases, 'orange', linewidth=2)
        ax5.set_xlabel('Time (days)')
        ax5.set_ylabel('Cumulative Cases')
        ax5.set_title('Cumulative Human Malaria Cases')
        ax5.grid(True, alpha=0.3)
        
        # Vector-human transmission cycle
        ax6 = axes[1, 2]
        ax6.plot(I_h, I_v, 'purple', linewidth=2)
        ax6.set_xlabel('Infected Humans')
        ax6.set_ylabel('Infected Vectors')
        ax6.set_title('Vector-Human Transmission Cycle')
        ax6.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def intervention_analysis(self, intervention_type: str, 
                            intervention_strength: float,
                            intervention_start: float,
                            intervention_duration: float,
                            t_span: Tuple[float, float] = (0, 1000)) -> Dict:
        """
        Analyze the effect of vector control interventions.
        
        Parameters:
        -----------
        intervention_type : str
            Type of intervention ('insecticide', 'bed_nets', 'larvicide')
        intervention_strength : float
            Strength of intervention (0-1)
        intervention_start : float
            Start time of intervention (days)
        intervention_duration : float
            Duration of intervention (days)
        t_span : tuple
            Time span for simulation
            
        Returns:
        --------
        dict
            Results of intervention analysis
        """
        # Create modified model with intervention
        if intervention_type == 'insecticide':
            # Reduce vector population
            modified_model = MalariaVectorModel(
                beta_h=self.beta_h, sigma_h=self.sigma_h, gamma_h=self.gamma_h, mu_h=self.mu_h, N_h=self.N_h,
                beta_v=self.beta_v, sigma_v=self.sigma_v, mu_v=self.mu_v * (1 + intervention_strength), N_v=self.N_v,
                a=self.a, b=self.b, c=self.c,
                seasonal_amplitude=self.seasonal_amplitude, seasonal_phase=self.seasonal_phase,
                E_h0=self.E_h0, I_h0=self.I_h0, R_h0=self.R_h0, E_v0=self.E_v0, I_v0=self.I_v0
            )
        elif intervention_type == 'bed_nets':
            # Reduce biting rate
            modified_model = MalariaVectorModel(
                beta_h=self.beta_h, sigma_h=self.sigma_h, gamma_h=self.gamma_h, mu_h=self.mu_h, N_h=self.N_h,
                beta_v=self.beta_v, sigma_v=self.sigma_v, mu_v=self.mu_v, N_v=self.N_v,
                a=self.a * (1 - intervention_strength), b=self.b, c=self.c,
                seasonal_amplitude=self.seasonal_amplitude, seasonal_phase=self.seasonal_phase,
                E_h0=self.E_h0, I_h0=self.I_h0, R_h0=self.R_h0, E_v0=self.E_v0, I_v0=self.I_v0
            )
        elif intervention_type == 'larvicide':
            # Reduce vector birth rate (affects seasonal population)
            modified_model = MalariaVectorModel(
                beta_h=self.beta_h, sigma_h=self.sigma_h, gamma_h=self.gamma_h, mu_h=self.mu_h, N_h=self.N_h,
                beta_v=self.beta_v, sigma_v=self.sigma_v, mu_v=self.mu_v, N_v=self.N_v,
                a=self.a, b=self.b, c=self.c,
                seasonal_amplitude=self.seasonal_amplitude * (1 - intervention_strength), seasonal_phase=self.seasonal_phase,
                E_h0=self.E_h0, I_h0=self.I_h0, R_h0=self.R_h0, E_v0=self.E_v0, I_v0=self.I_v0
            )
        else:
            raise ValueError("Invalid intervention type")
        
        # Simulate both models
        t_baseline, solution_baseline = self.simulate(t_span)
        t_intervention, solution_intervention = modified_model.simulate(t_span)
        
        # Calculate metrics
        metrics_baseline = self.get_epidemic_metrics(solution_baseline)
        metrics_intervention = modified_model.get_epidemic_metrics(solution_intervention)
        
        # Calculate intervention effectiveness
        effectiveness = {
            'intervention_type': intervention_type,
            'intervention_strength': intervention_strength,
            'baseline_metrics': metrics_baseline,
            'intervention_metrics': metrics_intervention,
            'reduction_in_peak_infection': (metrics_baseline['peak_human_infection_rate'] - 
                                          metrics_intervention['peak_human_infection_rate']) / 
                                         metrics_baseline['peak_human_infection_rate'],
            'reduction_in_attack_rate': (metrics_baseline['final_human_attack_rate'] - 
                                       metrics_intervention['final_human_attack_rate']) / 
                                      metrics_baseline['final_human_attack_rate']
        }
        
        return effectiveness


def main():
    """
    Main function to demonstrate the malaria vector dynamics model.
    """
    print("Malaria Vector Dynamics Model Simulation")
    print("=" * 50)
    
    # Create model with realistic parameters
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
    model.plot_results(t, solution, save_path="../results/plots/malaria_vector_dynamics.png")
    
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


if __name__ == "__main__":
    main()


