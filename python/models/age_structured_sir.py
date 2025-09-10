"""
Age-Structured SIR Model for Malaria Transmission
================================================

This module implements an age-structured SIR model for malaria transmission,
where the population is divided into different age groups with age-specific
transmission rates, recovery rates, and contact patterns.

Key features:
- Multiple age groups with different epidemiological parameters
- Age-specific contact matrices
- Age-dependent susceptibility and infectiousness
- Realistic population age distribution
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
from typing import Tuple, Dict, List, Optional
import seaborn as sns


class AgeStructuredMalariaSIR:
    """
    Age-structured SIR model for malaria transmission.
    """
    
    def __init__(self, 
                 age_groups: List[str] = None,
                 age_distribution: np.ndarray = None,
                 contact_matrix: np.ndarray = None,
                 beta_matrix: np.ndarray = None,
                 gamma: np.ndarray = None,
                 mu: np.ndarray = None,
                 N_total: int = 100000,
                 I0: np.ndarray = None):
        """
        Initialize the age-structured malaria SIR model.
        
        Parameters:
        -----------
        age_groups : list
            List of age group names
        age_distribution : np.ndarray
            Proportion of population in each age group
        contact_matrix : np.ndarray
            Contact rate matrix between age groups
        beta_matrix : np.ndarray
            Transmission rate matrix (age-specific)
        gamma : np.ndarray
            Recovery rate for each age group
        mu : np.ndarray
            Birth/death rate for each age group
        N_total : int
            Total population size
        I0 : np.ndarray
            Initial infected in each age group
        """
        
        # Default age groups (WHO standard)
        if age_groups is None:
            self.age_groups = ['0-4', '5-14', '15-29', '30-44', '45-59', '60+']
        else:
            self.age_groups = age_groups
            
        self.n_age_groups = len(self.age_groups)
        self.N_total = N_total
        
        # Default age distribution (typical for malaria-endemic regions)
        if age_distribution is None:
            self.age_distribution = np.array([0.15, 0.25, 0.25, 0.15, 0.10, 0.10])
        else:
            self.age_distribution = age_distribution
            
        # Ensure age distribution sums to 1
        self.age_distribution = self.age_distribution / np.sum(self.age_distribution)
        
        # Population in each age group
        self.N = self.age_distribution * N_total
        
        # Default contact matrix (symmetric, higher within-group contact)
        if contact_matrix is None:
            self.contact_matrix = self._create_default_contact_matrix()
        else:
            self.contact_matrix = contact_matrix
            
        # Default transmission rates (age-specific)
        if beta_matrix is None:
            self.beta_matrix = self._create_default_beta_matrix()
        else:
            self.beta_matrix = beta_matrix
            
        # Default recovery rates (age-specific)
        if gamma is None:
            self.gamma = np.array([0.08, 0.06, 0.05, 0.05, 0.04, 0.03])  # Slower recovery in older ages
        else:
            self.gamma = gamma
            
        # Default birth/death rates (age-specific)
        if mu is None:
            self.mu = np.array([0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001])
        else:
            self.mu = mu
            
        # Initial conditions
        if I0 is None:
            self.I0 = np.array([5, 10, 15, 10, 5, 5])  # Initial infected by age
        else:
            self.I0 = I0
            
        self.S0 = self.N - self.I0  # Initial susceptible
        self.R0 = np.zeros(self.n_age_groups)  # Initial recovered
        
        # Calculate age-specific R0
        self.R0_age = self._calculate_age_specific_r0()
        
    def _create_default_contact_matrix(self) -> np.ndarray:
        """Create a default contact matrix with realistic patterns."""
        # Higher contact within age groups, moderate between adjacent groups
        contact_matrix = np.ones((self.n_age_groups, self.n_age_groups))
        
        # Within-group contact is highest
        np.fill_diagonal(contact_matrix, 2.0)
        
        # Adjacent age groups have moderate contact
        for i in range(self.n_age_groups - 1):
            contact_matrix[i, i+1] = 1.5
            contact_matrix[i+1, i] = 1.5
            
        # Non-adjacent groups have lower contact
        for i in range(self.n_age_groups):
            for j in range(self.n_age_groups):
                if abs(i - j) > 1:
                    contact_matrix[i, j] = 0.5
                    
        return contact_matrix
    
    def _create_default_beta_matrix(self) -> np.ndarray:
        """Create age-specific transmission rate matrix."""
        # Base transmission rate
        beta_base = 0.3
        
        # Age-specific susceptibility and infectiousness
        susceptibility = np.array([1.2, 1.0, 0.8, 0.7, 0.6, 0.5])  # Higher in children
        infectiousness = np.array([1.1, 1.0, 0.9, 0.8, 0.7, 0.6])  # Higher in children
        
        # Create transmission matrix
        beta_matrix = np.outer(susceptibility, infectiousness) * beta_base
        
        return beta_matrix
    
    def _calculate_age_specific_r0(self) -> np.ndarray:
        """Calculate age-specific basic reproduction numbers."""
        R0_age = np.zeros(self.n_age_groups)
        
        for i in range(self.n_age_groups):
            # R0 for age group i
            R0_age[i] = (self.beta_matrix[i, i] * self.contact_matrix[i, i]) / (self.gamma[i] + self.mu[i])
            
        return R0_age
    
    def age_structured_sir_equations(self, y: np.ndarray, t: np.ndarray) -> np.ndarray:
        """
        Define the age-structured SIR differential equations.
        
        Parameters:
        -----------
        y : np.ndarray
            Current state [S1, S2, ..., Sn, I1, I2, ..., In, R1, R2, ..., Rn]
        t : np.ndarray
            Time points
            
        Returns:
        --------
        np.ndarray
            Derivatives for all compartments
        """
        # Reshape the state vector
        S = y[:self.n_age_groups]
        I = y[self.n_age_groups:2*self.n_age_groups]
        R = y[2*self.n_age_groups:]
        
        # Initialize derivatives
        dSdt = np.zeros(self.n_age_groups)
        dIdt = np.zeros(self.n_age_groups)
        dRdt = np.zeros(self.n_age_groups)
        
        # Calculate force of infection for each age group
        force_of_infection = np.zeros(self.n_age_groups)
        
        for i in range(self.n_age_groups):
            for j in range(self.n_age_groups):
                # Force of infection from age group j to age group i
                force_of_infection[i] += (self.beta_matrix[i, j] * 
                                        self.contact_matrix[i, j] * 
                                        I[j] / self.N[j])
        
        # SIR equations for each age group
        for i in range(self.n_age_groups):
            # Births enter susceptible compartment
            births = self.mu[i] * self.N[i]
            
            # Susceptible dynamics
            dSdt[i] = births - force_of_infection[i] * S[i] - self.mu[i] * S[i]
            
            # Infected dynamics
            dIdt[i] = force_of_infection[i] * S[i] - (self.gamma[i] + self.mu[i]) * I[i]
            
            # Recovered dynamics
            dRdt[i] = self.gamma[i] * I[i] - self.mu[i] * R[i]
        
        return np.concatenate([dSdt, dIdt, dRdt])
    
    def simulate(self, t_span: Tuple[float, float] = (0, 365), 
                 num_points: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simulate the age-structured SIR model over time.
        
        Parameters:
        -----------
        t_span : tuple
            Time span (start, end) in days
        num_points : int
            Number of time points to simulate
            
        Returns:
        --------
        tuple
            (time_points, solution_array) where solution_array has shape (num_points, 3*n_age_groups)
        """
        t = np.linspace(t_span[0], t_span[1], num_points)
        y0 = np.concatenate([self.S0, self.I0, self.R0])
        
        solution = odeint(self.age_structured_sir_equations, y0, t)
        
        return t, solution
    
    def get_age_specific_metrics(self, solution: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Calculate age-specific epidemic metrics.
        
        Parameters:
        -----------
        solution : np.ndarray
            Solution array from simulation
            
        Returns:
        --------
        dict
            Dictionary containing age-specific metrics
        """
        # Extract compartments
        S = solution[:, :self.n_age_groups]
        I = solution[:, self.n_age_groups:2*self.n_age_groups]
        R = solution[:, 2*self.n_age_groups:]
        
        metrics = {}
        
        # Peak infection by age
        metrics['peak_infection'] = np.max(I, axis=0)
        metrics['peak_infection_rate'] = metrics['peak_infection'] / self.N
        
        # Final attack rate by age
        metrics['final_attack_rate'] = R[-1, :] / self.N
        
        # Time to peak by age
        metrics['time_to_peak'] = np.argmax(I, axis=0)
        
        # Age-specific R0
        metrics['R0_age'] = self.R0_age
        
        return metrics
    
    def plot_age_structured_results(self, t: np.ndarray, solution: np.ndarray,
                                  save_path: str = None) -> None:
        """
        Plot the age-structured SIR model results.
        
        Parameters:
        -----------
        t : np.ndarray
            Time points
        solution : np.ndarray
            Solution array
        save_path : str, optional
            Path to save the plot
        """
        # Extract compartments
        S = solution[:, :self.n_age_groups]
        I = solution[:, self.n_age_groups:2*self.n_age_groups]
        R = solution[:, 2*self.n_age_groups:]
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # Plot 1: Total population dynamics
        ax1 = axes[0, 0]
        total_S = np.sum(S, axis=1)
        total_I = np.sum(I, axis=1)
        total_R = np.sum(R, axis=1)
        
        ax1.plot(t, total_S, 'b-', label='Susceptible', linewidth=2)
        ax1.plot(t, total_I, 'r-', label='Infected', linewidth=2)
        ax1.plot(t, total_R, 'g-', label='Recovered', linewidth=2)
        ax1.set_xlabel('Time (days)')
        ax1.set_ylabel('Number of individuals')
        ax1.set_title('Total Population Dynamics')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Age-specific infection rates
        ax2 = axes[0, 1]
        infection_rates = I / self.N * 100
        for i, age_group in enumerate(self.age_groups):
            ax2.plot(t, infection_rates[:, i], label=f'Age {age_group}', linewidth=2)
        ax2.set_xlabel('Time (days)')
        ax2.set_ylabel('Infection Rate (%)')
        ax2.set_title('Age-Specific Infection Rates')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Age-specific attack rates
        ax3 = axes[0, 2]
        final_attack_rates = R[-1, :] / self.N * 100
        bars = ax3.bar(self.age_groups, final_attack_rates, 
                      color=['skyblue', 'lightgreen', 'lightcoral', 'lightyellow', 'lightpink', 'lightgray'])
        ax3.set_xlabel('Age Group')
        ax3.set_ylabel('Final Attack Rate (%)')
        ax3.set_title('Age-Specific Final Attack Rates')
        ax3.tick_params(axis='x', rotation=45)
        
        # Add value labels on bars
        for i, v in enumerate(final_attack_rates):
            ax3.text(i, v + 0.5, f'{v:.1f}%', ha='center', va='bottom')
        
        # Plot 4: Contact matrix heatmap
        ax4 = axes[1, 0]
        sns.heatmap(self.contact_matrix, annot=True, fmt='.2f', cmap='Blues',
                   xticklabels=self.age_groups, yticklabels=self.age_groups, ax=ax4)
        ax4.set_title('Contact Matrix')
        ax4.set_xlabel('Age Group (Contactee)')
        ax4.set_ylabel('Age Group (Contactor)')
        
        # Plot 5: Transmission matrix heatmap
        ax5 = axes[1, 1]
        sns.heatmap(self.beta_matrix, annot=True, fmt='.3f', cmap='Reds',
                   xticklabels=self.age_groups, yticklabels=self.age_groups, ax=ax5)
        ax5.set_title('Transmission Rate Matrix')
        ax5.set_xlabel('Age Group (Infectious)')
        ax5.set_ylabel('Age Group (Susceptible)')
        
        # Plot 6: Age-specific R0
        ax6 = axes[1, 2]
        bars = ax6.bar(self.age_groups, self.R0_age, 
                      color=['red', 'orange', 'yellow', 'green', 'blue', 'purple'])
        ax6.set_xlabel('Age Group')
        ax6.set_ylabel('Basic Reproduction Number (R₀)')
        ax6.set_title('Age-Specific Basic Reproduction Numbers')
        ax6.tick_params(axis='x', rotation=45)
        
        # Add value labels on bars
        for i, v in enumerate(self.R0_age):
            ax6.text(i, v + 0.05, f'{v:.2f}', ha='center', va='bottom')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def compare_with_basic_sir(self, t_span: Tuple[float, float] = (0, 365)) -> None:
        """
        Compare age-structured model with equivalent basic SIR model.
        
        Parameters:
        -----------
        t_span : tuple
            Time span for comparison
        """
        # Create equivalent basic SIR model
        from basic_sir import MalariaSIRModel
        
        # Calculate population-weighted average parameters
        avg_beta = np.mean(np.diag(self.beta_matrix))
        avg_gamma = np.mean(self.gamma)
        avg_mu = np.mean(self.mu)
        
        basic_model = MalariaSIRModel(
            beta=avg_beta,
            gamma=avg_gamma,
            mu=avg_mu,
            N=self.N_total,
            I0=np.sum(self.I0),
            R0=0
        )
        
        # Simulate both models
        t_age, solution_age = self.simulate(t_span)
        t_basic, solution_basic = basic_model.simulate(t_span)
        
        # Extract total populations
        S_age = np.sum(solution_age[:, :self.n_age_groups], axis=1)
        I_age = np.sum(solution_age[:, self.n_age_groups:2*self.n_age_groups], axis=1)
        R_age = np.sum(solution_age[:, 2*self.n_age_groups:], axis=1)
        
        S_basic = solution_basic[:, 0]
        I_basic = solution_basic[:, 1]
        R_basic = solution_basic[:, 2]
        
        # Plot comparison
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # Susceptible comparison
        axes[0].plot(t_age, S_age, 'b-', label='Age-structured', linewidth=2)
        axes[0].plot(t_basic, S_basic, 'b--', label='Basic SIR', linewidth=2)
        axes[0].set_xlabel('Time (days)')
        axes[0].set_ylabel('Susceptible')
        axes[0].set_title('Susceptible Population')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Infected comparison
        axes[1].plot(t_age, I_age, 'r-', label='Age-structured', linewidth=2)
        axes[1].plot(t_basic, I_basic, 'r--', label='Basic SIR', linewidth=2)
        axes[1].set_xlabel('Time (days)')
        axes[1].set_ylabel('Infected')
        axes[1].set_title('Infected Population')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        
        # Recovered comparison
        axes[2].plot(t_age, R_age, 'g-', label='Age-structured', linewidth=2)
        axes[2].plot(t_basic, R_basic, 'g--', label='Basic SIR', linewidth=2)
        axes[2].set_xlabel('Time (days)')
        axes[2].set_ylabel('Recovered')
        axes[2].set_title('Recovered Population')
        axes[2].legend()
        axes[2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()


def main():
    """
    Main function to demonstrate the age-structured malaria SIR model.
    """
    print("Age-Structured Malaria SIR Model Simulation")
    print("=" * 50)
    
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
                                    save_path="../results/plots/malaria_age_structured.png")
    
    # Compare with basic SIR
    print("Comparing with basic SIR model...")
    model.compare_with_basic_sir()


if __name__ == "__main__":
    main()


