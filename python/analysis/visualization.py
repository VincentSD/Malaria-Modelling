"""
Advanced Visualization Tools for Malaria Models
==============================================

This module provides comprehensive visualization tools for malaria transmission models,
including comparison plots, parameter sensitivity analysis, and interactive visualizations.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import List, Dict, Tuple, Optional
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo


class MalariaVisualizer:
    """
    Advanced visualization class for malaria transmission models.
    """
    
    def __init__(self, style: str = 'seaborn-v0_8'):
        """
        Initialize the visualizer.
        
        Parameters:
        -----------
        style : str
            Matplotlib style to use
        """
        plt.style.use(style)
        sns.set_palette("husl")
        
    def plot_model_comparison(self, models: Dict[str, Dict], 
                            t_span: Tuple[float, float] = (0, 365),
                            save_path: str = None) -> None:
        """
        Compare multiple malaria models side by side.
        
        Parameters:
        -----------
        models : dict
            Dictionary with model names as keys and model objects as values
        t_span : tuple
            Time span for comparison
        save_path : str, optional
            Path to save the plot
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
        
        for i, (model_name, model_data) in enumerate(models.items()):
            model = model_data['model']
            t, solution = model.simulate(t_span)
            
            if hasattr(model, 'n_age_groups'):
                # Age-structured model
                S = np.sum(solution[:, :model.n_age_groups], axis=1)
                I = np.sum(solution[:, model.n_age_groups:2*model.n_age_groups], axis=1)
                R = np.sum(solution[:, 2*model.n_age_groups:], axis=1)
            else:
                # Basic model
                S, I, R = solution.T
            
            color = colors[i % len(colors)]
            
            # Plot susceptible
            axes[0, 0].plot(t, S, color=color, label=model_name, linewidth=2)
            
            # Plot infected
            axes[0, 1].plot(t, I, color=color, label=model_name, linewidth=2)
            
            # Plot recovered
            axes[1, 0].plot(t, R, color=color, label=model_name, linewidth=2)
            
            # Plot infection rate
            infection_rate = I / model.N * 100
            axes[1, 1].plot(t, infection_rate, color=color, label=model_name, linewidth=2)
        
        # Format plots
        axes[0, 0].set_title('Susceptible Population')
        axes[0, 0].set_xlabel('Time (days)')
        axes[0, 0].set_ylabel('Number of individuals')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        axes[0, 1].set_title('Infected Population')
        axes[0, 1].set_xlabel('Time (days)')
        axes[0, 1].set_ylabel('Number of individuals')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        axes[1, 0].set_title('Recovered Population')
        axes[1, 0].set_xlabel('Time (days)')
        axes[1, 0].set_ylabel('Number of individuals')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        axes[1, 1].set_title('Infection Rate')
        axes[1, 1].set_xlabel('Time (days)')
        axes[1, 1].set_ylabel('Infection Rate (%)')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def plot_parameter_sensitivity_heatmap(self, sensitivity_results: Dict[str, pd.DataFrame],
                                         save_path: str = None) -> None:
        """
        Create heatmap visualization of parameter sensitivity analysis.
        
        Parameters:
        -----------
        sensitivity_results : dict
            Dictionary with parameter names as keys and sensitivity DataFrames as values
        save_path : str, optional
            Path to save the plot
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        metrics = ['R0_basic', 'peak_infection_rate', 'final_attack_rate', 'time_to_peak']
        metric_titles = ['Basic R₀', 'Peak Infection Rate', 'Final Attack Rate', 'Time to Peak']
        
        for i, (metric, title) in enumerate(zip(metrics, metric_titles)):
            ax = axes[i // 2, i % 2]
            
            # Create pivot table for heatmap
            heatmap_data = []
            param_names = []
            
            for param_name, df in sensitivity_results.items():
                if metric in df.columns:
                    param_names.append(param_name)
                    heatmap_data.append(df[metric].values)
            
            if heatmap_data:
                heatmap_data = np.array(heatmap_data)
                
                # Create heatmap
                im = ax.imshow(heatmap_data, cmap='viridis', aspect='auto')
                
                # Set labels
                ax.set_xticks(range(len(sensitivity_results[param_names[0]]['value'])))
                ax.set_xticklabels([f'{val:.2f}' for val in sensitivity_results[param_names[0]]['value']])
                ax.set_yticks(range(len(param_names)))
                ax.set_yticklabels(param_names)
                
                # Add colorbar
                plt.colorbar(im, ax=ax)
                
                ax.set_title(f'{title} Sensitivity')
                ax.set_xlabel('Parameter Value')
                ax.set_ylabel('Parameter')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def plot_intervention_effectiveness(self, intervention_results: List[Dict],
                                      save_path: str = None) -> None:
        """
        Plot intervention effectiveness comparison.
        
        Parameters:
        -----------
        intervention_results : list
            List of intervention result dictionaries
        save_path : str, optional
            Path to save the plot
        """
        # Extract data for plotting
        intervention_types = []
        peak_reductions = []
        attack_reductions = []
        
        for result in intervention_results:
            intervention_types.append(result['intervention_type'])
            peak_reductions.append(result['reduction_in_peak_infection'] * 100)
            attack_reductions.append(result['reduction_in_attack_rate'] * 100)
        
        # Create bar plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        x = np.arange(len(intervention_types))
        width = 0.35
        
        bars1 = ax1.bar(x - width/2, peak_reductions, width, label='Peak Infection Reduction', 
                       color='skyblue', alpha=0.8)
        bars2 = ax2.bar(x - width/2, attack_reductions, width, label='Attack Rate Reduction', 
                       color='lightcoral', alpha=0.8)
        
        # Add value labels on bars
        for bar in bars1:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{height:.1f}%', ha='center', va='bottom')
        
        for bar in bars2:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{height:.1f}%', ha='center', va='bottom')
        
        # Format plots
        ax1.set_title('Intervention Effectiveness: Peak Infection Reduction')
        ax1.set_xlabel('Intervention Type')
        ax1.set_ylabel('Reduction (%)')
        ax1.set_xticks(x)
        ax1.set_xticklabels(intervention_types, rotation=45)
        ax1.grid(True, alpha=0.3)
        
        ax2.set_title('Intervention Effectiveness: Attack Rate Reduction')
        ax2.set_xlabel('Intervention Type')
        ax2.set_ylabel('Reduction (%)')
        ax2.set_xticks(x)
        ax2.set_xticklabels(intervention_types, rotation=45)
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def create_interactive_plot(self, model, t_span: Tuple[float, float] = (0, 365),
                              save_path: str = None) -> None:
        """
        Create interactive plotly visualization.
        
        Parameters:
        -----------
        model : MalariaModel
            Model to visualize
        t_span : tuple
            Time span for simulation
        save_path : str, optional
            Path to save the HTML file
        """
        t, solution = model.simulate(t_span)
        
        if hasattr(model, 'n_age_groups'):
            # Age-structured model
            S = np.sum(solution[:, :model.n_age_groups], axis=1)
            I = np.sum(solution[:, model.n_age_groups:2*model.n_age_groups], axis=1)
            R = np.sum(solution[:, 2*model.n_age_groups:], axis=1)
            
            # Age-specific data
            age_data = []
            for i in range(model.n_age_groups):
                age_data.append({
                    'time': t,
                    'infected': solution[:, model.n_age_groups + i],
                    'age_group': model.age_groups[i]
                })
        else:
            # Basic model
            S, I, R = solution.T
            age_data = None
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Population Dynamics', 'Infection Rate', 
                          'Phase Plane', 'Cumulative Cases'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        # Population dynamics
        fig.add_trace(
            go.Scatter(x=t, y=S, mode='lines', name='Susceptible', 
                      line=dict(color='blue', width=2)),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=t, y=I, mode='lines', name='Infected', 
                      line=dict(color='red', width=2)),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=t, y=R, mode='lines', name='Recovered', 
                      line=dict(color='green', width=2)),
            row=1, col=1
        )
        
        # Infection rate
        infection_rate = I / model.N * 100
        fig.add_trace(
            go.Scatter(x=t, y=infection_rate, mode='lines', name='Infection Rate (%)', 
                      line=dict(color='purple', width=2)),
            row=1, col=2
        )
        
        # Phase plane
        fig.add_trace(
            go.Scatter(x=S, y=I, mode='lines', name='Phase Plane', 
                      line=dict(color='orange', width=2)),
            row=2, col=1
        )
        
        # Cumulative cases
        cumulative_cases = R + I
        fig.add_trace(
            go.Scatter(x=t, y=cumulative_cases, mode='lines', name='Cumulative Cases', 
                      line=dict(color='brown', width=2)),
            row=2, col=2
        )
        
        # Update layout
        fig.update_layout(
            title_text="Interactive Malaria Model Visualization",
            showlegend=True,
            height=800
        )
        
        # Update axes labels
        fig.update_xaxes(title_text="Time (days)", row=1, col=1)
        fig.update_yaxes(title_text="Number of individuals", row=1, col=1)
        
        fig.update_xaxes(title_text="Time (days)", row=1, col=2)
        fig.update_yaxes(title_text="Infection Rate (%)", row=1, col=2)
        
        fig.update_xaxes(title_text="Susceptible", row=2, col=1)
        fig.update_yaxes(title_text="Infected", row=2, col=1)
        
        fig.update_xaxes(title_text="Time (days)", row=2, col=2)
        fig.update_yaxes(title_text="Cumulative Cases", row=2, col=2)
        
        # Show plot
        fig.show()
        
        # Save if path provided
        if save_path:
            fig.write_html(save_path)
    
    def plot_epidemic_curves_comparison(self, models_data: List[Dict],
                                      save_path: str = None) -> None:
        """
        Compare epidemic curves from different models.
        
        Parameters:
        -----------
        models_data : list
            List of dictionaries containing model data
        save_path : str, optional
            Path to save the plot
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        for model_data in models_data:
            model_name = model_data['name']
            t = model_data['time']
            solution = model_data['solution']
            model = model_data['model']
            
            if hasattr(model, 'n_age_groups'):
                I = np.sum(solution[:, model.n_age_groups:2*model.n_age_groups], axis=1)
            else:
                I = solution[:, 1]  # Infected compartment
            
            # Normalize by population
            I_rate = I / model.N * 100
            
            # Plot on all subplots
            for ax in axes.flat:
                ax.plot(t, I_rate, label=model_name, linewidth=2)
        
        # Format plots
        titles = ['Infection Rate Over Time', 'Log Scale Infection Rate', 
                 'Cumulative Infection Rate', 'Peak Comparison']
        
        for i, (ax, title) in enumerate(zip(axes.flat, titles)):
            ax.set_title(title)
            ax.set_xlabel('Time (days)')
            ax.set_ylabel('Infection Rate (%)')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            if i == 1:  # Log scale
                ax.set_yscale('log')
            elif i == 2:  # Cumulative
                for model_data in models_data:
                    model_name = model_data['name']
                    t = model_data['time']
                    solution = model_data['solution']
                    model = model_data['model']
                    
                    if hasattr(model, 'n_age_groups'):
                        R = np.sum(solution[:, 2*model.n_age_groups:], axis=1)
                    else:
                        R = solution[:, 2]  # Recovered compartment
                    
                    cumulative_rate = R / model.N * 100
                    ax.plot(t, cumulative_rate, label=f'{model_name} (Cumulative)', 
                           linewidth=2, linestyle='--')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()


def create_comprehensive_analysis(models: Dict[str, any], 
                                t_span: Tuple[float, float] = (0, 1000),
                                save_dir: str = "../results/plots/") -> None:
    """
    Create comprehensive analysis and visualization of multiple malaria models.
    
    Parameters:
    -----------
    models : dict
        Dictionary of model objects
    t_span : tuple
        Time span for analysis
    save_dir : str
        Directory to save plots
    """
    visualizer = MalariaVisualizer()
    
    # Prepare model data
    models_data = {}
    for name, model in models.items():
        t, solution = model.simulate(t_span)
        models_data[name] = {
            'model': model,
            'time': t,
            'solution': solution
        }
    
    # Model comparison
    print("Creating model comparison plots...")
    visualizer.plot_model_comparison(models_data, t_span, 
                                   save_path=f"{save_dir}model_comparison.png")
    
    # Interactive plots for each model
    for name, model in models.items():
        print(f"Creating interactive plot for {name}...")
        visualizer.create_interactive_plot(model, t_span, 
                                         save_path=f"{save_dir}{name}_interactive.html")
    
    # Epidemic curves comparison
    print("Creating epidemic curves comparison...")
    models_list = [{'name': name, **data} for name, data in models_data.items()]
    visualizer.plot_epidemic_curves_comparison(models_list, 
                                             save_path=f"{save_dir}epidemic_curves.png")
    
    print(f"Analysis complete! Plots saved to {save_dir}")


if __name__ == "__main__":
    # Example usage
    from models.basic_sir import MalariaSIRModel
    from models.seir_model import MalariaSEIRModel
    from models.vector_dynamics import MalariaVectorModel
    
    # Create sample models
    models = {
        'Basic SIR': MalariaSIRModel(beta=0.4, gamma=0.05, N=10000, I0=50),
        'SEIR': MalariaSEIRModel(beta=0.4, sigma=0.1, gamma=0.05, N=10000, I0=30, E0=20),
        'Vector Model': MalariaVectorModel(N_h=10000, N_v=50000, I_h0=30, E_h0=20)
    }
    
    # Run comprehensive analysis
    create_comprehensive_analysis(models)


