"""
Parameter Analysis Tools for Malaria Models
==========================================

This module provides comprehensive parameter analysis tools for malaria transmission models,
including sensitivity analysis, parameter estimation, and uncertainty quantification.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize, differential_evolution
from scipy.stats import uniform, norm
from typing import Dict, List, Tuple, Callable, Optional
import warnings
warnings.filterwarnings('ignore')


class ParameterAnalyzer:
    """
    Comprehensive parameter analysis for malaria transmission models.
    """
    
    def __init__(self, model_class, model_params: Dict):
        """
        Initialize the parameter analyzer.
        
        Parameters:
        -----------
        model_class : class
            Model class to analyze
        model_params : dict
            Base parameters for the model
        """
        self.model_class = model_class
        self.model_params = model_params
        self.results = {}
        
    def sensitivity_analysis(self, 
                           param_ranges: Dict[str, Tuple[float, float]],
                           num_points: int = 10,
                           t_span: Tuple[float, float] = (0, 365),
                           metrics: List[str] = None) -> pd.DataFrame:
        """
        Perform comprehensive sensitivity analysis.
        
        Parameters:
        -----------
        param_ranges : dict
            Dictionary with parameter names as keys and (min, max) tuples as values
        num_points : int
            Number of parameter values to test
        t_span : tuple
            Time span for simulation
        metrics : list
            List of metrics to calculate
            
        Returns:
        --------
        pd.DataFrame
            Results of sensitivity analysis
        """
        if metrics is None:
            metrics = ['peak_infection_rate', 'final_attack_rate', 'time_to_peak']
        
        results = []
        
        for param_name, (min_val, max_val) in param_ranges.items():
            print(f"Analyzing parameter: {param_name}")
            
            # Create parameter values
            param_values = np.linspace(min_val, max_val, num_points)
            
            for value in param_values:
                # Create model with modified parameter
                modified_params = self.model_params.copy()
                modified_params[param_name] = value
                
                try:
                    model = self.model_class(**modified_params)
                    t, solution = model.simulate(t_span)
                    
                    # Calculate metrics
                    if hasattr(model, 'get_epidemic_metrics'):
                        model_metrics = model.get_epidemic_metrics(solution)
                    else:
                        model_metrics = self._calculate_basic_metrics(model, solution)
                    
                    # Store results
                    result = {
                        'parameter': param_name,
                        'value': value,
                        'normalized_value': (value - min_val) / (max_val - min_val)
                    }
                    
                    for metric in metrics:
                        if metric in model_metrics:
                            result[metric] = model_metrics[metric]
                        else:
                            result[metric] = np.nan
                    
                    results.append(result)
                    
                except Exception as e:
                    print(f"Error with {param_name}={value}: {e}")
                    continue
        
        return pd.DataFrame(results)
    
    def _calculate_basic_metrics(self, model, solution: np.ndarray) -> Dict[str, float]:
        """
        Calculate basic epidemic metrics.
        
        Parameters:
        -----------
        model : MalariaModel
            Model object
        solution : np.ndarray
            Solution array
            
        Returns:
        --------
        dict
            Basic metrics
        """
        if hasattr(model, 'n_age_groups'):
            # Age-structured model
            I = np.sum(solution[:, model.n_age_groups:2*model.n_age_groups], axis=1)
            R = np.sum(solution[:, 2*model.n_age_groups:], axis=1)
        else:
            # Basic model
            I = solution[:, 1]  # Infected
            R = solution[:, 2]  # Recovered
        
        return {
            'peak_infection_rate': np.max(I) / model.N,
            'final_attack_rate': R[-1] / model.N,
            'time_to_peak': np.argmax(I)
        }
    
    def parameter_estimation(self, 
                           observed_data: Dict[str, np.ndarray],
                           param_bounds: Dict[str, Tuple[float, float]],
                           t_span: Tuple[float, float] = (0, 365),
                           method: str = 'differential_evolution') -> Dict:
        """
        Estimate model parameters from observed data.
        
        Parameters:
        -----------
        observed_data : dict
            Dictionary with 'time' and observed values
        param_bounds : dict
            Parameter bounds for estimation
        t_span : tuple
            Time span for simulation
        method : str
            Optimization method ('differential_evolution' or 'minimize')
            
        Returns:
        --------
        dict
            Estimation results
        """
        def objective_function(params):
            """Objective function for parameter estimation."""
            # Create model with current parameters
            modified_params = self.model_params.copy()
            for i, param_name in enumerate(param_bounds.keys()):
                modified_params[param_name] = params[i]
            
            try:
                model = self.model_class(**modified_params)
                t, solution = model.simulate(t_span)
                
                # Calculate model predictions
                if hasattr(model, 'n_age_groups'):
                    I_pred = np.sum(solution[:, model.n_age_groups:2*model.n_age_groups], axis=1)
                else:
                    I_pred = solution[:, 1]
                
                # Interpolate to match observed data time points
                I_pred_interp = np.interp(observed_data['time'], t, I_pred)
                
                # Calculate sum of squared errors
                sse = np.sum((observed_data['infected'] - I_pred_interp)**2)
                
                return sse
                
            except Exception as e:
                return 1e10  # Large penalty for invalid parameters
        
        # Prepare bounds for optimization
        bounds = list(param_bounds.values())
        param_names = list(param_bounds.keys())
        
        # Run optimization
        if method == 'differential_evolution':
            result = differential_evolution(objective_function, bounds, seed=42)
        else:
            # Initial guess (middle of bounds)
            x0 = [(bounds[i][0] + bounds[i][1]) / 2 for i in range(len(bounds))]
            result = minimize(objective_function, x0, bounds=bounds)
        
        # Create results dictionary
        estimated_params = {}
        for i, param_name in enumerate(param_names):
            estimated_params[param_name] = result.x[i]
        
        return {
            'estimated_parameters': estimated_params,
            'optimization_success': result.success,
            'objective_value': result.fun,
            'optimization_result': result
        }
    
    def uncertainty_quantification(self, 
                                 param_distributions: Dict[str, Tuple[str, float, float]],
                                 num_samples: int = 1000,
                                 t_span: Tuple[float, float] = (0, 365)) -> Dict:
        """
        Perform uncertainty quantification using Monte Carlo sampling.
        
        Parameters:
        -----------
        param_distributions : dict
            Dictionary with parameter names as keys and (distribution, param1, param2) tuples
        num_samples : int
            Number of Monte Carlo samples
        t_span : tuple
            Time span for simulation
            
        Returns:
        --------
        dict
            Uncertainty quantification results
        """
        results = []
        
        for i in range(num_samples):
            if i % 100 == 0:
                print(f"Monte Carlo sample {i}/{num_samples}")
            
            # Sample parameters
            sampled_params = self.model_params.copy()
            
            for param_name, (dist_type, param1, param2) in param_distributions.items():
                if dist_type == 'uniform':
                    value = uniform.rvs(param1, param2 - param1)
                elif dist_type == 'normal':
                    value = norm.rvs(param1, param2)
                else:
                    raise ValueError(f"Unknown distribution type: {dist_type}")
                
                sampled_params[param_name] = value
            
            try:
                # Create and simulate model
                model = self.model_class(**sampled_params)
                t, solution = model.simulate(t_span)
                
                # Calculate metrics
                if hasattr(model, 'get_epidemic_metrics'):
                    metrics = model.get_epidemic_metrics(solution)
                else:
                    metrics = self._calculate_basic_metrics(model, solution)
                
                # Store results
                result = {'sample': i}
                result.update(sampled_params)
                result.update(metrics)
                results.append(result)
                
            except Exception as e:
                print(f"Error in sample {i}: {e}")
                continue
        
        return pd.DataFrame(results)
    
    def plot_sensitivity_results(self, sensitivity_df: pd.DataFrame,
                               save_path: str = None) -> None:
        """
        Plot sensitivity analysis results.
        
        Parameters:
        -----------
        sensitivity_df : pd.DataFrame
            Results from sensitivity analysis
        save_path : str, optional
            Path to save the plot
        """
        # Get unique parameters and metrics
        parameters = sensitivity_df['parameter'].unique()
        metrics = [col for col in sensitivity_df.columns 
                  if col not in ['parameter', 'value', 'normalized_value']]
        
        # Create subplots
        n_params = len(parameters)
        n_metrics = len(metrics)
        
        fig, axes = plt.subplots(n_params, n_metrics, figsize=(5*n_metrics, 4*n_params))
        
        if n_params == 1:
            axes = axes.reshape(1, -1)
        if n_metrics == 1:
            axes = axes.reshape(-1, 1)
        
        for i, param in enumerate(parameters):
            param_data = sensitivity_df[sensitivity_df['parameter'] == param]
            
            for j, metric in enumerate(metrics):
                ax = axes[i, j]
                
                # Plot sensitivity curve
                ax.plot(param_data['normalized_value'], param_data[metric], 
                       'o-', linewidth=2, markersize=6)
                
                ax.set_title(f'{param} vs {metric}')
                ax.set_xlabel('Normalized Parameter Value')
                ax.set_ylabel(metric)
                ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def plot_uncertainty_results(self, uncertainty_df: pd.DataFrame,
                               save_path: str = None) -> None:
        """
        Plot uncertainty quantification results.
        
        Parameters:
        -----------
        uncertainty_df : pd.DataFrame
            Results from uncertainty quantification
        save_path : str, optional
            Path to save the plot
        """
        # Get metric columns
        metric_cols = [col for col in uncertainty_df.columns 
                      if col not in ['sample'] and col not in self.model_params.keys()]
        
        # Create subplots
        n_metrics = len(metric_cols)
        fig, axes = plt.subplots(2, (n_metrics + 1) // 2, figsize=(15, 10))
        axes = axes.flatten()
        
        for i, metric in enumerate(metric_cols):
            ax = axes[i]
            
            # Create histogram
            ax.hist(uncertainty_df[metric], bins=50, alpha=0.7, density=True)
            ax.axvline(uncertainty_df[metric].mean(), color='red', linestyle='--', 
                      label=f'Mean: {uncertainty_df[metric].mean():.3f}')
            ax.axvline(uncertainty_df[metric].median(), color='green', linestyle='--', 
                      label=f'Median: {uncertainty_df[metric].median():.3f}')
            
            ax.set_title(f'Distribution of {metric}')
            ax.set_xlabel(metric)
            ax.set_ylabel('Density')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for i in range(n_metrics, len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def correlation_analysis(self, uncertainty_df: pd.DataFrame) -> pd.DataFrame:
        """
        Perform correlation analysis between parameters and outcomes.
        
        Parameters:
        -----------
        uncertainty_df : pd.DataFrame
            Results from uncertainty quantification
            
        Returns:
        --------
        pd.DataFrame
            Correlation matrix
        """
        # Get parameter and metric columns
        param_cols = [col for col in uncertainty_df.columns 
                     if col in self.model_params.keys()]
        metric_cols = [col for col in uncertainty_df.columns 
                      if col not in ['sample'] and col not in param_cols]
        
        # Calculate correlations
        correlation_data = []
        
        for param in param_cols:
            for metric in metric_cols:
                corr = uncertainty_df[param].corr(uncertainty_df[metric])
                correlation_data.append({
                    'parameter': param,
                    'metric': metric,
                    'correlation': corr,
                    'abs_correlation': abs(corr)
                })
        
        correlation_df = pd.DataFrame(correlation_data)
        
        # Create correlation heatmap
        pivot_table = correlation_df.pivot(index='parameter', columns='metric', values='correlation')
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(pivot_table, annot=True, cmap='RdBu_r', center=0, 
                   fmt='.3f', cbar_kws={'label': 'Correlation Coefficient'})
        plt.title('Parameter-Metric Correlation Matrix')
        plt.tight_layout()
        plt.show()
        
        return correlation_df


def comprehensive_parameter_analysis(model_class, model_params: Dict,
                                   param_ranges: Dict[str, Tuple[float, float]],
                                   save_dir: str = "../results/analysis/") -> None:
    """
    Perform comprehensive parameter analysis.
    
    Parameters:
    -----------
    model_class : class
        Model class to analyze
    model_params : dict
        Base parameters for the model
    param_ranges : dict
        Parameter ranges for sensitivity analysis
    save_dir : str
        Directory to save results
    """
    analyzer = ParameterAnalyzer(model_class, model_params)
    
    print("Starting comprehensive parameter analysis...")
    
    # Sensitivity analysis
    print("1. Performing sensitivity analysis...")
    sensitivity_results = analyzer.sensitivity_analysis(param_ranges)
    analyzer.plot_sensitivity_results(sensitivity_results, 
                                    save_path=f"{save_dir}sensitivity_analysis.png")
    
    # Save sensitivity results
    sensitivity_results.to_csv(f"{save_dir}sensitivity_results.csv", index=False)
    
    # Uncertainty quantification
    print("2. Performing uncertainty quantification...")
    param_distributions = {}
    for param, (min_val, max_val) in param_ranges.items():
        # Use uniform distribution centered on base value
        base_val = model_params[param]
        param_distributions[param] = ('uniform', min_val, max_val)
    
    uncertainty_results = analyzer.uncertainty_quantification(param_distributions, num_samples=500)
    analyzer.plot_uncertainty_results(uncertainty_results, 
                                    save_path=f"{save_dir}uncertainty_analysis.png")
    
    # Save uncertainty results
    uncertainty_results.to_csv(f"{save_dir}uncertainty_results.csv", index=False)
    
    # Correlation analysis
    print("3. Performing correlation analysis...")
    correlation_results = analyzer.correlation_analysis(uncertainty_results)
    correlation_results.to_csv(f"{save_dir}correlation_analysis.csv", index=False)
    
    print(f"Analysis complete! Results saved to {save_dir}")


if __name__ == "__main__":
    # Example usage
    from models.basic_sir import MalariaSIRModel
    
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
    comprehensive_parameter_analysis(MalariaSIRModel, model_params, param_ranges)


