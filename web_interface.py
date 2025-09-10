#!/usr/bin/env python3
"""
Malaria Model Simulation - Web Interface
=======================================

Interactive web-based dashboard for malaria transmission models.
Built with Streamlit for easy parameter adjustment and visualization.
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import sys
import os

# Add the python directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), 'python'))

from models.basic_sir import MalariaSIRModel
from models.seir_model import MalariaSEIRModel
from models.age_structured_sir import AgeStructuredMalariaSIR
from models.vector_dynamics import MalariaVectorModel

# Page configuration
st.set_page_config(
    page_title="🦟 Malaria Model Simulation Dashboard",
    page_icon="🦟",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .model-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 1rem 0;
    }
    .metric-box {
        background-color: #e8f4fd;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
    }
    .parameter-explanation {
        background-color: #fff3cd;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #ffc107;
    }
</style>
""", unsafe_allow_html=True)

def main():
    # Header
    st.markdown('<h1 class="main-header">🦟 Malaria Model Simulation Dashboard</h1>', unsafe_allow_html=True)
    
    # Sidebar for model selection and parameters
    st.sidebar.title("🎛️ Model Controls")
    
    # Model selection
    model_type = st.sidebar.selectbox(
        "Select Model Type",
        ["Basic SIR", "SEIR", "Age-Structured SIR", "Vector Dynamics"],
        help="Choose the complexity level of the malaria transmission model"
    )
    
    st.sidebar.markdown("---")
    
    # Model-specific parameters
    if model_type == "Basic SIR":
        render_basic_sir_interface()
    elif model_type == "SEIR":
        render_seir_interface()
    elif model_type == "Age-Structured SIR":
        render_age_structured_interface()
    elif model_type == "Vector Dynamics":
        render_vector_dynamics_interface()

def render_basic_sir_interface():
    st.sidebar.subheader("📊 Basic SIR Parameters")
    
    # Parameter inputs
    beta = st.sidebar.slider(
        "Transmission Rate (β)", 
        min_value=0.1, max_value=1.0, value=0.4, step=0.05,
        help="Rate of infection from infected to susceptible individuals per day"
    )
    
    gamma = st.sidebar.slider(
        "Recovery Rate (γ)", 
        min_value=0.01, max_value=0.3, value=0.05, step=0.01,
        help="Rate of recovery from infection per day (1/infectious period)"
    )
    
    mu = st.sidebar.slider(
        "Birth/Death Rate (μ)", 
        min_value=0.00001, max_value=0.001, value=0.0001, step=0.00001,
        help="Population turnover rate per day"
    )
    
    N = st.sidebar.number_input(
        "Population Size (N)", 
        min_value=1000, max_value=100000, value=10000, step=1000,
        help="Total population size"
    )
    
    I0 = st.sidebar.number_input(
        "Initial Infected (I₀)", 
        min_value=1, max_value=1000, value=50, step=1,
        help="Number of initially infected individuals"
    )
    
    simulation_days = st.sidebar.slider(
        "Simulation Duration (days)", 
        min_value=100, max_value=2000, value=1000, step=100
    )
    
    # Create model and simulate
    if st.sidebar.button("🚀 Run Simulation", type="primary"):
        with st.spinner("Running simulation..."):
            model = MalariaSIRModel(beta=beta, gamma=gamma, mu=mu, N=N, I0=I0, R0=0)
            t, solution = model.simulate(t_span=(0, simulation_days), num_points=1000)
            metrics = model.get_epidemic_metrics(solution)
            
            # Store results in session state
            st.session_state['model'] = model
            st.session_state['t'] = t
            st.session_state['solution'] = solution
            st.session_state['metrics'] = metrics
            st.session_state['model_type'] = 'Basic SIR'
    
    # Display results if available
    if 'model' in st.session_state and st.session_state['model_type'] == 'Basic SIR':
        display_basic_sir_results()

def render_seir_interface():
    st.sidebar.subheader("📊 SEIR Parameters")
    
    beta = st.sidebar.slider("Transmission Rate (β)", 0.1, 1.0, 0.4, 0.05)
    sigma = st.sidebar.slider("Incubation Rate (σ)", 0.05, 0.5, 0.1, 0.01)
    gamma = st.sidebar.slider("Recovery Rate (γ)", 0.01, 0.3, 0.05, 0.01)
    mu = st.sidebar.slider("Birth/Death Rate (μ)", 0.00001, 0.001, 0.0001, 0.00001)
    N = st.sidebar.number_input("Population Size (N)", 1000, 100000, 10000, 1000)
    E0 = st.sidebar.number_input("Initial Exposed (E₀)", 1, 1000, 20, 1)
    I0 = st.sidebar.number_input("Initial Infected (I₀)", 1, 1000, 30, 1)
    seasonal_amp = st.sidebar.slider("Seasonal Amplitude", 0.0, 0.8, 0.3, 0.1)
    simulation_days = st.sidebar.slider("Simulation Duration (days)", 100, 2000, 1000, 100)
    
    if st.sidebar.button("🚀 Run SEIR Simulation", type="primary"):
        with st.spinner("Running SEIR simulation..."):
            model = MalariaSEIRModel(
                beta=beta, sigma=sigma, gamma=gamma, mu=mu, N=N,
                E0=E0, I0=I0, R0=0, seasonal_amplitude=seasonal_amp
            )
            t, solution = model.simulate(t_span=(0, simulation_days), num_points=1000)
            metrics = model.get_epidemic_metrics(solution)
            
            st.session_state['model'] = model
            st.session_state['t'] = t
            st.session_state['solution'] = solution
            st.session_state['metrics'] = metrics
            st.session_state['model_type'] = 'SEIR'
    
    if 'model' in st.session_state and st.session_state['model_type'] == 'SEIR':
        display_seir_results()

def render_age_structured_interface():
    st.sidebar.subheader("📊 Age-Structured Parameters")
    
    N_total = st.sidebar.number_input("Total Population", 10000, 500000, 100000, 10000)
    
    st.sidebar.markdown("**Initial Infected by Age Group:**")
    age_groups = ['0-4', '5-14', '15-29', '30-44', '45-59', '60+']
    I0_age = []
    for i, age_group in enumerate(age_groups):
        I0_age.append(st.sidebar.number_input(f"Age {age_group}", 0, 1000, [10, 20, 30, 20, 10, 10][i], 1))
    
    simulation_days = st.sidebar.slider("Simulation Duration (days)", 100, 2000, 1000, 100)
    
    if st.sidebar.button("🚀 Run Age-Structured Simulation", type="primary"):
        with st.spinner("Running age-structured simulation..."):
            model = AgeStructuredMalariaSIR(N_total=N_total, I0=np.array(I0_age))
            t, solution = model.simulate(t_span=(0, simulation_days), num_points=1000)
            metrics = model.get_age_specific_metrics(solution)
            
            st.session_state['model'] = model
            st.session_state['t'] = t
            st.session_state['solution'] = solution
            st.session_state['metrics'] = metrics
            st.session_state['model_type'] = 'Age-Structured SIR'
    
    if 'model' in st.session_state and st.session_state['model_type'] == 'Age-Structured SIR':
        display_age_structured_results()

def render_vector_dynamics_interface():
    st.sidebar.subheader("📊 Vector Dynamics Parameters")
    
    # Human parameters
    st.sidebar.markdown("**Human Population:**")
    N_h = st.sidebar.number_input("Human Population (N_h)", 5000, 100000, 10000, 1000)
    
    # Vector parameters
    st.sidebar.markdown("**Vector Population:**")
    N_v = st.sidebar.number_input("Vector Population (N_v)", 10000, 500000, 50000, 5000)
    
    # Biting parameters
    st.sidebar.markdown("**Transmission Parameters:**")
    a = st.sidebar.slider("Biting Rate (a)", 0.1, 0.8, 0.3, 0.05)
    b = st.sidebar.slider("Human→Vector Transmission (b)", 0.1, 0.9, 0.5, 0.05)
    c = st.sidebar.slider("Vector→Human Transmission (c)", 0.1, 0.9, 0.5, 0.05)
    
    # Initial conditions
    st.sidebar.markdown("**Initial Conditions:**")
    E_h0 = st.sidebar.number_input("Initial Exposed Humans", 1, 1000, 20, 1)
    I_h0 = st.sidebar.number_input("Initial Infected Humans", 1, 1000, 30, 1)
    E_v0 = st.sidebar.number_input("Initial Exposed Vectors", 100, 10000, 500, 100)
    I_v0 = st.sidebar.number_input("Initial Infected Vectors", 100, 10000, 1000, 100)
    
    simulation_days = st.sidebar.slider("Simulation Duration (days)", 100, 2000, 1000, 100)
    
    if st.sidebar.button("🚀 Run Vector Dynamics Simulation", type="primary"):
        with st.spinner("Running vector dynamics simulation..."):
            model = MalariaVectorModel(
                N_h=N_h, N_v=N_v, a=a, b=b, c=c,
                E_h0=E_h0, I_h0=I_h0, E_v0=E_v0, I_v0=I_v0
            )
            t, solution = model.simulate(t_span=(0, simulation_days), num_points=1000)
            metrics = model.get_epidemic_metrics(solution)
            
            st.session_state['model'] = model
            st.session_state['t'] = t
            st.session_state['solution'] = solution
            st.session_state['metrics'] = metrics
            st.session_state['model_type'] = 'Vector Dynamics'
    
    if 'model' in st.session_state and st.session_state['model_type'] == 'Vector Dynamics':
        display_vector_dynamics_results()

def display_basic_sir_results():
    model = st.session_state['model']
    t = st.session_state['t']
    solution = st.session_state['solution']
    metrics = st.session_state['metrics']
    
    # Main results section
    st.header("📊 Basic SIR Model Results")
    
    # Key metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Basic R₀", f"{metrics['R0_basic']:.2f}")
    with col2:
        st.metric("Peak Infection Rate", f"{metrics['peak_infection_rate']*100:.1f}%")
    with col3:
        st.metric("Final Attack Rate", f"{metrics['final_attack_rate']*100:.1f}%")
    with col4:
        st.metric("Time to Peak", f"{metrics['time_to_peak']} days")
    
    # Explanation
    st.markdown("""
    <div class="parameter-explanation">
    <h4>📖 What These Results Mean:</h4>
    <ul>
    <li><strong>Basic R₀ = {:.2f}</strong>: Each infected person infects {:.1f} others on average. 
    Values > 1 indicate epidemic potential.</li>
    <li><strong>Peak Infection Rate = {:.1f}%</strong>: Maximum percentage of population infected simultaneously.</li>
    <li><strong>Final Attack Rate = {:.1f}%</strong>: Total percentage of population that becomes infected.</li>
    <li><strong>Time to Peak = {} days</strong>: How quickly the epidemic reaches its maximum.</li>
    </ul>
    </div>
    """.format(metrics['R0_basic'], metrics['R0_basic'], 
               metrics['peak_infection_rate']*100, metrics['final_attack_rate']*100, 
               metrics['time_to_peak']), unsafe_allow_html=True)
    
    # Interactive plots
    S, I, R = solution.T
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Population Dynamics', 'Infection Rate Over Time', 
                       'Phase Plane (S vs I)', 'Cumulative Cases'),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # Population dynamics
    fig.add_trace(go.Scatter(x=t, y=S, mode='lines', name='Susceptible', 
                            line=dict(color='blue', width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=I, mode='lines', name='Infected', 
                            line=dict(color='red', width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=R, mode='lines', name='Recovered', 
                            line=dict(color='green', width=2)), row=1, col=1)
    
    # Infection rate
    infection_rate = I / model.N * 100
    fig.add_trace(go.Scatter(x=t, y=infection_rate, mode='lines', name='Infection Rate (%)', 
                            line=dict(color='purple', width=2)), row=1, col=2)
    
    # Phase plane
    fig.add_trace(go.Scatter(x=S, y=I, mode='lines', name='Phase Plane', 
                            line=dict(color='orange', width=2)), row=2, col=1)
    
    # Cumulative cases
    cumulative_cases = R + I
    fig.add_trace(go.Scatter(x=t, y=cumulative_cases, mode='lines', name='Cumulative Cases', 
                            line=dict(color='brown', width=2)), row=2, col=2)
    
    # Update layout
    fig.update_layout(height=800, showlegend=True, title_text="Malaria SIR Model Dynamics")
    
    # Update axes labels
    fig.update_xaxes(title_text="Time (days)", row=1, col=1)
    fig.update_yaxes(title_text="Number of individuals", row=1, col=1)
    fig.update_xaxes(title_text="Time (days)", row=1, col=2)
    fig.update_yaxes(title_text="Infection Rate (%)", row=1, col=2)
    fig.update_xaxes(title_text="Susceptible", row=2, col=1)
    fig.update_yaxes(title_text="Infected", row=2, col=1)
    fig.update_xaxes(title_text="Time (days)", row=2, col=2)
    fig.update_yaxes(title_text="Cumulative Cases", row=2, col=2)
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Parameter sensitivity analysis
    st.subheader("🔬 Parameter Sensitivity Analysis")
    
    if st.button("Run Sensitivity Analysis"):
        with st.spinner("Analyzing parameter sensitivity..."):
            # Test different transmission rates
            beta_values = np.linspace(0.1, 0.8, 10)
            results = []
            
            for beta_val in beta_values:
                temp_model = MalariaSIRModel(beta=beta_val, gamma=model.gamma, 
                                           mu=model.mu, N=model.N, I0=model.I0, R0=0)
                t_temp, solution_temp = temp_model.simulate(t_span=(0, 1000), num_points=1000)
                metrics_temp = temp_model.get_epidemic_metrics(solution_temp)
                results.append({
                    'Beta': beta_val,
                    'R0': metrics_temp['R0_basic'],
                    'Peak Infection Rate': metrics_temp['peak_infection_rate'] * 100,
                    'Final Attack Rate': metrics_temp['final_attack_rate'] * 100
                })
            
            df = pd.DataFrame(results)
            
            # Create sensitivity plot
            fig_sens = make_subplots(
                rows=1, cols=3,
                subplot_titles=('R₀ vs Transmission Rate', 'Peak Infection vs Transmission Rate', 
                               'Final Attack Rate vs Transmission Rate')
            )
            
            fig_sens.add_trace(go.Scatter(x=df['Beta'], y=df['R0'], mode='lines+markers', 
                                        name='R₀'), row=1, col=1)
            fig_sens.add_trace(go.Scatter(x=df['Beta'], y=df['Peak Infection Rate'], 
                                        mode='lines+markers', name='Peak Infection Rate'), row=1, col=2)
            fig_sens.add_trace(go.Scatter(x=df['Beta'], y=df['Final Attack Rate'], 
                                        mode='lines+markers', name='Final Attack Rate'), row=1, col=3)
            
            fig_sens.update_layout(height=400, showlegend=False)
            fig_sens.update_xaxes(title_text="Transmission Rate (β)", row=1, col=1)
            fig_sens.update_yaxes(title_text="R₀", row=1, col=1)
            fig_sens.update_xaxes(title_text="Transmission Rate (β)", row=1, col=2)
            fig_sens.update_yaxes(title_text="Peak Infection Rate (%)", row=1, col=2)
            fig_sens.update_xaxes(title_text="Transmission Rate (β)", row=1, col=3)
            fig_sens.update_yaxes(title_text="Final Attack Rate (%)", row=1, col=3)
            
            st.plotly_chart(fig_sens, use_container_width=True)
            
            # Display results table
            st.subheader("📋 Sensitivity Analysis Results")
            st.dataframe(df.round(3))

def display_seir_results():
    model = st.session_state['model']
    t = st.session_state['t']
    solution = st.session_state['solution']
    metrics = st.session_state['metrics']
    
    st.header("📊 SEIR Model Results")
    
    # Key metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Basic R₀", f"{metrics['R0_basic']:.2f}")
    with col2:
        st.metric("Effective R₀", f"{metrics['R0_effective']:.2f}")
    with col3:
        st.metric("Peak Infection Rate", f"{metrics['peak_infection_rate']*100:.1f}%")
    with col4:
        st.metric("Latent Period", f"{metrics['latent_period']:.1f} days")
    
    # SEIR-specific explanation
    st.markdown("""
    <div class="parameter-explanation">
    <h4>📖 SEIR Model Explanation:</h4>
    <p>The SEIR model includes an <strong>Exposed (E)</strong> compartment representing individuals 
    who are infected but not yet infectious. This creates a more realistic latent period 
    before individuals become infectious.</p>
    <ul>
    <li><strong>Latent Period</strong>: Time from infection to becoming infectious</li>
    <li><strong>Effective R₀</strong>: Accounts for the latent period in transmission</li>
    <li><strong>Seasonal Effects</strong>: Transmission varies throughout the year</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # SEIR plots
    S, E, I, R = solution.T
    
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('SEIR Population Dynamics', 'Infection and Exposure Rates', 
                       'Seasonal Transmission Rate', 'Phase Plane (S vs I)')
    )
    
    # Population dynamics
    fig.add_trace(go.Scatter(x=t, y=S, mode='lines', name='Susceptible', 
                            line=dict(color='blue', width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=E, mode='lines', name='Exposed', 
                            line=dict(color='orange', width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=I, mode='lines', name='Infected', 
                            line=dict(color='red', width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=R, mode='lines', name='Recovered', 
                            line=dict(color='green', width=2)), row=1, col=1)
    
    # Infection and exposure rates
    infection_rate = I / model.N * 100
    exposure_rate = E / model.N * 100
    fig.add_trace(go.Scatter(x=t, y=infection_rate, mode='lines', name='Infection Rate (%)', 
                            line=dict(color='red', width=2)), row=1, col=2)
    fig.add_trace(go.Scatter(x=t, y=exposure_rate, mode='lines', name='Exposure Rate (%)', 
                            line=dict(color='orange', width=2)), row=1, col=2)
    
    # Seasonal transmission rate
    beta_seasonal = [model.seasonal_transmission_rate(time) for time in t]
    fig.add_trace(go.Scatter(x=t, y=beta_seasonal, mode='lines', name='Seasonal β', 
                            line=dict(color='purple', width=2)), row=2, col=1)
    
    # Phase plane
    fig.add_trace(go.Scatter(x=S, y=I, mode='lines', name='Phase Plane', 
                            line=dict(color='orange', width=2)), row=2, col=2)
    
    fig.update_layout(height=800, showlegend=True, title_text="Malaria SEIR Model Dynamics")
    
    # Update axes
    fig.update_xaxes(title_text="Time (days)", row=1, col=1)
    fig.update_yaxes(title_text="Number of individuals", row=1, col=1)
    fig.update_xaxes(title_text="Time (days)", row=1, col=2)
    fig.update_yaxes(title_text="Rate (%)", row=1, col=2)
    fig.update_xaxes(title_text="Time (days)", row=2, col=1)
    fig.update_yaxes(title_text="Transmission Rate (β)", row=2, col=1)
    fig.update_xaxes(title_text="Susceptible", row=2, col=2)
    fig.update_yaxes(title_text="Infected", row=2, col=2)
    
    st.plotly_chart(fig, use_container_width=True)

def display_age_structured_results():
    model = st.session_state['model']
    t = st.session_state['t']
    solution = st.session_state['solution']
    metrics = st.session_state['metrics']
    
    st.header("📊 Age-Structured SIR Model Results")
    
    # Age-specific metrics
    st.subheader("🎯 Age-Specific Attack Rates")
    
    age_data = []
    for i, age_group in enumerate(model.age_groups):
        age_data.append({
            'Age Group': age_group,
            'Peak Infection Rate (%)': metrics['peak_infection_rate'][i] * 100,
            'Final Attack Rate (%)': metrics['final_attack_rate'][i] * 100,
            'Time to Peak (days)': metrics['time_to_peak'][i],
            'Age-Specific R₀': metrics['R0_age'][i]
        })
    
    df_age = pd.DataFrame(age_data)
    
    # Display metrics
    col1, col2 = st.columns(2)
    
    with col1:
        st.dataframe(df_age.round(2))
    
    with col2:
        # Age-specific R₀ bar chart
        fig_r0 = px.bar(df_age, x='Age Group', y='Age-Specific R₀', 
                       title='Age-Specific Basic Reproduction Numbers',
                       color='Age-Specific R₀', color_continuous_scale='Reds')
        st.plotly_chart(fig_r0, use_container_width=True)
    
    # Age-specific infection rates over time
    st.subheader("📈 Age-Specific Infection Rates Over Time")
    
    S = solution[:, :model.n_age_groups]
    I = solution[:, model.n_age_groups:2*model.n_age_groups]
    R = solution[:, 2*model.n_age_groups:]
    
    fig_age = go.Figure()
    
    for i, age_group in enumerate(model.age_groups):
        infection_rate = I[:, i] / model.N[i] * 100
        fig_age.add_trace(go.Scatter(
            x=t, y=infection_rate, mode='lines', 
            name=f'Age {age_group}', line=dict(width=2)
        ))
    
    fig_age.update_layout(
        title='Age-Specific Infection Rates Over Time',
        xaxis_title='Time (days)',
        yaxis_title='Infection Rate (%)',
        height=500
    )
    
    st.plotly_chart(fig_age, use_container_width=True)
    
    # Contact matrix visualization
    st.subheader("🤝 Contact Matrix")
    
    contact_df = pd.DataFrame(
        model.contact_matrix,
        index=model.age_groups,
        columns=model.age_groups
    )
    
    fig_contact = px.imshow(
        contact_df, 
        title='Contact Matrix Between Age Groups',
        color_continuous_scale='Blues',
        aspect='auto'
    )
    
    st.plotly_chart(fig_contact, use_container_width=True)

def display_vector_dynamics_results():
    model = st.session_state['model']
    t = st.session_state['t']
    solution = st.session_state['solution']
    metrics = st.session_state['metrics']
    
    st.header("📊 Vector Dynamics Model Results")
    
    # Key metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Basic R₀", f"{metrics['R0_basic']:.2f}")
    with col2:
        st.metric("Peak Human Infection", f"{metrics['peak_human_infection_rate']*100:.1f}%")
    with col3:
        st.metric("Peak Vector Infection", f"{metrics['peak_vector_infection_rate']*100:.1f}%")
    with col4:
        st.metric("Vector:Human Ratio", f"{metrics['vector_human_ratio']:.1f}")
    
    # Vector dynamics explanation
    st.markdown("""
    <div class="parameter-explanation">
    <h4>📖 Vector Dynamics Model Explanation:</h4>
    <p>This model includes both human and mosquito (vector) populations, representing the 
    complete malaria transmission cycle. Mosquitoes bite infected humans, become infected, 
    and then transmit the parasite to susceptible humans.</p>
    <ul>
    <li><strong>Human-Vector Cycle</strong>: Complete transmission cycle between humans and mosquitoes</li>
    <li><strong>Seasonal Effects</strong>: Vector populations vary with climate</li>
    <li><strong>Biting Rates</strong>: Frequency of mosquito-human contact</li>
    <li><strong>Transmission Probabilities</strong>: Likelihood of transmission in each direction</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # Vector dynamics plots
    S_h, E_h, I_h, R_h, S_v, E_v, I_v = solution.T
    
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Human Population Dynamics', 'Vector Population Dynamics', 
                       'Infection Rates Comparison', 'Seasonal Vector Population')
    )
    
    # Human dynamics
    fig.add_trace(go.Scatter(x=t, y=S_h, mode='lines', name='Susceptible Humans', 
                            line=dict(color='blue', width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=E_h, mode='lines', name='Exposed Humans', 
                            line=dict(color='orange', width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=I_h, mode='lines', name='Infected Humans', 
                            line=dict(color='red', width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=R_h, mode='lines', name='Recovered Humans', 
                            line=dict(color='green', width=2)), row=1, col=1)
    
    # Vector dynamics
    fig.add_trace(go.Scatter(x=t, y=S_v, mode='lines', name='Susceptible Vectors', 
                            line=dict(color='lightblue', width=2)), row=1, col=2)
    fig.add_trace(go.Scatter(x=t, y=E_v, mode='lines', name='Exposed Vectors', 
                            line=dict(color='yellow', width=2)), row=1, col=2)
    fig.add_trace(go.Scatter(x=t, y=I_v, mode='lines', name='Infected Vectors', 
                            line=dict(color='darkred', width=2)), row=1, col=2)
    
    # Infection rates comparison
    human_infection_rate = I_h / model.N_h * 100
    vector_infection_rate = I_v / model.N_v * 100
    fig.add_trace(go.Scatter(x=t, y=human_infection_rate, mode='lines', name='Human Infection Rate (%)', 
                            line=dict(color='red', width=2)), row=2, col=1)
    fig.add_trace(go.Scatter(x=t, y=vector_infection_rate, mode='lines', name='Vector Infection Rate (%)', 
                            line=dict(color='purple', width=2)), row=2, col=1)
    
    # Seasonal vector population
    N_v_seasonal = [model.N_v * model.seasonal_vector_population(time) for time in t]
    fig.add_trace(go.Scatter(x=t, y=N_v_seasonal, mode='lines', name='Seasonal Vector Population', 
                            line=dict(color='green', width=2)), row=2, col=2)
    
    fig.update_layout(height=800, showlegend=True, title_text="Malaria Vector Dynamics Model")
    
    # Update axes
    fig.update_xaxes(title_text="Time (days)", row=1, col=1)
    fig.update_yaxes(title_text="Number of humans", row=1, col=1)
    fig.update_xaxes(title_text="Time (days)", row=1, col=2)
    fig.update_yaxes(title_text="Number of vectors", row=1, col=2)
    fig.update_xaxes(title_text="Time (days)", row=2, col=1)
    fig.update_yaxes(title_text="Infection Rate (%)", row=2, col=1)
    fig.update_xaxes(title_text="Time (days)", row=2, col=2)
    fig.update_yaxes(title_text="Vector Population", row=2, col=2)
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Intervention analysis
    st.subheader("🛡️ Intervention Analysis")
    
    if st.button("Analyze Interventions"):
        with st.spinner("Analyzing intervention effectiveness..."):
            interventions = [
                ('insecticide', 0.5, 200, 100),
                ('bed_nets', 0.7, 200, 100),
                ('larvicide', 0.6, 200, 100)
            ]
            
            intervention_results = []
            
            for intervention_type, strength, start, duration in interventions:
                effectiveness = model.intervention_analysis(
                    intervention_type, strength, start, duration
                )
                intervention_results.append({
                    'Intervention': intervention_type.title(),
                    'Strength': strength,
                    'Peak Reduction (%)': effectiveness['reduction_in_peak_infection'] * 100,
                    'Attack Rate Reduction (%)': effectiveness['reduction_in_attack_rate'] * 100
                })
            
            df_interventions = pd.DataFrame(intervention_results)
            
            # Display results
            col1, col2 = st.columns(2)
            
            with col1:
                st.dataframe(df_interventions.round(1))
            
            with col2:
                fig_interventions = px.bar(
                    df_interventions, 
                    x='Intervention', 
                    y=['Peak Reduction (%)', 'Attack Rate Reduction (%)'],
                    title='Intervention Effectiveness',
                    barmode='group'
                )
                st.plotly_chart(fig_interventions, use_container_width=True)

if __name__ == "__main__":
    main()


