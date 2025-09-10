# 🧪 Experimental Setups for Malaria Model Simulation

This document provides detailed experimental setups and scenarios you can explore using the web interface.

## 🚀 Quick Start

1. **Start the web interface:**
   ```bash
   ./start_web_interface.sh
   ```

2. **Open your browser to:** `http://localhost:8501`

3. **Choose a model type** from the sidebar dropdown

4. **Adjust parameters** using the sliders and inputs

5. **Click "Run Simulation"** to see results

## 🧬 Model Types & Experimental Setups

### 1. Basic SIR Model Experiments

#### 🎯 **Experiment 1: Understanding R₀ (Basic Reproduction Number)**

**Objective:** Understand how transmission rate affects epidemic potential

**Setup:**
- **Transmission Rate (β):** Start with 0.2, then try 0.4, 0.6, 0.8
- **Recovery Rate (γ):** Keep at 0.05 (20-day infectious period)
- **Population Size:** 10,000
- **Initial Infected:** 50

**What to Observe:**
- How R₀ changes with transmission rate
- Peak infection rates
- Time to epidemic peak
- Final attack rates

**Expected Results:**
- R₀ < 1: No epidemic (disease dies out)
- R₀ > 1: Epidemic occurs
- Higher R₀ = faster, more severe epidemics

#### 🎯 **Experiment 2: Recovery Rate Impact**

**Objective:** Study how treatment affects epidemic dynamics

**Setup:**
- **Transmission Rate (β):** Keep at 0.4
- **Recovery Rate (γ):** Try 0.02, 0.05, 0.1, 0.2
- **Population Size:** 10,000
- **Initial Infected:** 50

**What to Observe:**
- How faster recovery (better treatment) affects epidemic size
- Changes in peak timing
- Impact on final attack rates

**Expected Results:**
- Faster recovery = smaller epidemics
- Better treatment = lower peak infections

#### 🎯 **Experiment 3: Population Size Effects**

**Objective:** Understand how population size affects epidemic patterns

**Setup:**
- **Transmission Rate (β):** 0.4
- **Recovery Rate (γ):** 0.05
- **Population Size:** Try 1,000, 5,000, 10,000, 50,000
- **Initial Infected:** 50

**What to Observe:**
- Epidemic curves at different population sizes
- Peak infection rates
- Time to peak

**Expected Results:**
- Larger populations = longer epidemics
- Similar peak infection rates (as percentages)

### 2. SEIR Model Experiments

#### 🎯 **Experiment 4: Latent Period Effects**

**Objective:** Study how incubation period affects epidemic dynamics

**Setup:**
- **Transmission Rate (β):** 0.4
- **Incubation Rate (σ):** Try 0.05, 0.1, 0.2, 0.3
- **Recovery Rate (γ):** 0.05
- **Population Size:** 10,000
- **Initial Exposed:** 20
- **Initial Infected:** 30

**What to Observe:**
- How latent period affects epidemic timing
- Differences between Basic R₀ and Effective R₀
- Peak timing and magnitude

**Expected Results:**
- Longer latent period = delayed epidemic peak
- Effective R₀ accounts for latent period

#### 🎯 **Experiment 5: Seasonal Transmission**

**Objective:** Explore seasonal malaria patterns

**Setup:**
- **Transmission Rate (β):** 0.4
- **Incubation Rate (σ):** 0.1
- **Recovery Rate (γ):** 0.05
- **Seasonal Amplitude:** Try 0.0, 0.2, 0.4, 0.6
- **Simulation Duration:** 2000 days (to see multiple seasons)

**What to Observe:**
- Seasonal patterns in transmission
- Annual epidemic cycles
- Peak timing relative to seasons

**Expected Results:**
- Higher seasonal amplitude = more pronounced seasonal patterns
- Annual epidemic cycles

### 3. Age-Structured Model Experiments

#### 🎯 **Experiment 6: Age-Specific Transmission**

**Objective:** Study how different age groups affect transmission

**Setup:**
- **Total Population:** 100,000
- **Initial Infected by Age:**
  - 0-4 years: 20
  - 5-14 years: 40
  - 15-29 years: 60
  - 30-44 years: 40
  - 45-59 years: 20
  - 60+ years: 20

**What to Observe:**
- Age-specific attack rates
- Contact matrix effects
- Age-specific R₀ values

**Expected Results:**
- Children typically have higher attack rates
- Age-specific R₀ varies by age group

#### 🎯 **Experiment 7: Age-Targeted Interventions**

**Objective:** Test age-specific intervention strategies

**Setup:**
- **Scenario A:** High initial infection in children (0-4: 50, others: 5)
- **Scenario B:** High initial infection in adults (15-29: 50, others: 5)
- **Scenario C:** Uniform initial infection (all ages: 20)

**What to Observe:**
- How initial age distribution affects overall epidemic
- Age-specific epidemic curves
- Final attack rates by age

**Expected Results:**
- Children as initial source = faster epidemic spread
- Adults as initial source = different epidemic pattern

### 4. Vector Dynamics Model Experiments

#### 🎯 **Experiment 8: Vector-Human Ratio Effects**

**Objective:** Study how mosquito density affects transmission

**Setup:**
- **Human Population:** 10,000
- **Vector Population:** Try 10,000, 25,000, 50,000, 100,000
- **Biting Rate (a):** 0.3
- **Transmission Probabilities:** b=0.5, c=0.5

**What to Observe:**
- How vector density affects R₀
- Human infection rates
- Vector infection rates

**Expected Results:**
- Higher vector density = higher R₀
- More mosquitoes = more human infections

#### 🎯 **Experiment 9: Biting Rate Impact**

**Objective:** Study how mosquito behavior affects transmission

**Setup:**
- **Human Population:** 10,000
- **Vector Population:** 50,000
- **Biting Rate (a):** Try 0.1, 0.2, 0.3, 0.5
- **Transmission Probabilities:** b=0.5, c=0.5

**What to Observe:**
- R₀ changes with biting rate
- Epidemic dynamics
- Peak infection rates

**Expected Results:**
- Higher biting rate = higher transmission
- Biting rate has quadratic effect on R₀

#### 🎯 **Experiment 10: Intervention Effectiveness**

**Objective:** Compare different vector control strategies

**Setup:**
- **Baseline:** Standard parameters
- **Intervention A:** Insecticide (50% vector mortality increase)
- **Intervention B:** Bed nets (70% biting rate reduction)
- **Intervention C:** Larvicide (60% seasonal amplitude reduction)

**What to Observe:**
- Reduction in peak infections
- Reduction in final attack rates
- Relative effectiveness of interventions

**Expected Results:**
- All interventions reduce transmission
- Different interventions have different effectiveness
- Combined interventions may be most effective

## 📊 Advanced Analysis Features

### Parameter Sensitivity Analysis

1. **Run the Basic SIR model**
2. **Click "Run Sensitivity Analysis"**
3. **Observe how different transmission rates affect:**
   - R₀ values
   - Peak infection rates
   - Final attack rates

### Model Comparison

1. **Run the same scenario with different models:**
   - Basic SIR
   - SEIR
   - Age-Structured SIR
   - Vector Dynamics

2. **Compare results:**
   - Epidemic curves
   - Peak timing
   - Attack rates
   - R₀ values

### Interactive Visualizations

- **Hover over plots** for detailed values
- **Zoom and pan** on interactive charts
- **Toggle data series** on/off
- **Download plots** as images

## 🎓 Educational Scenarios

### Scenario 1: Malaria Outbreak in a Village
- **Population:** 5,000
- **Initial infected:** 10
- **High transmission rate:** 0.6
- **Study:** How quickly does the outbreak spread?

### Scenario 2: Seasonal Malaria in a Region
- **Population:** 50,000
- **High seasonal variation:** 0.7
- **Study:** Annual epidemic patterns

### Scenario 3: Age-Specific Malaria Risk
- **Focus on children:** High initial infection in 0-4 age group
- **Study:** How age affects transmission patterns

### Scenario 4: Vector Control Program
- **Before intervention:** High vector population
- **After intervention:** Reduced vector population
- **Study:** Effectiveness of control measures

## 🔬 Research Questions to Explore

1. **What is the minimum R₀ needed for an epidemic?**
2. **How does population size affect epidemic duration?**
3. **Which age groups are most important for transmission?**
4. **What is the most effective intervention strategy?**
5. **How do seasonal patterns affect long-term transmission?**
6. **What is the relationship between vector density and human infection?**
7. **How does latent period affect epidemic timing?**
8. **What are the trade-offs between different control strategies?**

## 💡 Tips for Effective Experimentation

1. **Start simple:** Begin with Basic SIR model
2. **Change one parameter at a time:** Understand individual effects
3. **Use realistic ranges:** Check parameter ranges in `data/malaria_parameters.csv`
4. **Compare scenarios:** Run multiple simulations with different parameters
5. **Take notes:** Document your observations and insights
6. **Explore edge cases:** Try extreme parameter values
7. **Use sensitivity analysis:** Understand parameter importance

## 📚 Further Reading

- **Malaria Epidemiology:** WHO Malaria Reports
- **Mathematical Modeling:** Anderson & May (1991) "Infectious Diseases of Humans"
- **Vector Dynamics:** Smith et al. (2012) "Ross, Macdonald, and a Theory for the Dynamics and Control of Mosquito-Transmitted Pathogens"

Happy experimenting! 🧬🦟📊


