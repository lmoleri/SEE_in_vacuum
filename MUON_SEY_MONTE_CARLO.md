# Monte Carlo Calculation for Secondary Electron Emission from Muon Interactions

## Overview

This document describes the implementation of a custom Monte Carlo method to calculate the probability of extracting secondary electrons from muon interactions in a thin shell material. The method uses energy deposition data from Geant4 simulations and applies a probabilistic model to estimate secondary electron emission.

## Physical Model

### Step 1: Energy Deposition from Geant4

The Geant4 simulation provides energy deposition $\Delta E$ per muon-shell interaction event. This is the total energy deposited in the shell material by a muon and its secondaries.

### Step 2: Internal Free Electron Production

The average energy required to produce one internal free electron is:

$$\epsilon = 27~\text{eV}$$

This is a conservative (high) estimate, meaning it takes at least 27 eV to produce one free electron.

The average number of free electrons produced internally per event is:

$$N_{\rm int} = \frac{\Delta E}{\epsilon}$$

where:
- $N_{\rm int}$ = number of internally produced free electrons
- $\Delta E$ = energy deposited in the shell (from Geant4)
- $\epsilon$ = average energy per internal free electron production (27 eV)

### Step 3: Production Depth Assumption

**Initial assumption**: All free electrons are produced at a fixed depth:

$$z = 2.5~\text{nm}$$

This is a simplified model. In a more sophisticated implementation, the production depth could be sampled from a distribution based on the muon's energy loss profile.

### Step 4: Escape Probability

The probability for a free electron produced at depth $z$ to reach and escape the surface is:

$$\langle P_{\rm esc}\rangle_z = B\,e^{-\alpha z}$$

where:
- $B = 0.46$ (escape probability at surface)
- $\alpha = 0.0075~\text{Å}^{-1}$ (attenuation coefficient)
- $z$ = production depth (in Å, where 1 nm = 10 Å)

**Note**: For $z = 2.5~\text{nm} = 25~\text{Å}$:
$$\langle P_{\rm esc}\rangle_{2.5~\text{nm}} = 0.46 \times e^{-0.0075 \times 25} = 0.46 \times e^{-0.1875} \approx 0.38$$

### Step 5: Average Secondary Electron Yield per Interaction

The average (expected) number of secondary electrons emitted per muon-shell interaction is:

$$\mu^{\rm z}_{N_{\rm int}} = N_{\rm int} \times \langle P_{\rm esc}\rangle_z$$

Substituting the expressions:
$$\mu^{\rm z}_{N_{\rm int}} = \frac{\Delta E}{\epsilon} \times B\,e^{-\alpha z}$$

### Step 6: Poisson Distribution for Secondary Electron Emission

The probability of emitting exactly $k$ secondary electrons in a single muon-shell interaction follows a Poisson distribution:

$$\Pr\!\big(N_{\rm SE}^{\rm shell} = k\big) = \frac{\left(\mu^{\rm z}_{N_{\rm int}}\right)^k}{k!}\,e^{-\mu^{\rm z}_{N_{\rm int}}}$$

Substituting the full expression:
$$\Pr\!\big(N_{\rm SE}^{\rm shell} = k\big) = \frac{\left(\dfrac{\Delta E}{\epsilon}\,\langle P_{\text{esc}}\rangle_z\right)^k}{k!}\,\exp\!\left[-\dfrac{\Delta E}{\epsilon}\,\langle P_{\text{esc}}\rangle_z \right]$$

### Step 7: Monte Carlo Sampling for Individual Interactions

For the $i$-th ionizing muon interaction with the shell:

1. Calculate $\mu^{\rm z}_{N_{\rm int}}^{(i)}$ from the energy deposition $\Delta E^{(i)}$:
   $$\mu^{\rm z}_{N_{\rm int}}^{(i)} = \frac{\Delta E^{(i)}}{\epsilon} \times B\,e^{-\alpha z}$$

2. Generate a uniform random number $u \in (0,1)$

3. Sample $N_{\rm SE}^{(i)}$ from the Poisson distribution using the cumulative distribution function (CDF):
   $$N_{\rm SE}^{(i)} = k \quad \text{such that} \quad \sum_{n=0}^{k}\frac{(\mu^{\rm z}_{N_{\rm int}}^{(i)})^n}{n!}e^{-\mu^{\rm z}_{N_{\rm int}}^{(i)}} \ge u$$

   This means finding the smallest integer $k$ where the cumulative probability exceeds $u$.

### Step 8: Total Secondary Electron Yield

For $n$ muon-shell interactions, the total number of secondary electrons emitted is:

$$N_{\rm SE} = \sum_{i=1}^n N_{\rm SE}^{(i)}$$

## Implementation Steps

### Input Data

1. **Energy deposition histogram**: From Geant4 simulation (`EdepPrimary` histogram)
   - Contains energy deposition values $\Delta E^{(i)}$ for each event
   - Each entry represents one muon-shell interaction

2. **Physical parameters**:
   - $\epsilon = 27~\text{eV}$ (energy per free electron)
   - $B = 0.46$ (surface escape probability)
   - $\alpha = 0.0075~\text{Å}^{-1}$ (attenuation coefficient)
   - $z = 2.5~\text{nm} = 25~\text{Å}$ (production depth)

### Algorithm

```python
# Pseudocode for Monte Carlo calculation

# Constants
epsilon = 27.0  # eV
B = 0.46
alpha = 0.0075  # Å^-1
z = 25.0  # Å (2.5 nm)

# Calculate escape probability at depth z
P_esc = B * exp(-alpha * z)

# Load energy deposition data from Geant4 ROOT file
# For each event i:
total_SE = 0
for event in events:
    delta_E = event.energy_deposition  # eV
    
    # Calculate average number of free electrons
    N_int = delta_E / epsilon
    
    # Calculate Poisson parameter
    mu = N_int * P_esc
    
    # Sample from Poisson distribution
    u = random.uniform(0, 1)
    k = 0
    cumulative = exp(-mu)  # P(k=0)
    
    while cumulative < u:
        k += 1
        cumulative += (mu**k / factorial(k)) * exp(-mu)
    
    N_SE_i = k
    total_SE += N_SE_i

# Result: total_SE is the total number of secondary electrons
```

### Implementation Details

#### 1. Poisson CDF Sampling

The cumulative distribution function for Poisson is:
$$F(k) = \sum_{n=0}^{k} \frac{\lambda^n}{n!} e^{-\lambda}$$

where $\lambda = \mu^{\rm z}_{N_{\rm int}}^{(i)}$.

To sample $k$ given random $u$:
- Start with $k=0$, $F(0) = e^{-\lambda}$
- Increment $k$ and add terms until $F(k) \ge u$
- Return $k$ as the sampled value

**Optimization**: For large $\lambda$, use normal approximation or other efficient methods.

#### 2. Handling Edge Cases

- **Very small $\Delta E$**: If $\Delta E < \epsilon$, then $N_{\rm int} < 1$, and $\mu$ is very small. Most samples will be $k=0$.
- **Very large $\Delta E$**: For large energy depositions, $\mu$ can be large. Use efficient Poisson sampling algorithms.
- **Zero energy deposition**: Skip events with $\Delta E = 0$.

#### 3. Statistical Analysis

After processing all events:
- **Mean SEY**: $\langle N_{\rm SE} \rangle = \frac{1}{n}\sum_{i=1}^n N_{\rm SE}^{(i)}$
- **Variance**: $\sigma^2 = \frac{1}{n-1}\sum_{i=1}^n (N_{\rm SE}^{(i)} - \langle N_{\rm SE} \rangle)^2$
- **Distribution**: Histogram of $N_{\rm SE}$ values

## Expected Output

1. **Per-event secondary electron counts**: $N_{\rm SE}^{(i)}$ for each event
2. **Total secondary electron yield**: $N_{\rm SE}$ summed over all events
3. **Average SEY**: Mean number of secondary electrons per muon interaction
4. **Distribution histogram**: Histogram of $N_{\rm SE}$ values
5. **Comparison with Geant4**: Compare Monte Carlo results with direct Geant4 tracking (if available)

## Future Enhancements

1. **Depth-dependent production**: Sample production depth $z$ from energy loss profile instead of fixed value
2. **Energy-dependent parameters**: Make $\epsilon$, $B$, and $\alpha$ functions of electron energy
3. **Multiple depth layers**: Consider electrons produced at different depths with different escape probabilities
4. **Angular distribution**: Include angular dependence of escape probability
5. **Material-specific parameters**: Adjust parameters based on material properties (Al2O3, SiO2, etc.)

## Validation

1. **Consistency check**: Verify that $\langle N_{\rm SE} \rangle \approx \mu^{\rm z}_{N_{\rm int}}$ averaged over all events
2. **Poisson statistics**: Check that the distribution of $N_{\rm SE}$ follows expected Poisson behavior
3. **Energy scaling**: Verify that SEY scales approximately linearly with energy deposition for small $\mu$
4. **Comparison with literature**: Compare results with experimental or theoretical values

## Code Structure

### Recommended Implementation

1. **ROOT macro** (`calculate_muon_sey.C`):
   - Read energy deposition histogram from Geant4 output
   - Implement Monte Carlo sampling
   - Generate output histograms and statistics

2. **Python script** (alternative):
   - Use `uproot` or `ROOT` Python bindings to read ROOT files
   - Use `numpy` for random number generation and statistics
   - Use `scipy.stats` for Poisson distribution functions

3. **Integration with existing code**:
   - Could be added as a post-processing step in `RunAction::EndOfRunAction()`
   - Or as a separate analysis script

## References

- Energy per free electron: $\epsilon = 27~\text{eV}$ (conservative estimate)
- Escape probability parameters: $B = 0.46$, $\alpha = 0.0075~\text{Å}^{-1}$
- Production depth: $z = 2.5~\text{nm}$ (initial assumption)
