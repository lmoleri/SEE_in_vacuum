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

3. Sample $N_{\rm SE}^{(i)}$ from the Poisson distribution. Mathematically this is inverse-CDF sampling; the implementation uses **ROOT's TRandom::Poisson($\mu$)** for correct Poisson variates.

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
    
    # Sample from Poisson distribution (implementation uses ROOT's TRandom::Poisson(mu))
    N_SE_i = Poisson(mu)
    total_SE += N_SE_i

# Result: total_SE is the total number of secondary electrons
```

### Implementation Details

#### 1. Poisson sampling

The script samples $N_{\rm SE} \sim \text{Poisson}(\mu)$ using **ROOT's TRandom::Poisson($\mu$)**. This ensures correct Poisson variates and consistency: the mean of the sampled $N_{\rm SE}$ over all events equals $\langle \mu \rangle = (\langle \Delta E \rangle/\epsilon)\,P_{\rm esc}$, so **Mean SEY (MC)** matches **Expected (from sampled Edep)** up to statistics. The script reports both **Expected (histogram)** (from the histogram's GetMean()) and **Expected (from sampled Edep)**; the latter must match Mean SEY.

#### 2. Handling Edge Cases

- **Very small $\Delta E$**: If $\Delta E < \epsilon$, then $N_{\rm int} < 1$, and $\mu$ is very small. Most samples will be $k=0$.
- **Very large $\Delta E$**: For large energy depositions, $\mu$ can be large. Use efficient Poisson sampling algorithms.
- **Zero energy deposition**: Skip events with $\Delta E = 0$.

#### 3. Statistical Analysis

After processing all events:
- **Mean SEY**: $\langle N_{\rm SE} \rangle = \frac{1}{n}\sum_{i=1}^n N_{\rm SE}^{(i)}$
- **Variance**: $\sigma^2 = \frac{1}{n-1}\sum_{i=1}^n (N_{\rm SE}^{(i)} - \langle N_{\rm SE} \rangle)^2$
- **Distribution**: Histogram of $N_{\rm SE}$ values

### Decomposition of results and theoretical P(≥1 SE | Edep > 0)

The script reports a **decomposition** that relates the fraction of events with at least one secondary electron (SEE) to the energy-deposition histogram and the escape probability.

#### Definition of “Edep > 0”

“Edep > 0” is defined by **excluding the histogram bin that contains 0** (the “zero bin”). The number of events with Edep > 0 is therefore the histogram integral minus the content of that bin. This matches the idea that only events that deposit energy *above* the zero bin contribute to the conditional statistics.

#### Quantities from the histogram

- **Fraction of events with Edep > 0**:  
  $f_{\rm dep} = (\text{integral over all bins} - \text{content of zero bin}) / \text{integral over all bins}$.

- **Mean Edep given Edep > 0**:  
  $\langle \Delta E \rangle_{\Delta E>0} = \frac{\sum_i n_i\,E_i}{\text{integral\_positive}}$  
  where the sum runs over all bins *except* the zero bin, $n_i$ is the bin content and $E_i$ the bin center (eV).

#### Theoretical P(≥1 SE | Edep > 0)

For a single event with energy deposition $\Delta E$, the number of secondary electrons is $N_{\rm SE} \sim \text{Poisson}(\mu)$ with $\mu = (\Delta E/\epsilon)\,\langle P_{\rm esc}\rangle_z$. So the probability of at least one SE is:

$$\Pr(N_{\rm SE} \ge 1 \mid \Delta E) = 1 - e^{-\mu} = 1 - \exp\!\left(-\frac{\Delta E}{\epsilon}\,\langle P_{\rm esc}\rangle_z\right).$$

The **theoretical P(≥1 SE | Edep > 0)** is the average of this quantity over the “Edep > 0” part of the histogram, **weighted by bin content**:

$$\Pr(N_{\rm SE} \ge 1 \mid \Delta E > 0)_{\rm theory}
= \frac{\sum_i n_i\,\bigl(1 - \exp(-(E_i/\epsilon)\,\langle P_{\rm esc}\rangle_z)\bigr)}{\text{integral\_positive}},$$

where the sum is over all bins *except* the zero bin, $n_i$ is the content of bin $i$, and $E_i$ is the bin center (eV). So each bin contributes with weight $n_i$; the numerator is the sum of “content × probability of at least one SE at that bin center”, and the denominator is the total number of events with Edep > 0.

This is exactly what the script computes and reports as “P(≥1 SE | Edep > 0) [theoretical from histogram]” (e.g. 52.7%): it loops over bins (skipping the zero bin), and for each bin adds $n_i \times (1 - \exp(-\mu_i))$ to the numerator, with $\mu_i = (E_i/\epsilon)\,P_{\rm esc}$.

#### Expected vs actual fraction with SEE

The **expected** fraction of events with at least one SEE is computed so it is directly comparable to the **actual** MC fraction (same population: all events; same definition: fraction with $N_{\rm SE} \ge 1$). Sum over *all* histogram bins (including the bin containing zero):

$$\text{expected fraction with SEE} = \sum_i \frac{n_i}{N_{\rm total}}\, \Pr(N_{\rm SE} \ge 1 \mid E_i) = \sum_i \frac{n_i}{N_{\rm total}}\, \bigl(1 - e^{-(E_i/\epsilon)\,\langle P_{\rm esc}\rangle_z}\bigr),$$

where $n_i$ is the content of bin $i$, $E_i$ is the bin center (eV), and $N_{\rm total}$ is the histogram integral. For bins with $E_i \le 0$, $\Pr(N_{\rm SE} \ge 1 \mid E_i) = 0$. This matches the MC procedure: each event is drawn from a bin (with probability $n_i/N_{\rm total}$) and then $N_{\rm SE}$ is Poisson with mean $(E_i/\epsilon)\,P_{\rm esc}$ when using bin centers, or the same in expectation when using GetRandom() within bins.

The **actual** fraction with SEE (from the Monte Carlo) is the fraction of events for which the sampled $N_{\rm SE} \ge 1$ . **Expected and actual are directly comparable** and should agree within MC sampling noise. The script reports both; small differences are due to finite statistics and (when using GetRandom()) sampling within bins rather than at bin centers.

#### Single ionization (Edep ≈ ε)

For a single “ionization” (ΔE ≈ ε = 27 eV), $\mu = \langle P_{\rm esc}\rangle_z \approx 0.38$, so:

$$\Pr(N_{\rm SE} \ge 1 \mid \Delta E \approx \epsilon) \approx 1 - e^{-0.38} \approx 32\%.$$

So **32%** is the probability of at least one SE when the deposition is about one ionization. The **theoretical P(≥1 SE | Edep > 0)** (e.g. 52.7%) is higher because it is an average over the full “Edep > 0” distribution, which has a tail to higher ΔE; those higher depositions have larger μ and hence larger $\Pr(N_{\rm SE} \ge 1)$.

## Expected Output

1. **Per-event secondary electron counts**: $N_{\rm SE}^{(i)}$ for each event
2. **Total secondary electron yield**: $N_{\rm SE}$ summed over all events
3. **Average SEY**: Mean number of secondary electrons per muon interaction
4. **Distribution histogram**: Histogram of $N_{\rm SE}$ values
5. **N_int histogram**: Distribution of $N_{\rm int} = \Delta E/\epsilon$ per event (filled once per event when using histogram sampling; Entries = number of events with Edep > 0)
6. **Edep debug plot** (when using histogram sampling): Fine-binned (0.5 eV) histogram of MC-sampled energy deposition values from EdepPrimary, saved as PDF and ROOT (`*_Edep_debug.pdf`, `*_Edep_debug.root`), for checking behavior near zero
7. **Consistency**: Mean SEY (MC) ≈ Expected (from sampled Edep); expected and actual fraction with SEE agree within MC noise. The PDF plot displays **Check: Mean SEY = Expected** and **Check: fraction with SEE** boxes so these checks are visible at a glance (see [PLOT_LEGEND_EXPLANATION.md](PLOT_LEGEND_EXPLANATION.md)).

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

2. **Python script** (current implementation):
   - Use ROOT Python bindings to read ROOT files and for Poisson sampling (TRandom::Poisson)
   - Use `numpy` for statistics
   - Poisson variates: ROOT's TRandom::Poisson($\mu$) for correct mean and consistency with Expected

3. **Integration with existing code**:
   - Could be added as a post-processing step in `RunAction::EndOfRunAction()`
   - Or as a separate analysis script

## References

- Energy per free electron: $\epsilon = 27~\text{eV}$ (conservative estimate)
- Escape probability parameters: $B = 0.46$, $\alpha = 0.0075~\text{Å}^{-1}$
- Production depth: $z = 2.5~\text{nm}$ (initial assumption)
