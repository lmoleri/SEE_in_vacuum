# Plot Text Explanation

## Understanding the Monte Carlo SEY Plot Statistics

### P_esc (Escape Probability)

**P_esc** is the probability that a free electron produced at depth $z$ will reach and escape the surface.

**Formula:**
$$P_{\rm esc}(z) = B \times e^{-\alpha z}$$

**Parameters:**
- $B = 0.46$: Surface escape probability (probability at $z = 0$)
- $\alpha = 0.0075$ Å⁻¹: Attenuation coefficient (describes how escape probability decreases with depth)
- $z = 25$ Å = 2.5 nm: Production depth (assumed fixed depth for all electrons)

**Example calculation:**
For $z = 2.5$ nm = 25 Å:
$$P_{\rm esc}(2.5~\text{nm}) = 0.46 \times e^{-0.0075 \times 25} = 0.46 \times e^{-0.1875} \approx 0.38$$

**Physical meaning:**
- P_esc ≈ 0.38 means that about 38% of free electrons produced at 2.5 nm depth can escape to the surface
- The remaining 62% are absorbed or lose energy before reaching the surface
- This is a key parameter determining the efficiency of secondary electron emission

### Expected (Theoretical Mean SEY)

**Expected** is the theoretical expected value of the mean secondary electron yield per event, calculated from the average energy deposition.

**Formula:**
$$\text{Expected} = \frac{\langle \Delta E \rangle}{\epsilon} \times P_{\rm esc}$$

Where:
- $\langle \Delta E \rangle$: Mean energy deposition per event (from Geant4 histogram)
- $\epsilon = 27$ eV: Average energy required to produce one free electron
- $P_{\rm esc}$: Escape probability at production depth

**Physical meaning:**
- This is the **theoretical prediction** of the mean SEY, assuming:
  - All events have the mean energy deposition
  - All free electrons are produced at the same depth (2.5 nm)
  - Escape probability is constant
- **Mean SEY** (from Monte Carlo) is the **actual simulated result** with:
  - Variable energy deposition per event
  - Poisson statistics for secondary electron emission
  - Random sampling from the distribution

**Comparison:**
- **Expected**: Theoretical value (deterministic calculation)
- **Mean SEY**: Monte Carlo result (includes statistical fluctuations)

If the Monte Carlo is working correctly, the **Mean SEY** should be close to the **Expected** value, with some statistical variation due to:
- Poisson fluctuations in secondary electron emission
- Distribution of energy depositions (not all events have the mean value)
- Finite number of events

### Example Interpretation

For a typical muon simulation:
- **Mean SEY: 0.0314**: Average number of secondary electrons per event (from Monte Carlo)
- **Expected: 0.0323**: Theoretical prediction based on mean energy deposition
- **P_esc = 0.3814**: 38.14% of free electrons escape from 2.5 nm depth

The close agreement between Mean SEY and Expected validates the Monte Carlo implementation.
