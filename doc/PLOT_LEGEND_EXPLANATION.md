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

The plot may show two expected values:

- **Expected (histogram)**: $\langle \Delta E \rangle_{\rm hist}/\epsilon \times P_{\rm esc}$, where $\langle \Delta E \rangle_{\rm hist}$ is the histogram's GetMean(). This is the theoretical mean SEY from the histogram distribution.
- **Expected (from sampled Edep)**: Mean of $\mu$ over the actual MC-sampled $\Delta E$ values. **This must match Mean SEY** up to Poisson statistics; the script uses ROOT's Poisson so they agree in expectation.

**Formula (both):**
$$\text{Expected} = \frac{\langle \Delta E \rangle}{\epsilon} \times P_{\rm esc}$$

**Physical meaning:**
- **Mean SEY** (from Monte Carlo) is the average number of secondary electrons per event.
- **Expected (from sampled Edep)** uses the same population as the MC (the sampled $\Delta E$ values), so Mean SEY and this Expected should be close; small differences are Poisson fluctuations.
- **Expected (histogram)** uses the histogram's mean; it can differ slightly from Mean SEY if GetRandom() sampling and GetMean() differ (e.g. binning).

**Comparison:**
- **Mean SEY**: Monte Carlo result (includes statistical fluctuations)
- **Expected (from sampled Edep)**: Must match Mean SEY (validates correct Poisson sampling and $P_{\rm esc}$ application)
- **Expected (histogram)**: Theoretical from histogram; close to Mean SEY when sampling matches the histogram

### Example Interpretation

For a typical muon simulation:
- **Mean SEY: 0.0363**: Average number of secondary electrons per event (from Monte Carlo)
- **Expected (from sampled Edep): 0.0358**: Must match Mean SEY; close agreement validates the implementation
- **Expected (histogram): 0.0367**: From histogram GetMean()
- **P_esc = 0.3814**: 38.14% of free electrons escape from 2.5 nm depth

The close agreement between Mean SEY and Expected (from sampled Edep) validates the Monte Carlo implementation.

### Check: fraction with SEE

The plot shows a **Check: fraction with SEE** box with two numbers:

- **Expected (all bins Σ)**: Theoretical fraction of events with at least one secondary electron, computed as \(\sum_i (n_i/N)\,\Pr(N_{\rm SE} \ge 1 \mid E_i)\) over *all* histogram bins (including zero). Same population and definition as the MC.
- **Actual (MC)**: Fraction of Monte Carlo events for which the sampled \(N_{\rm SE} \ge 1\).

Expected and Actual are directly comparable and should agree within MC sampling noise. See [MUON_SEY_MONTE_CARLO.md](MUON_SEY_MONTE_CARLO.md) for the formula and derivation.
