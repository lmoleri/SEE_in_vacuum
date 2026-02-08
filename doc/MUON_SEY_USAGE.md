# Monte Carlo SEY Calculation - Usage Guide

## Quick Start

```bash
# Basic usage with default parameters
conda run -n geant4 python calculate_muon_sey.py <input_root_file>

# Example with actual data
conda run -n geant4 python calculate_muon_sey.py \
  results/scan_thick5nm_particlemu-_energy4GeV_events10000_modelPAI/SEE_in_vacuum_thick5nm_particlemu-_energy4000MeV_events10000.root
```

## Command Line Options

```bash
python calculate_muon_sey.py <input_file> [options]
python calculate_muon_sey.py --config config/toy_model/toy_model_config.json [options]

Options:
  --config, -c CONFIG          Path to JSON config (provides input file, histogram, seed, epsilon, B, alpha, depth, bin_by_bin). CLI overrides config.
  --histogram, -H HISTOGRAM    Histogram name (default: EdepPrimary)
  --seed, -s SEED              Random seed (default: 42)
  --epsilon, -e EPSILON        Energy per free electron in eV (default: 27.0)
  --B B                        Surface escape probability (default: 0.46)
  --alpha, -a ALPHA            Attenuation coefficient in Å^-1 (default: 0.0075)
  --depth, -d DEPTH            Production depth in Å (default: 25.0 = 2.5 nm)
  --bin-by-bin                 Process bin-by-bin instead of sampling
```

## Processing Modes

## Geant4 Energy-Deposition Assumption

This workflow uses Geant4 **only for energy deposition**; secondary electrons are computed in the custom MC.
In that case, **large production cuts are desirable** because they suppress Geant4 secondary creation and leave
energy as local deposition. Step size still matters for thin layers, but production cuts do not need to be small.

### 1. Histogram Sampling (Default)

Samples energy deposition values from the histogram distribution. This treats the histogram as a probability distribution and samples from it.

**Use when**: You want to use the histogram as a distribution and generate Monte Carlo samples.

```bash
python calculate_muon_sey.py input.root
```

### 2. Bin-by-Bin Processing

Processes each histogram bin as representing actual events. For each bin, processes all events in that bin with the bin center energy.

**Use when**: You want to process the actual histogram data structure directly.

```bash
python calculate_muon_sey.py input.root --bin-by-bin
```

## Output Files

The script generates several output files:

1. **`*_SEY_MonteCarlo.root`**: ROOT file containing:
   - `SEY_MonteCarlo`: Histogram of secondary electron yield per event
   - `EdepPrimary_original`: Copy of input energy deposition histogram
   - `Summary`: Text summary with statistics

2. **`*_SEY_MonteCarlo.pdf`**: PDF plot showing:
   - SEY distribution histogram
   - **Results**: Mean SEY, Expected (histogram), Expected (from sampled Edep)
   - **Check: Mean SEY = Expected**: Validates that Poisson sampling and P_esc are correct (Mean SEY ≈ Expected from sampled Edep)
   - **Check: fraction with SEE**: Expected (all bins Σ) vs Actual (MC); directly comparable, should agree within MC noise
   - Physical parameters used

3. **`*_SEY_MonteCarlo_plot.root`**: ROOT file with the canvas (for further editing)

4. **`*_Edep_debug.pdf`** and **`*_Edep_debug.root`** (when using histogram sampling): Debug plot of MC-sampled energy deposition values from EdepPrimary, with fine binning (0.5 eV) from 0 to 250 eV to check behavior near zero. One entry per sampled event.

5. **N_int histogram** (in the PDF and ROOT output): Distribution of $N_{\rm int} = \Delta E/\epsilon$ per event. Filled **once per event** when using histogram sampling (Entries = number of events with Edep > 0), so it matches the MC run.

## Example Results

For a 4 GeV muon simulation with 10,000 events:

```
Physical Parameters:
  ε (energy per free electron) = 27.0 eV
  B (surface escape prob) = 0.46
  α (attenuation coeff) = 0.0075 Å^-1
  z (production depth) = 25.0 Å (2.5 nm)
  P_esc(z) = 0.3814

Results:
  Total events processed: 10000
  Total secondary electrons: ~360 (depends on seed)
  Mean SEY per event: ~0.036
  Expected mean (histogram): ~0.037
  Expected mean from sampled Edep: ~0.036  (must match Mean SEY)
  Standard deviation: ~0.45
  Median SEY: 0.0
  Min SEY: 0
  Max SEY: 21
  
  Mean energy deposition: 2.29 eV
  Mean free electrons per event: 0.08
```

## Understanding the Results

### Mean SEY vs Expected Mean

The **Mean SEY** is the average number of secondary electrons per event from the Monte Carlo simulation. The script reports two expected values:

- **Expected (histogram)**: $\langle \Delta E \rangle_{\rm hist}/\epsilon \times P_{\text{esc}}$, where $\langle \Delta E \rangle_{\rm hist}$ is the histogram's GetMean().
- **Expected (from sampled Edep)**: Mean of $\mu$ over the actual sampled $\Delta E$ values, i.e. the same population as the MC. **This must match Mean SEY** up to Poisson statistics; the script checks this and warns if they differ.

Poisson sampling uses **ROOT's TRandom::Poisson($\mu$)** so that Mean SEY equals Expected (from sampled Edep) in expectation. Small differences are due to:
- Poisson statistics (variance = mean)
- Finite number of events

### Decomposition: fraction with SEE and P(≥1 SE | Edep > 0)

The script reports a **decomposition** that relates the fraction of events with at least one SEE to the energy-deposition histogram:

- **Fraction of events with Edep > 0**: fraction of events that deposit energy above the histogram’s “zero bin” (events in the zero bin are excluded).
- **Mean Edep (given Edep > 0)**: average energy deposition among those events.
- **P(≥1 SE | Edep > 0) [theoretical from histogram]**: probability of at least one secondary electron, *averaged* over the “Edep > 0” part of the histogram (weighted by bin content). See [MUON_SEY_MONTE_CARLO.md](MUON_SEY_MONTE_CARLO.md) for the full formula and derivation.
- **Expected fraction with SEE** (from histogram) = \(\sum_i (n_i/N)\,\Pr(N_{\rm SE} \ge 1 \mid E_i)\) over *all* bins. This is directly comparable to the actual MC fraction (same population and definition).
- **Actual fraction with SEE** (from the Monte Carlo) is the fraction of events for which the sampled \(N_{\rm SE} \ge 1\). Expected and actual should agree within MC sampling noise; the script reports both.

For a single ionization (Edep ≈ ε = 27 eV), P(≥1 SE) ≈ 32%. The theoretical P(≥1 SE | Edep > 0) (e.g. ~53%) is higher because it averages over the full Edep > 0 distribution, which has a tail to higher depositions.

### Interpreting the Distribution

- **Median SEY = 0.0**: Most events produce zero secondary electrons (typical for low energy deposition)
- **Max SEY**: Maximum number of secondaries in a single event (rare high-energy deposition events)
- **Standard deviation**: Spread of the distribution (Poisson variance)

## Customizing Physical Parameters

### Example: Different Material Properties

```bash
# For a material with different escape probability
python calculate_muon_sey.py input.root --B 0.5 --alpha 0.008

# For different production depth
python calculate_muon_sey.py input.root --depth 30  # 3.0 nm

# For different energy per free electron
python calculate_muon_sey.py input.root --epsilon 25
```

## Validation

The script includes built-in validation:

1. **Consistency check**: Mean SEY (MC) must match **Expected (from sampled Edep)** up to Poisson statistics; the script warns if they differ by more than 10%. Poisson sampling uses ROOT's TRandom::Poisson so this holds in expectation.
2. **Expected vs actual fraction with SEE**: Expected fraction (sum over all bins) and actual MC fraction with at least one SE should agree within MC noise.
3. **Poisson statistics**: The distribution of $N_{\rm SE}$ should follow Poisson behavior.
4. **Energy scaling**: SEY should scale with energy deposition.

## Troubleshooting

### "Histogram not found"

Make sure the histogram name matches. Check available histograms:
```bash
conda run -n geant4 root -l input.root
# Then in ROOT: .ls
```

### Very low or zero SEY

- Check energy deposition values (may be too low)
- Verify physical parameters (especially ε and P_esc)
- Check if energy deposition histogram has entries

### Large differences between sampling methods

- Histogram sampling: treats histogram as distribution
- Bin-by-bin: uses actual bin structure
- Differences are expected, especially for sparse histograms

## Integration with Analysis Workflow

1. **Run Geant4 simulation** → Generate ROOT files with energy deposition
2. **Run Monte Carlo calculation** → Generate SEY distribution
3. **Compare with other models** → Use different physics models (PAI, Livermore, Penelope)
4. **Parameter studies** → Vary physical parameters to match experimental data

## Related: Toy model (multiple crossings per event)

For studies where one "event" has many shell crossings, use **`run_toy_events.py`** (see README, section "Toy events"). It samples ΔE from an EdepPrimary histogram and draws SEE per crossing from Poisson(μ(ΔE)); total SE per event is the sum over crossings.

**Config:** `config/toy_model/toy_model_config.json` — keys include `crossings_per_event` (single value or list, e.g. `[100, 200]`), `n_events`, `edep_root_file`, and the same physics parameters (epsilon, B, alpha, depth).

**Outputs:**
- **Results:** `results/toy_events_*/toy_events_SE_per_event.root`, `toy_events_summary.txt`
- **Plots:** Per crossings value: `TotalSE_per_event_{N}_crossings.pdf` and `.root` (histogram of total SE per event; each plot shows mean, std, and **efficiency** = fraction of events with ≥1 SE)
- **Summary overlay:** `TotalSE_per_event_summary.pdf` and `.root` — all crossing counts overlaid on one canvas (when more than one crossings value is in the config)

## References

See [MUON_SEY_MONTE_CARLO.md](MUON_SEY_MONTE_CARLO.md) for detailed physical model description and mathematical formulation.
