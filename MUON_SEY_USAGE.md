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

Options:
  --histogram, -H HISTOGRAM    Histogram name (default: EdepPrimary)
  --seed, -s SEED              Random seed (default: 42)
  --epsilon, -e EPSILON        Energy per free electron in eV (default: 27.0)
  --B B                        Surface escape probability (default: 0.46)
  --alpha, -a ALPHA            Attenuation coefficient in Å^-1 (default: 0.0075)
  --depth, -d DEPTH            Production depth in Å (default: 25.0 = 2.5 nm)
  --bin-by-bin                 Process bin-by-bin instead of sampling
```

## Processing Modes

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
   - Statistics (mean, std, expected value)
   - Physical parameters used

3. **`*_SEY_MonteCarlo_plot.root`**: ROOT file with the canvas (for further editing)

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
  Total secondary electrons: 805
  Mean SEY per event: 0.0805
  Expected mean (theoretical): 0.0323
  Standard deviation: 0.7121
  Median SEY: 0.0
  Min SEY: 0
  Max SEY: 21
  
  Mean energy deposition: 2.29 eV
  Mean free electrons per event: 0.08
```

## Understanding the Results

### Mean SEY vs Expected Mean

The **Mean SEY** is the average number of secondary electrons per event from the Monte Carlo simulation. The **Expected mean (theoretical)** is calculated as:

$$\text{Expected mean} = \frac{\langle \Delta E \rangle}{\epsilon} \times P_{\text{esc}}$$

where $\langle \Delta E \rangle$ is the mean energy deposition from the histogram.

**Note**: The Monte Carlo mean may differ from the expected mean due to:
- Poisson statistics (variance = mean)
- Sampling method (histogram sampling vs bin-by-bin)
- Finite number of events

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

1. **Consistency check**: Compares Monte Carlo mean with theoretical expected mean
2. **Poisson statistics**: The distribution should follow Poisson behavior
3. **Energy scaling**: SEY should scale with energy deposition

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

## References

See `MUON_SEY_MONTE_CARLO.md` for detailed physical model description and mathematical formulation.
