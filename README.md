# GEANT4 SEE in Vacuum Simulation

This GEANT4 simulation models a primary particle gun (electron or muon) shooting at an Al2O3
(aluminum oxide) layer. Further documentation (Monte Carlo model, usage, plot explanations) is in the **`doc/`** folder.

## Geometry

- **Target**: Al2O3 layer
  - Thickness: 20 nm
  - Diameter: 100 nm
  - Shape: Cylindrical disk

- **Primary Gun**:
  - Particle: Electron (e-) or Muon (mu-)
  - Default energy: 1 MeV
  - Default position: 1 micron above the Al2O3 layer
  - Default direction: Downward along +z axis

## Requirements

- Geant4 (via conda environment `geant4`)
- ROOT (installed in the same conda environment)

## Building the Project

1. Create a build directory:
```bash
mkdir build
cd build
```

2. Configure with CMake:
```bash
conda activate geant4
cmake ..
```

3. Compile:
```bash
make
```

## Running the Simulation

The simulation can be run in several ways:

### 1. Interactive Mode (with visualization):
```bash
./SEE_in_vacuum
```
Opens the Geant4 visualization interface for interactive exploration.

### 2. Batch Mode (with macro file):
```bash
./SEE_in_vacuum ../run.mac
```
Runs a single simulation using a Geant4 macro file.

### 3. Parametric Scan (JSON)

Configuration files live in `config/`: **Geant4 parametric scans** use `config/geant4/*.json`; **toy model** scripts use `config/toy_model/*.json`. Create or edit a JSON file there with arrays of thickness and energy values, then pass it to the executable:

```bash
./SEE_in_vacuum ../config/geant4/scan.json
```

Example `config/geant4/scan.json`:
```json
{
  "sample_thickness_nm": [10, 20, 50],
  "primary_energy_MeV": [0.5, 1.0, 2.0],
  "em_model": "PAI",
  "primary_particle": "e-",
  "pai_enabled": true,
  "livermore_atomic_deexcitation": true,
  "events": 100000,
  "output_dir": "results/scan_thick10-20-50nm_particlee-_energy0p5-1-2MeV_events100000"
}
```

`primary_particle` can be `"e-"` or `"mu-"` (default is `"e-"`).

`em_model` can be:
- `"PAI"` (default): uses `G4EmStandardPhysics_option4` + PAI in Al2O3
- `"G4EmLivermorePhysics"` (or `"livermore"`)
- `"G4EmPenelopePhysics"` (or `"penelope"`)

**Important: Model applicability by particle type:**
- **PAI model**: Applied to both **electrons** and **muons** in the Al2O3 region. The PAI (Photoabsorption Ionization) model is manually attached to electron and muon ionization processes.
- **Livermore model**: Only applies to **electrons and photons**. For muons, Livermore physics constructors fall back to standard GEANT4 muon models (Bethe-Bloch based). This means Livermore and standard physics will produce identical results for muons.
- **Penelope model**: Only applies to **electrons and photons**. For muons, Penelope physics constructors fall back to standard GEANT4 muon models (Bethe-Bloch based). This means Penelope and standard physics will produce identical results for muons.

**Note**: When comparing muon simulations, only the PAI model will show differences from standard physics. Livermore and Penelope muon results will be identical to each other and to standard GEANT4 muon physics, as they all use the same underlying Bethe-Bloch models for muon ionization.

`pai_enabled` is only used when `em_model` is `"PAI"`.

`livermore_atomic_deexcitation` is only used when `em_model` is `"G4EmLivermorePhysics"`
(or `"livermore"`), and toggles atomic deexcitation (Fluo/Auger/PIXE).

Output files are created inside `output_dir` for each combination. When `output_dir`
is relative, it is resolved from the project root (not `build/`), e.g.:
`results/scan_thick10-20-50nm_particlee-_energy0p5-1-2MeV_events100000/SEE_in_vacuum_thick20nm_particlee-_energy1MeV_events100000.root`

## Outputs

- `SEE_in_vacuum.root` (ROOT file with histograms and canvases)
  - `EdepPrimary`: total energy deposition in Al2O3 (eV) from the primary particle and all its descendants (secondaries, delta rays, etc.). This represents the total energy that originated from the primary particle. The histogram is created with fine binning (1 eV per bin) and a reasonable range (0-10000 eV), then automatically optimized after the simulation to trim the range after the last non-empty bin, ensuring proper x-axis range and binning that matches the data distribution.
  - `EdepInteractions`: number of energy-depositing steps in Al2O3 per event
  - `EdepStep`: energy deposition per step in Al2O3 (eV) from any particle (primary or secondary)
  - `PAITransfer`: per-step energy transfer proxy in Al2O3 (eV)
  - `PrimaryResidualEnergy`: primary residual kinetic energy at end of event (eV)
  - `PrimaryEndVolume`: end volume category (Al2O3/World/OutOfWorld/Other)
  - `StepLengthAl2O3`: step length distribution in Al2O3 (nm)
  - `ResidualEnergyVsEndVolume`: 2D residual energy vs end volume category
  - `ResidualEnergyVsLastProcess`: 2D residual energy vs last process category
  - `ResidualEnergyVsStopStatus`: 2D residual energy vs stop status category
  - `RunMeta`: ntuple with `primaryEnergyMeV`, `sampleThicknessNm`,
  `maxPrimaryEnergyMeV`, `paiEnabled`, `primaryParticle`, `emModel`,
  and `livermoreAtomicDeexcitation`
  - `EdepPrimaryCanvas`, `EdepInteractionsCanvas`, and other canvases saved with annotations

## Plotting

Run the ROOT macro to draw and save canvases (with annotations) into the ROOT file. From the project root:

```bash
root -l draw_histo.C
```

When run on scan outputs, canvases are also exported to `plots/` as `.root` and `.pdf`,
mirroring the `results/` subfolder structure.

To plot a specific file from a parametric scan:
```bash
root -l 'draw_histo.C("SEE_in_vacuum_thick20nm_energy1MeV.root")'
```

To process all scan outputs at once:
```bash
bash scripts/run_draw_all.sh results/scan_thick5-10-15-20-25nm_energy1MeV_events100000
```

To create summary ROOT files with overlapping histograms and legends (one per particle):
```bash
bash scripts/run_draw_summary.sh results/scan_thick5-10-15-20-25nm_energy1MeV_events100000
```

To compare models per energy (one canvas per energy value):
```bash
bash scripts/run_draw_summary_models.sh \
  results/scan_*_modelPAI \
  results/scan_*_modelLivermore \
  results/scan_*_modelPenelope \
  results/summary_models.root
```

Model comparison canvases are also exported to `plots/` alongside the summary root file.

### 4. Geometry Toy Model (Shell Crossings)

This repository also includes a geometry-only toy model for shell crossings. It supports analytical scans and a 2D slice Monte Carlo (circles cut from 3D spheres).

**Config file:** `config/geometry/geometry_config.json`

**Key fields:**
- `reference`: baseline parameters (`L_um`, `a_nm`, `d_nm`, `phi`)
- `scan_ranges`: arrays for `phi`, `d_nm`, `a_nm`, `L_um`
- `montecarlo_2d`: MC settings (`W_um`, `n_rays`, `seed`, `scan_slice_plots`)

**Run analytical scans:**
```bash
conda run -n geant4 python geometry_analytical.py
```

**Run 2D MC scans:**
```bash
conda run -n geant4 python geometry_montecarlo_2d.py --scan
```

**Regenerate plots without re-running MC:**
```bash
conda run -n geant4 python geometry_montecarlo_2d.py --plot-only
```

**Outputs:**
- `results/geometry/geometry_mc_2d_scan_summary.json`: cached scan results (used by `--plot-only`)
- `plots/geometry/mc_2d_scan_*.pdf|root`: scan summaries
- `plots/geometry/slice_configs_2d/*_zoom.pdf|root`: zoomed slice configurations
- `plots/geometry/distributions_2d/crossings_dist_*.pdf|root`: per-scan crossing histograms
- `plots/geometry/crossings_vs_core_radius_mc_vs_ana_2d.pdf|root`: 2D MC vs analytical comparison

**Filling factor convention:** The MC uses the 3D volume filling factor `phi` to set sphere density. The achieved value is estimated from the slice area (Delesse principle) and is always reported as a **volume** filling factor. Plots annotate this as `Achieved #phi (vol)` and the `phi` scan uses achieved `phi` on the x-axis.

See `doc/GEOMETRY_2D_MC.md` for more details.

### 5. Monte Carlo Post-Processing (Toy Model for Muon SEY)

For muon simulations, a custom Monte Carlo method can be applied to calculate secondary electron emission (SEY) from energy deposition data. This post-processing step uses a probabilistic model to estimate SEY based on:

- Energy deposition per event (from Geant4)
- Physical parameters (energy per free electron, escape probability, production depth)
- Poisson statistics for secondary electron emission

**Usage:**
```bash
conda run -n geant4 python calculate_muon_sey.py <input_root_file> [options]
conda run -n geant4 python calculate_muon_sey.py --config config/toy_model/toy_model_config.json [options]
```

**Example:**
```bash
conda run -n geant4 python calculate_muon_sey.py \
  results/scan_thick5nm_particlemu-_energy4GeV_events10000_modelPAI/SEE_in_vacuum_thick5nm_particlemu-_energy4000MeV_events10000.root \
  --bin-by-bin
conda run -n geant4 python calculate_muon_sey.py --config config/toy_model/toy_model_config.json --bin-by-bin
```

**Options:**
- `--config, -c`: Path to toy model config JSON (provides input file as `edep_root_file` or `input_file`, plus histogram, seed, epsilon, B, alpha, depth, bin_by_bin). CLI overrides config.
- `--histogram, -H`: Histogram name (default: `EdepPrimary` or from config)
- `--seed, -s`: Random seed (default: 42 or from config)
- `--epsilon, -e`: Energy per free electron in eV (default: 27.0 or from config)
- `--B`: Surface escape probability (default: 0.46 or from config)
- `--alpha, -a`: Attenuation coefficient in Å⁻¹ (default: 0.0075 or from config)
- `--depth, -d`: Production depth in Å (default: 25.0 = 2.5 nm, or from config)
- `--bin-by-bin`: Process histogram bin-by-bin instead of sampling (overrides config if set)

**Output:**
- ROOT file: `*_SEY_MonteCarlo.root` (contains SEY histogram, N_int histogram, and statistics)
- PDF plot: Saved in `plots/` folder with logarithmic Y-axis; includes **Results** (Mean SEY, Expected values), **Check: Mean SEY = Expected** (validates Poisson sampling), and **Check: fraction with SEE** (Expected vs Actual, comparable over all events)
- Edep debug: `*_Edep_debug.pdf` and `*_Edep_debug.root` (histogram sampling only) — MC-sampled energy deposition with fine binning near zero
- Statistics: Mean SEY, expected (histogram and from sampled Edep), standard deviation; consistency checks that Mean SEY ≈ Expected (from sampled Edep) and expected vs actual fraction with SEE

**Processing Modes:**
- **Histogram sampling** (default): Samples energy deposition from the histogram distribution
- **Bin-by-bin** (`--bin-by-bin`): Processes each histogram bin as actual events (more accurate)

**Physical Model:**
The Monte Carlo method implements:
- Internal free electron production: $N_{\rm int} = \Delta E / \epsilon$ (where $\epsilon = 27$ eV)
- Escape probability: $P_{\rm esc} = B e^{-\alpha z}$ (at depth $z = 2.5$ nm)
- Poisson sampling for secondary electron emission per event
- Total SEY calculation across all events

See [doc/MUON_SEY_MONTE_CARLO.md](doc/MUON_SEY_MONTE_CARLO.md) for detailed documentation of the physical model and implementation.

**Workflow:**
1. Run Geant4 muon simulation (e.g., with PAI model)
2. Extract energy deposition histogram (`EdepPrimary`)
3. Apply Monte Carlo post-processing to calculate SEY
4. Compare results with experimental data or other models

**6. Toy events (many crossings per event)**

For studies where one "event" has many primary crossings of the shell: energy depositions are sampled from an EdepPrimary histogram (e.g. from a 100k Geant4 run), and SEE per crossing is drawn from Poisson(μ(ΔE)). Total SE per event is the sum over crossings. All parameters are set via a JSON config file.

**Config file:** `config/toy_model/toy_model_config.json` (or path passed as first argument). Keys: `edep_root_file`, `histogram`, `n_events`, `crossings_per_event` (single value or list for scan), `seed`, `epsilon`, `B`, `alpha`, `depth`, `output_dir`.

**Run:**
```bash
conda run -n geant4 python run_toy_events.py
```
(Defaults to `config/toy_model/toy_model_config.json`; or pass a path, e.g. `config/toy_model/toy_model_config.json`.)

**Example:** 4 GeV muons, 100 and 200 crossings per event, 10k events — config uses `edep_root_file` pointing to the 100k-event ROOT file and `crossings_per_event: [100, 200]`.

**Output:**
- **Results:** `results/toy_events_4GeV_muons/toy_events_SE_per_event.root` (histograms of total SE per event per crossings value), `toy_events_summary.txt`, `toy_model_config_used.json`.
- **Plots (per crossings value):** `plots/toy_events_4GeV_muons/TotalSE_per_event_{N}_crossings.pdf` and `.root` — histogram of total SE per event, with mean, std, and **efficiency** (fraction of events with ≥1 SE).
- **Summary overlay:** `TotalSE_per_event_summary.pdf` and `.root` — all crossing counts overlaid on one canvas (when `crossings_per_event` has more than one value).
- **Text summary** includes mean, std, min, max, and efficiency per crossings value.

Notes:
- The energy deposition plot uses log Y scale by default.
- Annotations read `RunMeta` from the ROOT file.
- Plot titles/legends include the selected EM model (PAI/Livermore/Penelope).
- EM low-energy cutoffs are set to 0.1 eV (see `src/PhysicsList.cc`).

## Simulation Methods Summary

| Method | Use Case | Input | Output |
|--------|----------|-------|--------|
| **Interactive** | Exploration, debugging | None (GUI) | Single ROOT file |
| **Batch (Macro)** | Single simulation run | `.mac` file | Single ROOT file |
| **Parametric Scan** | Systematic studies | `config/geant4/*.json` | Multiple ROOT files |
| **Monte Carlo Post-Processing** | Muon SEY calculation | ROOT file from scan | SEY plots, N_int and Edep debug, statistics |
| **Toy events** | Many crossings per event | `config/toy_model/*.json` | SE per event histograms, summary overlay, efficiency |

## Customization

You can modify the electron gun parameters in `src/PrimaryGeneratorAction.cc`:
- Energy: Change `SetParticleEnergy()`
- Position: Change `SetParticlePosition()`
- Direction: Change `SetParticleMomentumDirection()`

You can also modify the Al2O3 geometry in `src/DetectorConstruction.cc`:
- Thickness: Change the `thickness` variable (or use JSON config)
- Diameter: Change the `radius` variable

For parametric scans, all parameters are configurable via JSON in `config/geant4/` (see example above).

## Example Macro File (run.mac)

Create a `run.mac` file to run simulations:

```
# Set number of events
/run/beamOn 100000

# Optional: Change electron energy
/gun/energy 5 MeV

# Optional: Change electron position
/gun/position 0 0 -1 um

# Optional: Change electron direction
/gun/direction 0 0 1
```

## Additional Documentation

- **`doc/MUON_SEY_MONTE_CARLO.md`**: Detailed documentation of the Monte Carlo physical model and implementation for muon secondary electron emission
- **`doc/MUON_SEY_USAGE.md`**: Usage guide for the Monte Carlo SEY calculation script with examples and troubleshooting
- **`doc/PLOT_LEGEND_EXPLANATION.md`**: Explanation of P_esc and Expected values on the SEY plot
