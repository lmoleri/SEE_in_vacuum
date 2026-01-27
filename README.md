# GEANT4 SEE in Vacuum Simulation

This GEANT4 simulation models a primary particle gun (electron or muon) shooting at an Al2O3
(aluminum oxide) layer.

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

### Interactive Mode (with visualization):
```bash
./SEE_in_vacuum
```

### Batch Mode:
```bash
./SEE_in_vacuum ../run.mac
```

### Parametric Scan (JSON)

Create a JSON file with arrays of thickness and energy values, then pass it to the executable:

```bash
./SEE_in_vacuum ../scan.json
```

Example `scan.json`:
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

Run the ROOT macro to draw and save canvases (with annotations) into the ROOT file:

```bash
root -l /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/draw_histo.C
```

When run on scan outputs, canvases are also exported to `plots/` as `.root` and `.pdf`,
mirroring the `results/` subfolder structure.

To plot a specific file from a parametric scan:
```bash
root -l '/Users/luca/Documents/software/GEANT4/SEE_in_vacuum/draw_histo.C("SEE_in_vacuum_thick20nm_energy1MeV.root")'
```

To process all scan outputs at once:
```bash
bash /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/scripts/run_draw_all.sh /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/results/scan_thick5-10-15-20-25nm_energy1MeV_events100000
```

To create summary ROOT files with overlapping histograms and legends (one per particle):
```bash
bash /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/scripts/run_draw_summary.sh /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/results/scan_thick5-10-15-20-25nm_energy1MeV_events100000
```

To compare models per energy (one canvas per energy value):
```bash
bash /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/scripts/run_draw_summary_models.sh \
  /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/results/scan_*_modelPAI \
  /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/results/scan_*_modelLivermore \
  /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/results/scan_*_modelPenelope \
  /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/results/summary_models.root
```

Model comparison canvases are also exported to `plots/` alongside the summary root file.

Notes:
- The energy deposition plot uses log Y scale by default.
- Annotations read `RunMeta` from the ROOT file.
- Plot titles/legends include the selected EM model (PAI/Livermore/Penelope).
- EM low-energy cutoffs are set to 0.1 eV (see `src/PhysicsList.cc`).

## Customization

You can modify the electron gun parameters in `src/PrimaryGeneratorAction.cc`:
- Energy: Change `SetParticleEnergy()`
- Position: Change `SetParticlePosition()`
- Direction: Change `SetParticleMomentumDirection()`

You can also modify the Al2O3 geometry in `src/DetectorConstruction.cc`:
- Thickness: Change the `thickness` variable
- Diameter: Change the `radius` variable

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
