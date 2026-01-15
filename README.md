# GEANT4 SEE in Vacuum Simulation

This GEANT4 simulation models an electron gun shooting electrons at an Al2O3 (aluminum oxide) layer.

## Geometry

- **Target**: Al2O3 layer
  - Thickness: 20 nm
  - Diameter: 100 nm
  - Shape: Cylindrical disk

- **Electron Gun**:
  - Particle: Electron (e-)
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
  "events": 100000,
  "output_dir": "results/scan_thick10-20-50nm_energy0p5-1-2MeV_events100000"
}
```

Output files are created inside `output_dir` for each combination. When `output_dir`
is relative, it is resolved from the project root (not `build/`), e.g.:
`results/scan_thick10-20-50nm_energy0p5-1-2MeV_events100000/SEE_in_vacuum_thick20nm_energy1MeV_events100000.root`

## Outputs

- `SEE_in_vacuum.root` (ROOT file with histograms and canvases)
  - `EdepPrimary`: primary e- energy deposition in Al2O3 (eV)
  - `EdepInteractions`: number of energy-depositing steps in Al2O3 per event
  - `RunMeta`: ntuple with `primaryEnergyMeV` and `sampleThicknessNm`
  - `EdepPrimaryCanvas`, `EdepInteractionsCanvas`: canvases saved with annotations

## Plotting

Run the ROOT macro to draw and save canvases (with annotations) into the ROOT file:

```bash
root -l /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/draw_histo.C
```

To plot a specific file from a parametric scan:
```bash
root -l '/Users/luca/Documents/software/GEANT4/SEE_in_vacuum/draw_histo.C("SEE_in_vacuum_thick20nm_energy1MeV.root")'
```

To process all scan outputs at once:
```bash
bash /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/scripts/run_draw_all.sh /Users/luca/Documents/software/GEANT4/SEE_in_vacuum/results/scan_thick5-10-15-20-25nm_energy1MeV_events100000
```

Notes:
- The energy deposition plot uses log Y scale by default.
- Overflow entries are folded into the last visible bin for `EdepPrimary`.
- Annotations read `RunMeta` from the ROOT file.

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
