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

## Building the Project

1. Create a build directory:
```bash
mkdir build
cd build
```

2. Configure with CMake:
```bash
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
./SEE_in_vacuum run.mac
```

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
/run/beamOn 1000

# Optional: Change electron energy
/gun/energy 5 MeV

# Optional: Change electron position
/gun/position 0 0 -1 um

# Optional: Change electron direction
/gun/direction 0 0 1
```
