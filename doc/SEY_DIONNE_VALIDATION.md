# SEY Validation vs Dionne Model (Fig. 9)

This validation compares Geant4 Monte Carlo SEY (emitted electrons per primary) against the
analytic Dionne model used in the attached AIP Advances paper (Fig. 9).

**Scope:** The Dionne analytic model here is intended **only for electron primaries** (`e-`/`e+`).

## Dionne model

The SEY is modeled as:

$$
\delta = \frac{B}{E_a}\left(\frac{A n}{\alpha}\right)^{1/n}(\alpha d)^{1/n-1}\left(1-e^{-\alpha d}\right)
$$

$$
d = \frac{E_p^n}{A^n}
$$


- `E_p` is the primary electron energy (eV).
- `E_a` is the inner SE excitation energy (eV).
- `A`, `B`, `n`, `α` are material parameters (see Table II/III in the paper).
- `d` is the penetration depth (nm) implied by the model (Eq. 5).
- In the paper, `α` is listed without explicit units. The validation script defaults to
  interpreting `α` in **1/Å** and converts it to **1/nm** by multiplying by 10.
  Use `--alpha-unit nm` (or `--alpha-scale 1`) to disable this scaling if needed.

## Escape probability vs depth (toy MC, option 2)

When we do not track the exact depth of each energy-depositing step, the toy MC
approximates the escape probability by assuming an exponential depth distribution:

$$
P_{\mathrm{esc}}(z) = B e^{-\alpha z}
$$

$$
f(z) = \frac{1}{d} e^{-z/d}
$$

$$
\langle P_{\mathrm{esc}} \rangle
  = \int_0^\infty P_{\mathrm{esc}}(z) f(z)\, dz
  = \frac{B}{1+\alpha d}
$$

- `B` is the surface escape probability (z = 0).
- `α` is the attenuation coefficient.
- `d` is the mean penetration depth (nm) from Dionne Eq. 5.

This gives an energy-dependent escape probability through `d(E)`, which naturally
produces a peak and a high-energy falloff in the SEY curve. The alternative (option 1)
is to weight each Geant4 step by `exp(-α z)` using its actual depth, which avoids this
assumption.

## Parameter presets

The validation script ships with presets taken from the paper:

- **Table II (intrinsic):**
  - Al2O3: `B=0.46`, `A=37`, `n=1.61`, `α=0.0075`
  - Si: `B=0.26`, `A=20`, `n=1.38`, `α=0.040`
- **Table III (mixed material / thickness):**
  - 1 nm: `B=0.382`, `A=25`, `n=1.45`, `α=0.030`
  - 3 nm: `B=0.425`, `A=35`, `n=1.57`, `α=0.013`
  - 5 nm: `B=0.450`, `A=36`, `n=1.60`, `α=0.0090`

`E_a` is not explicitly listed in the paper; the script uses a default of `E_a = 1.0 eV`
and allows overrides via CLI.

## Running the validation

1. **Run an electron scan** to produce SEY metadata in the ROOT output:
   - Example config: `config/geant4/scan_dionne_validation_20nm.json`
   - Run:
     ```bash
     conda run -n geant4 ./build/SEE_in_vacuum config/geant4/scan_dionne_validation_20nm.json
     ```

2. **Generate the validation plot**:
   ```bash
   conda run -n geant4 python scripts/validate_sey_dionne.py \
     --results-dir results/scan_thick20nm_particlee-_energy100-1000eV_step100eV_events10000_modelPenelope \
     --material al2o3 \
     --Ea 1.0
   ```

Output is saved in `plots/MC_electrons_on_shell_dionne-model/mc_vs_dionne/`.

## 5 nm validation (Table III)

For the 5 nm Al2O3 case (Table III parameters), you can validate against different EM models
by running the scan with a chosen `em_model` in the config and passing `--em-model` to the validator.

```bash
conda run -n geant4 ./build/SEE_in_vacuum config/geant4/scan_dionne_validation_5nm.json
conda run -n geant4 python scripts/validate_sey_dionne.py \
  --results-dir results/scan_dionne_validation_5nm_particlee-_energy100-1000eV_step100eV_events10000_modelPenelope \
  --material al2o3_5nm \
  --em-model Penelope \
  --fit-fig9 \
  --fig9-ea-on-material
```

### Dual-fit overlay (Fig. 9 + MC peak)

To overlay both analytic curves (Fig. 9 fit and MC-peak fit) on the same plot:

```bash
conda run -n geant4 python scripts/validate_sey_dionne.py \
  --results-dir results/scan_dionne_validation_5nm_particlee-_energy100-1000eV_step100eV_events10000_modelPenelope \
  --material al2o3_5nm \
  --mc-source toy \
  --em-model Penelope \
  --fit-fig9 \
  --fit-mc-peak
```

Note: the SEY stored in `RunMeta` is computed from **emitted electrons** (all e- leaving the Al2O3 layer,
including primaries), matching the usual definition of total yield.
