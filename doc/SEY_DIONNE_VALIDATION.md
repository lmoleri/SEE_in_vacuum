# SEY Validation vs Dionne Model (Fig. 9)

This validation compares Geant4 Monte Carlo SEY (emitted electrons per primary) against the
analytic Dionne model used in the attached AIP Advances paper (Fig. 9).

**Scope:** The Dionne analytic model here is intended **only for electron primaries** (`e-`/`e+`).

## Dionne model

The SEY is modeled as:

```
δ = (B / E_a) * (A n / α)^(1/n) * (α d)^(1/n - 1) * (1 - e^(-α d))
d = E_p^n / A^n
```

- `E_p` is the primary electron energy (eV).
- `E_a` is the inner SE excitation energy (eV).
- `A`, `B`, `n`, `α` are material parameters (see Table II/III in the paper).
- `d` is the penetration depth (nm) implied by the model.

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

Output is saved in `plots/validation/`.

Note: the SEY stored in `RunMeta` is computed from **emitted electrons** (all e- leaving the Al2O3 layer,
including primaries), matching the usual definition of total yield.
