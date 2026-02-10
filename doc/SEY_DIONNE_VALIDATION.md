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
- In the paper, `α` is listed without explicit units. The validation script now defaults to
  interpreting `α` in **1/nm**. Use `--alpha-unit angstrom` (or `--alpha-scale 10`)
  if you want to apply an Å → nm conversion.

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

`E_a` is not explicitly listed in the paper; the script uses a default of `E_a = 10 eV`
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
     --Ea 10
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
  --em-model Penelope
```

### Step-level depth weighting (electrons only)

When `sey_alpha_inv_nm` is set in the JSON scan config, Geant4 will also produce
`EdepPrimaryWeighted`, which applies the per-step escape weighting:

$$
\Delta E_{\mathrm{weighted}} = \sum_i \Delta E_i e^{-\alpha z_i}
$$

Here `z_i` is measured from the **entrance side** (the side the primary impinges on).

To use this weighted histogram in the toy MC, pass `--histogram EdepPrimaryWeighted`
and `--depth-model weighted` (or `--weighted-edep`):

```bash
conda run -n geant4 python calculate_muon_sey.py \
  --input-dir results/scan_dionne_validation_5nm_particlee-_energy100-1000eV_step100eV_events10000_modelPenelope \
  --histogram EdepPrimaryWeighted \
  --weighted-edep
```

Note: the SEY stored in `RunMeta` is computed from **emitted electrons** (all e- leaving the Al2O3 layer,
including primaries), matching the usual definition of total yield.

## Effective production cuts and step limits (Al2O3 region)

We dumped the **effective production cuts** actually used in the Al2O3 region. With the current
defaults, the range cut is **1 mm**, which corresponds to **hundreds of keV** for electrons in Al2O3.
That means low-energy secondaries are not produced in the Al2O3 region at all.

If you are using Geant4 **only for energy deposition** and compute secondaries in a custom MC,
this is actually desirable: **large production cuts suppress Geant4 secondaries and leave energy
as local deposition**.

Example output (Al2O3 region, 5 nm scan):

```

### Max step size in Al2O3

To get meaningful **depth-resolved** energy deposition in a 5 nm layer, you may need
to limit the **maximum step length** inside Al2O3. This is controlled by `max_step_nm`
in the scan JSON:

```json
{
  "max_step_nm": 0.1
}
```

This **does not create secondaries**. It only forces Geant4 to break the primary
track into smaller steps so you can resolve the depth profile.

### Suppressing Geant4‑generated secondaries (electrons)

If you are using Geant4 **only for energy deposition**, you typically want **no secondary
electrons produced by Geant4**, because you handle multiplication in the custom MC.

Two knobs matter:

1. **Production cuts (range cuts)** — large cuts suppress low‑energy secondaries.
2. **Atomic de‑excitation (Fluo/Auger)** — can still generate secondaries even when
   cuts are large if `DeexcitationIgnoreCut` is true.

Recommended settings in the scan JSON:

```json
{
  "disable_deexcitation": true,
  "deexcitation_ignore_cut": false
}
```

These settings disable Fluo/Auger for all EM models and ensure any remaining de‑excitation
respects cuts.

## Baseline transport diagnostics (5 nm Al2O3, no substrate)

To sanity-check Geant4 transport before fitting SEY, we run a **baseline scan**
for 5 nm Al2O3 with no substrate and inspect several diagnostics.

**Baseline result folders:**
- Penelope: `results/scan_dionne_validation_5nm_baseline_penelope_sub0nm_r100nm`
- Livermore: `results/scan_dionne_validation_5nm_baseline_livermore_sub0nm_r100nm`

**Key plots (all in `plots/MC_electrons_on_shell_dionne-model/diagnostics_edep_depth/`):**
- `edep_depth_vs_energy_baseline_penelope(.pdf/.png)`  
  `edep_depth_vs_energy_baseline_livermore(.pdf/.png)`  
  Heatmap of **energy deposition vs depth** (primary e- only). Bright near 0 depth
  indicates deposition happens close to the entrance surface.
- `edep_depth_weighted_vs_energy_baseline_*(.pdf/.png)`  
  Same as above, but weighted by escape probability `exp(-alpha z)`.
- `tracklen_depth_vs_energy_*_baseline.png`  
  **Primary track length vs depth** (path-length weighted). This answers
  “where the primary travels,” not just where it deposits energy.
- `mean_depth_tracklen_vs_energy_*_baseline.png`  
  Mean depth of deposition and mean primary track length vs energy,
  with a 5 nm reference line.
- `full_edep_fraction_vs_energy_*_baseline.png`  
  Fraction of events with `EdepPrimary >= 0.98 E0` (a proxy for “full stopping”).
- `end_volume_fraction_vs_energy_*_baseline.png`  
  Fraction of events ending in **Al2O3 / World / OutOfWorld / Other**.

**dE/dx diagnostics (in `plots/MC_electrons_on_shell_dionne-model/diagnostics_dedx/`):**
- `dedx_vs_energy_baseline_penelope.pdf`
- `dedx_vs_energy_baseline_livermore.pdf`

### Notes on depth vs track-length diagnostics

`EdepDepthPrimary` is **energy deposition vs depth** (energy-weighted).
It is normal for this to peak near the entrance surface because low-energy
electrons deposit energy quickly.

`PrimaryTrackLengthDepth` is **track length vs depth** (path-length weighted),
and is the correct diagnostic for “where the primary actually travels” inside
the film.
gamma cut: 1 mm -> ~7.1 keV
e- cut:    1 mm -> ~826 keV
e+ cut:    1 mm -> ~791 keV
proton:    1 mm -> 100 keV
```

Global EM step settings (shown by Geant4):

```
Step function for e+-: (0.2, 0.01 mm)
MscRangeFactor (e-/e+): 0.08
MscLambdaLimit: 1 mm
```

Implication: these defaults are **orders of magnitude larger than the 5 nm layer thickness** and can
push the effective energy-deposition peak to higher energies. To study the ~300 eV peak region, we
will need to tighten production cuts and step limits specifically for the Al2O3 region.
