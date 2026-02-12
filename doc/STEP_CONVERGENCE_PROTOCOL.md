# Step-Size Convergence Protocol for Electron Transport

## Purpose

This protocol defines how to choose `max_step_nm` for low-energy electron transport in thin Al2O3 layers (5 nm), using objective convergence criteria instead of manual tuning.

The immediate motivation is the observed sensitivity of `EdepPrimary` bimodality and exit fractions to the step limit.

## Scope

- Geometry: Al2O3 thickness `5 nm`, substrate `0 nm`, radius `100 nm`
- Primary: `e-` normal incidence
- Energy range for convergence: `350-900 eV` (sensitive region)
- Physics mode: energy-deposition workflow (no Geant4 secondaries in model output)

## Required settings

Use these settings in all convergence runs:

- `disable_deexcitation: true`
- `deexcitation_ignore_cut: false`
- fixed EM model per campaign (run Penelope and Livermore separately)
- fixed random seed policy (or high statistics)

## Step grid and statistics

Run the same scan for:

- `max_step_nm = 0.4, 0.2, 0.1, 0.05, 0.025`
- energies every `25 eV` in `350-900 eV`
- statistics: `>= 100000` events per energy point

## Quantities of interest (QoIs)

Extract the following per energy:

- `f_stop`, `f_ent`, `f_opp`
- `mean(EdepPrimary)`, `mode(EdepPrimary)`
- opposite-exit residual-energy quantiles: `q50`, `q90` of `PrimaryExitEnergyOpposite`
- final toy-model SEY (same post-processing settings for all step values)

Notes for `PrimaryExitClass`:

- class values are `1=entrance`, `2=opposite`, `3=lateral`
- with histogram range `[0,4]` and 4 bins, these appear in ROOT bins `2,3,4`

## Pairwise convergence test

Compare each step pair `(s, s/2)` for all energies and QoIs.

For fraction-like metrics use relative variation:

$$
\Delta_{\mathrm{rel}}(x; s, s/2)=
\frac{|x(s)-x(s/2)|}{\max(|x(s/2)|,\epsilon)}, \quad \epsilon=10^{-6}
$$

For energy-like metrics use absolute variation:

$$
\Delta_{\mathrm{abs}}(x; s, s/2)=|x(s)-x(s/2)|
$$

Acceptance thresholds:

- fractions (`f_stop`, `f_ent`, `f_opp`): `\Delta_rel < 2%`
- energy metrics (`mean`, `mode`, quantiles): `\Delta_abs < 10 eV`
- toy SEY: `\Delta_rel < 2%`

## Step selection rule

Choose the **largest** step `s*` that passes thresholds against `s*/2` for:

- all QoIs
- all energies in `350-900 eV`

If no step passes globally:

- report non-convergence
- use smallest tested step as temporary fallback
- flag model-form uncertainty

## Robustness after convergence

After choosing `s*`, repeat the scan with alternative electron-transport options (msc/single-scattering choices) to separate:

- numerical step-size effects
- physics-model effects

If QoIs still move significantly, the uncertainty is model-form dominated.

## Benchmark requirement before SEY fitting

Before using Geant4 `EdepPrimary` for final SEY validation:

- benchmark transmission/backscatter observables against thin-film electron transport data
- verify that the simulated threshold behavior vs thickness and energy is realistic

## Recommended run workflow

1. Generate one JSON config per step value with identical settings except `max_step_nm` and `output_dir`.
2. Run all scans.
3. Extract QoIs to a single table with columns:
   `step_nm, E0_eV, f_stop, f_ent, f_opp, mean_edep, mode_edep, q50_eexit_opp, q90_eexit_opp, sey`.
4. Compute pairwise deltas and pass/fail flags.
5. Select `s*` with the rule above and record it in docs/config.

## Output artifacts to keep

- Raw ROOT scan outputs for each step
- Consolidated CSV table of QoIs
- CSV of pairwise deltas and pass/fail
- One summary plot per QoI vs energy with one curve per step
- One final note documenting chosen `s*` and residual uncertainty
