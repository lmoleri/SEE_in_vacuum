# TODO

- Validation of MC simulation for electrons on slab (reproduce existing SEE data).
- Full MC simulation of MIP detection including electron multiplication (no backscattering; all electrons absorbed in shell; electric field histogram from Comsol; geometrical mean free path histogram).
- Move all python scripts into a dedicated folder.
- Reorganize Python scripts with proper documentation (clear structure, purpose, inputs/outputs, and usage examples for each script).
- Verify attenuation factor in depth-weighted model: current attenuation over 5 nm appears too small (~10%); check against thickness dependence in data, which suggests saturation beyond about 5 nm.
- Run step-size convergence scan for electron transport in 5 nm Al2O3 (`max_step_nm = 0.4, 0.2, 0.1, 0.05, 0.025`) over the sensitive energy range (`350-900 eV`) with high statistics (`>=100k` events/point).
- Extract and compare convergence QoIs across step sizes: `f_stop`, `f_opp`, `f_ent`, `mode(EdepPrimary)`, `mean(EdepPrimary)`, median/q90 of `PrimaryExitEnergyOpposite`, and final toy-model SEY.
- Define and apply objective convergence acceptance criteria between successive step halvings (`<2%` for class fractions; `<10 eV` for energy metrics), then choose the largest converged step.
- Repeat the converged-step scan with alternative electron transport settings (msc / single-scattering choices) to separate step-size effects from model-form effects.
- Benchmark transmission/backscatter observables against thin-film electron transport data before using Geant4 energy-deposition outputs for SEY fitting.
- Add event-level diagnostics ntuple for primary electrons (`EdepPrimary`, `Eexit`, `exitClass`, `stopStatus`, `endVolume`, `nSteps`, `trackLength`, `maxDepth`, boundary-crossing counters, direction reversals, first/last process).
- Add per-event process-resolved energy budget (`Edep_by_eIoni`, `Edep_by_msc`, `Edep_by_other`) plus `Edep_firstStep`, `Edep_maxStep`, and `depth_firstEdep`.
- Add class-conditioned `EdepPrimary` histograms (`stop`, `entrance-exit`, `opposite-exit`, `lateral-exit`) to separate class mixing from straggling within one class.
- Add explicit boundary-step diagnostics for primary crossings (`Al2O3->World` and `World->Al2O3`: crossing position, kinetic energy, process, step status, step length).
- Add event-level energy-closure check (`E0`, `EdepPrimary`, primary exit energy, escaped-secondary energy if present) and histogram `deltaE = E0 - accounted_energy`.
- Add configurable step-level sampled ntuple for primary steps (`stepIndex`, `preE`, `postE`, `edep`, `stepLen`, `depth`, `process`, `stepStatus`) to inspect transport pathologies.
- Add scattering diagnostics (`cos(theta_z)` per step, cumulative angular deflection, max deflection) and compare by exit class.
- Add path-length and max-depth diagnostics split by exit class to test whether bimodality comes from distinct transport topologies.
