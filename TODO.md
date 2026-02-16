# TODO

- Validation of MC simulation for electrons on slab (reproduce existing SEE data).
- Evaluate and integrate Geant4-MicroElec extension for low-energy electron transport/energy-deposition studies relevant to SEE and REELS validation.
- Run a thickness scan at fixed primary energy and validate transmission law from `PrimaryExitClass` (`class 2` fraction): fit `T(d) = T_inf + (T0 - T_inf) exp(-d/Lambda_eff)`, check log-linearity (`ln(T)` or `ln(T-T_inf)` vs `d`), and compare `Lambda_eff` across `max_step_nm = 0.1, 0.2, 0.3`.
- Full MC simulation of MIP detection including electron multiplication (no backscattering; all electrons absorbed in shell; electric field histogram from Comsol; geometrical mean free path histogram).
- Reorganize Python scripts with proper documentation (clear structure, purpose, inputs/outputs, and usage examples for each script).
- Verify attenuation factor in depth-weighted model: current attenuation over 5 nm appears too small (~10%); check against thickness dependence in data, which suggests saturation beyond about 5 nm.
- Run step-size convergence scan for electron transport in 5 nm Al2O3 (`max_step_nm = 0.4, 0.2, 0.1, 0.05, 0.025`) over the sensitive energy range (`350-900 eV`) with high statistics (`>=100k` events/point).
- Extract and compare convergence QoIs across step sizes: `f_stop`, `f_opp`, `f_ent`, `mode(EdepPrimary)`, `mean(EdepPrimary)`, median/q90 of `PrimaryExitEnergyOpposite`, and final toy-model SEY.
- Define and apply objective convergence acceptance criteria between successive step halvings (`<2%` for class fractions; `<10 eV` for energy metrics), then choose the largest converged step.
- Repeat the converged-step scan with alternative electron transport settings (msc / single-scattering choices) to separate step-size effects from model-form effects.
- Expose msc transport knobs in JSON config and wire them to `G4EmParameters`: `MscStepLimitType`, `MscRangeFactor`, `MscGeomFactor`, `MscSafetyFactor`, `MscLambdaLimit`, `MscSkin`, `MscThetaLimit`, `MscEnergyLimit`, `TransportationWithMsc`, `SingleScatteringType`, plus `StepFunction`.
- Benchmark transmission/backscatter observables against thin-film electron transport data before using Geant4 energy-deposition outputs for SEY fitting.
- Add event-level energy-closure check (`E0`, `EdepPrimary`, primary exit energy, escaped-secondary energy if present) and histogram `deltaE = E0 - accounted_energy`.
- Extend boundary-step diagnostics from first outward crossing to full crossing history (`Al2O3->World` and `World->Al2O3`: crossing position, kinetic energy, process, step status, step length).
- Extend scattering diagnostics with cumulative/max angular-deflection observables (beyond current per-step `deltaThetaDeg`/`dirZ`) and compare by exit class.
- Complete path-length and max-depth diagnostics as explicit split-by-class summary plots (all classes).
