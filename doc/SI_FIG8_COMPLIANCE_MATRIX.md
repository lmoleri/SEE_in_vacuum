# Si Fig. 3-14 (Paper Fig. 8) Compliance Matrix

## Scope

This checklist compares:

- thesis requirements for Fig. 3-14 (paper Fig. 8), and
- the current implementation + JSON configurations in this repository.

Primary thesis source used: `/tmp/codex_docs/thesis.pdf`.
This matrix was extended after a full-thesis paragraph sweep (pages 1-349), keeping only implementation points that can affect the Si Fig. 3-14 benchmark.

Status convention:

- `PASS`: requirement is implemented and present in the selected config path.
- `PARTIAL`: capability exists but default/selected path is mixed or not fully aligned.
- `FAIL`: current default path does not satisfy the requirement.

## Matrix

| ID | Thesis requirement (Fig. 3-14) | Expected setup | Evidence in repo | Status | Action |
|---|---|---|---|---|---|
| R1 | Silicon TEY benchmark on flat target, normal incidence | Si target; electron beam; normal incidence | Si-only geometry is supported via zero film + Si substrate: `src/DetectorConstruction.cc:125`, `src/DetectorConstruction.cc:156`; auto/explicit Si scoring: `src/main.cc:1183`, `src/main.cc:1196`; normal-incidence default (`+z`) when no incidence override: `src/main.cc:997`, `src/main.cc:998` | PASS | None |
| R2 | Compare three modes (a/b/c) | (a) all corrections; (b) no weakly-bound initial energy; (c) no initial energy + no surface | Mode parsing and mapping: `src/main.cc:621`, `src/main.cc:634`, `src/main.cc:637`, `src/main.cc:640` | PASS | None |
| R3 | Mode (b): weakly-bound initial energy disabled | `WeaklyBoundInitialEnergy = 0` in MicroElec data | Data override machinery: `src/main.cc:246`, `src/main.cc:332`, activation in scans: `src/main.cc:728` | PASS | None |
| R4 | Mode (c): no surface potential barrier | Surface process disabled | Surface process gating: `src/PhysicsList.cc:471`, `src/PhysicsList.cc:535`; mode (c) forces `microelec_surface_enabled=false`: `src/main.cc:640` | PASS | None |
| R5 | Si material constants (Table 3-7) | Surface barrier `4.05 eV`; weakly-bound initial energy `1.73 eV` | EMLOW Si data file has these values: `/Users/luca/miniforge3/pkgs/geant4-data-emlow-8.6.1-hd8ed1ab_0/share/Geant4/data/EMLOW8.6.1/microelec/Structure/Data_Si.dat` | PASS | Keep data version pinned/documented |
| R6 | MicroElec low-energy transport in Si for TEY | MicroElec elastic + inelastic + surface in Si region | Region/process configuration: `src/PhysicsList.cc:447`, `src/PhysicsList.cc:519`, `src/PhysicsList.cc:524`, `src/PhysicsList.cc:535`, model attachment: `src/PhysicsList.cc:564`, `src/PhysicsList.cc:586` | PASS | None |
| R7 | TEY count for irradiated surface (entrance side) for Fig. 3-14 comparison | Entrance-side first-exit accounting per track; TEY from emitted tracks / primaries | Paper counters: `src/RunAction.cc:382`, `src/RunAction.cc:385`, `src/RunAction.cc:1223`; registration path: `src/SteppingAction.cc:364`, `src/EventAction.cc:249` | PASS | Use `teyEntrancePaper` for Fig. 3-14 overlays |
| R8 | SEY/BEY 50 eV split (used in thesis context) | `<50 eV` vs `>=50 eV` split from entrance side | Implemented in paper counters: `src/RunAction.cc:1233`, `src/RunAction.cc:1239` | PASS | None |
| R9 | Run energy domain comparable to Fig. 3-14 | include ~20-1600 eV | Frozen baseline scans include 20-1600 eV in all three configs: `config/geant4/scan_si_fig8_paper_mode_a_sub5000nm_r100000nm.json`, `config/geant4/scan_si_fig8_paper_mode_b_sub5000nm_r100000nm.json`, `config/geant4/scan_si_fig8_paper_mode_c_sub5000nm_r100000nm.json` | PASS | None |
| R10 | Benchmark config must stay aligned (deexcitation/transport assumptions) | consistent a/b/c configs used for paper plot | Frozen baseline a/b/c configs keep `disable_deexcitation=false`, `deexcitation_ignore_cut=false`, `microelec_regime=\"sey_benchmark\"`, `microelec_fig8_strict_mode=true`, `microelec_lo_phonon_enabled=false`, `microelec_capture_enabled=false` | PASS | Use only the frozen baseline config family for Fig. 3-14 runs. |
| R11 | Avoid adding non-thesis Si processes by default | no Si LO-phonon/capture unless intentionally tested | LO phonon and capture are now opt-in only: no strict-mode auto-enable in `src/main.cc`; process insertion remains gated in `src/PhysicsList.cc:529`; paper configs explicitly pin `microelec_lo_phonon_enabled=false` and `microelec_capture_enabled=false` | PASS | None |
| R12 | Semi-infinite Si substrate approximation | effectively thick substrate | Re-check done for mode (a) with larger geometry/statistics (`r=100 um`, `50k` events/point) across `1/5/10 um`, summarized in `plots/si_fig8_substrate_sensitivity_modea_sub1-5-10um_teyEntrancePaper_r100um_events50k/si_fig8_substrate_sensitivity_modea_summary.txt`: max relative differences `1um vs 5um = 2.99%`, `1um vs 10um = 2.98%`, `5um vs 10um = 1.09%`. Baseline is now pinned to `sub=5 um` for paper-mode configs. | PASS | None |
| R13 | Plotting pipeline should default to paper metric for Fig. 3-14 | default metric = paper TEY | Benchmark script default is `teyEntrancePaper`: `scripts/benchmark_si_microelec_fig8_modes.py:310` | PASS | None |
| R14 | Legacy configs should not silently diverge from paper baseline | avoid disabling deexcitation in paper benchmark | Legacy marker configs were removed during cleanup; only paper-mode a/b/c configs are retained under `config/geant4/` | PASS | None |
| R15 | Geant4 baseline consistency with thesis/paper MicroElec context | same Geant4 generation as benchmark reference, or explicit equivalence proof | Current runs in this repo are Geant4 `11.3.p02`: `results/_logs/scan_si_fig8_paper_mode_a_sub5000nm_r100000_baseline.log:6` | PARTIAL | Keep this as an explicit caveat (latest Geant4 used here; version-gap may exist vs thesis baseline). |
| R16 | Surface model fidelity to Eq. 3-26 implementation details | verify exponential barrier transmission law (`a=0.5e-10 m`), specular reflection, and boundary energy bookkeeping match thesis assumptions | Thesis model description: pages 93-95 (Eq. 3-26 and interface handling). Repo delegates this to upstream `G4MicroElecSurface` (no local implementation): `src/PhysicsList.cc:535-538` | PARTIAL | Keep using upstream process, but add/retain source-level audit notes for the exact Geant4 version used in benchmark claims. |
| R17 | Stopping rule at surface barrier (tracking cut semantics) | electrons tracked until energy below material threshold, then removed (Si threshold/work function context in thesis) | Si data carries `WorkFunction=4.05 eV`: `/Users/luca/miniforge3/pkgs/geant4-data-emlow-8.6.1-hd8ed1ab_0/share/Geant4/data/EMLOW8.6.1/microelec/Structure/Data_Si.dat:3`; surface process active in mode (a): `results/_logs/scan_si_fig8_paper_mode_a_sub5000nm_r100000_baseline.log:676` | PARTIAL | Add a dedicated runtime diagnostic counter of stop reasons (`below barrier` vs other) for Si mode-(a) benchmark runs. |
| R18 | Weakly-bound initial-energy correction details | mode (a)/(b) should follow thesis logic for weakly-bound initial energy treatment in inelastic sampling | Mode toggles are implemented (`microelec_mode` a/b): `src/main.cc:634-640`, with Data_Si override in (b): `src/main.cc:728-745`; thesis low-energy safeguard details described on pages 90-92 | PARTIAL | Audit upstream `G4MicroElecInelasticModel_new` behavior vs thesis text (especially low-energy fallback when transfer exceeds available energy). |
| R19 | Auger/deexcitation contribution in benchmark regime | shell-vacancy deexcitation path must be active and not inadvertently suppressed by cuts | Runtime confirms Auger on and `DeexcitationIgnoreCut=0`: `results/_logs/scan_si_fig8_paper_mode_a_sub5000nm_r100000_baseline.log:139-141`; effective cuts are low (0.1 eV) in Si region: `...baseline.log:640-643` | PARTIAL | Keep current low cuts; add an explicit Auger-yield toggle check (`deexcitation_ignore_cut` true/false) at one probe energy to bound sensitivity. |
| R20 | Low-energy process composition in Si region | benchmark should not be dominated by standard low-energy `msc/eIoni` in Si when MicroElec is targeted | In SiRegion, `eIoni` is suppressed below 10 MeV and `msc` below 100 MeV (`src/PhysicsList.cc:553-584`); strict mode also removes extra standard processes from e- stack (`src/PhysicsList.cc:501-507`) and disables `StepLimiter` injection (`src/PhysicsList.cc:346-363`). Baseline configs pin `microelec_fig8_strict_mode=true`. | PASS | None |
| R21 | Vacuum-side material semantics at the emitting boundary | benchmark boundary should match thesis vacuum/material assumptions without hidden aliasing differences | Frozen baseline a/b/c configs are now explicit `world_material=\"Vacuum\"`; default fallback remains `G4_Galactic` in code (`src/main.cc:596`, `src/DetectorConstruction.cc:26`) | PARTIAL | Keep one side-by-side sanity run (`Vacuum` vs `G4_Galactic` alias) to quantify any residual bias. |

## Recommended "paper-faithful" config family

Use these as baseline for Fig. 3-14 reproduction:

- `config/geant4/scan_si_fig8_paper_mode_a_sub5000nm_r100000nm.json`
- `config/geant4/scan_si_fig8_paper_mode_b_sub5000nm_r100000nm.json`
- `config/geant4/scan_si_fig8_paper_mode_c_sub5000nm_r100000nm.json`

These configs are frozen with explicit:

- `world_material: "Vacuum"`
- `incidence_angle_normal_deg: 0.0`
- `sample_thickness_nm: [0]` and `substrate_thickness_nm: [5000]`
- `microelec_fig8_strict_mode: true`
- `microelec_lo_phonon_enabled: false`, `microelec_capture_enabled: false`

with comparison metric set explicitly to:

- `teyEntrancePaper`

## Runtime compliance probe (strict baseline)

Latest probe run (300 eV, 2000 events/mode) using the frozen baseline configs:

- report: `results/si_fig8_runtime_compliance_probe_300eV_events2000/si_fig8_runtime_compliance_checks.md`
- table: `results/si_fig8_runtime_compliance_probe_300eV_events2000/si_fig8_runtime_compliance_checks.csv`
- launcher script: `scripts/probe_si_fig8_runtime_compliance.py`

Result: all checks passed (`12/12`) for modes `a`, `b`, `c`:

- strict profile active (`microelec_fig8_strict_mode=true`)
- mode mapping enforced (`disable_initial_energy`, `surface`)
- process stack is strict-clean (`StepLimiter`, `eBrem`, `ePairProd`, `CoulombScat` absent)
- MicroElec Si processes present (`e-_G4Dielectrics`, `e-_G4MicroElecElastic`)
- mode-dependent surface process presence correct (`a/b` on, `c` off)

## Full strict baseline rerun (a/b/c)

After freezing configs and passing runtime compliance checks, full scans were rerun for all three modes:

- `results/scan_si_fig8_paper_modea_sub5000nm_r100000nm_particlee-_energy20-1600eV_events10000_refprofile_strict_inc90surf_activeSi_mefig8a_strict`
- `results/scan_si_fig8_paper_modeb_sub5000nm_r100000nm_particlee-_energy20-1600eV_events10000_refprofile_strict_inc90surf_activeSi_mefig8b_strict`
- `results/scan_si_fig8_paper_modec_sub5000nm_r100000nm_particlee-_energy20-1600eV_events10000_refprofile_strict_inc90surf_activeSi_mefig8c_strict`

Updated comparison overlay:

- `plots/si_fig8_mode_benchmark_paper_baseline_teyEntrancePaper_refprofile_strict_sub5um_r100um_events10k/si_microelec_fig8_mode_comparison_abs.pdf`
- `plots/si_fig8_mode_benchmark_paper_baseline_teyEntrancePaper_refprofile_strict_sub5um_r100um_events10k/si_microelec_fig8_mode_comparison_logy.pdf`
- `plots/si_fig8_mode_benchmark_paper_baseline_teyEntrancePaper_refprofile_strict_sub5um_r100um_events10k/si_microelec_fig8_mode_comparison.root`

## Highest-priority fixes

1. Resolve R15/R16/R18 equivalence risk (version/process-logic gap) before declaring Fig. 3-14 reproduction.
2. Close R17/R19/R21 with short probe diagnostics (barrier-stop reason split, Auger sensitivity toggle, `Vacuum` vs `G4_Galactic` boundary check).
