# Geometry 2D Slice Monte Carlo

This note summarizes the 2D slice Monte Carlo and analytical scans for shell crossings.

## Overview

The 2D slice MC models a thin slice through a 3D slab filled with non-overlapping core-shell spheres. Each sphere appears as a circle whose radius depends on where the slice intersects the sphere. The MC traces vertical rays through the slice and counts shell crossings.

## Configuration

The main configuration file is `config/geometry/geometry_config.json`.

Key fields:
- `reference`: baseline parameters (`L_um`, `a_nm`, `d_nm`, `phi`)
- `scan_ranges`: arrays for `phi`, `d_nm`, `a_nm`, `L_um`
- `montecarlo_2d`: MC settings (`W_um`, `n_rays`, `seed`, `scan_slice_plots`)

## Running

Analytical scans:
```bash
conda run -n geant4 python geometry_analytical.py
```

2D MC scans:
```bash
conda run -n geant4 python geometry_montecarlo_2d.py --scan
```

Regenerate plots from cached results (no MC run):
```bash
conda run -n geant4 python geometry_montecarlo_2d.py --plot-only
```

## Outputs

Scan summary cache:
- `results/geometry/geometry_mc_2d_scan_summary.json`

Scan plots:
- `plots/geometry/mc_2d_scan_phi.pdf|root`
- `plots/geometry/mc_2d_scan_d_nm.pdf|root`
- `plots/geometry/mc_2d_scan_a_nm.pdf|root`
- `plots/geometry/mc_2d_scan_L_um.pdf|root`
- `plots/geometry/crossings_vs_core_radius_mc_vs_ana_2d.pdf|root`

Slice configuration examples:
- `plots/geometry/slice_configs_2d/slice_config_*_zoom.pdf|root`

Crossing distributions per scan value:
- `plots/geometry/distributions_2d/crossings_dist_*.pdf|root`

## Filling Factor Convention

The input `phi` is always the **3D volume filling factor** used to set the sphere density. The achieved filling factor is estimated from the 2D slice area (Delesse principle) and is reported as a **volume** filling factor. Plots annotate this as `Achieved #phi (vol)` and the `phi` scan uses achieved `phi` on the x-axis while listing the target values in the text box.
