#!/usr/bin/env python3
"""
2D Slice Monte Carlo for shell crossing geometry.

This simulates a true 2D slice through a 3D slab filled with non-overlapping
core-shell spheres. Circles in the slice have varying radii depending on
where the slice cuts through each sphere.

Key features:
- Circle outer radius: r_outer = sqrt(R² - h²) where h is slice height
- Circle core radius: r_core = sqrt(a² - h²) if |h| < a, else no core visible
- Volume fraction phi is the 3D parameter (determines sphere density)
- Non-overlapping circles (guaranteed by 3D sphere non-overlap)
"""

import json
import os
import sys
import array
import time
import numpy as np
from geometry_analytical import expected_crossings


def format_param(value):
    """
    Format a parameter value for filenames (e.g., 0.5 -> 0p5).
    """
    if abs(value - round(value)) < 1e-9:
        out = f"{int(round(value))}"
    else:
        out = f"{value:.6g}"
    return out.replace(".", "p")


def format_seconds(value):
    if value < 1.0:
        return f"{value * 1000.0:.1f} ms"
    if value < 60.0:
        return f"{value:.2f} s"
    return f"{value / 60.0:.2f} min"


def select_plot_indices(n_values, mode):
    """
    Select which indices to plot in a scan.
    mode can be:
      - "all"
      - "none"
      - "first"
      - "last"
      - "first_last"
      - int (number of evenly spaced plots)
    """
    if n_values <= 0:
        return set()
    if mode is None:
        mode = "all"
    if isinstance(mode, int):
        count = max(0, mode)
        if count <= 0:
            return set()
        if count == 1:
            return {0}
        if count >= n_values:
            return set(range(n_values))
        indices = set()
        for i in range(count):
            idx = int(round(i * (n_values - 1) / (count - 1)))
            indices.add(idx)
        return indices
    if isinstance(mode, str):
        key = mode.strip().lower()
        if key == "all":
            return set(range(n_values))
        if key == "none":
            return set()
        if key == "first":
            return {0}
        if key == "last":
            return {n_values - 1}
        if key == "first_last":
            return {0, n_values - 1} if n_values > 1 else {0}
    return set(range(n_values))


def generate_circles_2d(L_um, W_um, a_nm, d_nm, phi, seed=None, max_attempts_per_circle=5000):
    """
    Generate circles representing a 2D slice through 3D sphere distribution.
    
    Parameters:
        L_um: slab thickness in um (z-direction, vertical)
        W_um: slice width in um (x-direction, horizontal)
        a_nm: core radius in nm
        d_nm: shell thickness in nm
        phi: 3D volume fraction
        seed: random seed
        max_attempts_per_circle: max placement attempts per circle
    
    Returns:
        circles: list of (x, z, r_outer, r_core) tuples in nm
        actual_phi: achieved volume fraction (approximate)
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Convert to nm
    L_nm = L_um * 1000.0
    W_nm = W_um * 1000.0
    R_nm = a_nm + d_nm  # total sphere radius
    
    # 3D sphere density from volume fraction
    # n_3D = phi / V_sphere = 3*phi / (4*pi*R³)
    V_sphere = (4.0 / 3.0) * np.pi * R_nm**3
    n_3D = phi / V_sphere
    
    # Expected circles in slice: spheres with centers within R of slice plane
    # n_2D = n_3D * 2R (density per unit area of slice)
    n_2D = n_3D * 2.0 * R_nm
    
    # Expected number of circles in W x L area
    N_target = int(n_2D * W_nm * L_nm)
    print(f"  2D slice: W={W_um} um, L={L_um} um")
    print(f"  Target: {N_target} circles for phi={phi}")
    sys.stdout.flush()
    
    circles = []
    consecutive_failures = 0

    # Use a simple cell list to accelerate overlap checks.
    # Max circle radius is R, so max overlap distance is 2R.
    # With cell size = 2R, any overlapping circles must be in the same
    # or neighboring cells (3x3 neighborhood).
    cell_size = 2.0 * R_nm
    nx = max(1, int(np.ceil(W_nm / cell_size)))
    nz = max(1, int(np.ceil(L_nm / cell_size)))
    grid = [[[] for _ in range(nz)] for _ in range(nx)]

    def cell_index(x, z):
        ix = int(x / cell_size)
        iz = int(z / cell_size)
        ix = min(max(ix, 0), nx - 1)
        iz = min(max(iz, 0), nz - 1)
        return ix, iz
    
    # 2D slice geometry (slice plane at y = 0):
    #
    #        y
    #        ^          sphere center at y = h
    #        |                   *
    #        |                .-' '-.
    #  y=0 --+-------------- /-------\  <-- slice plane
    #        |               \-------/
    #        |                '-.__.-'
    #
    # For each sphere center offset h (uniform in [-R, R]):
    #   r_outer = sqrt(R^2 - h^2)
    #   r_core  = sqrt(a^2 - h^2) if |h| < a, else no core visible.
    #
    while len(circles) < N_target:
        placed = False
        for _ in range(max_attempts_per_circle):
            # Random position in slice
            x = np.random.uniform(0, W_nm)
            z = np.random.uniform(R_nm, L_nm - R_nm)  # keep circles inside
            
            # Random slice height through sphere (uniform in [-R, R])
            h = np.random.uniform(-R_nm, R_nm)
            
            # Outer circle radius from slice
            r_outer = np.sqrt(R_nm**2 - h**2)
            
            # Core circle radius (only if slice passes through core)
            if abs(h) < a_nm:
                r_core = np.sqrt(a_nm**2 - h**2)
            else:
                r_core = 0.0  # no core visible
            
            # Check overlap with nearby circles only (cell list)
            overlaps = False
            ix, iz = cell_index(x, z)
            for dix in (-1, 0, 1):
                nix = (ix + dix) % nx  # periodic in x
                for diz in (-1, 0, 1):
                    niz = iz + diz
                    if niz < 0 or niz >= nz:
                        continue
                    for cx, cz, cr_out, cr_core in grid[nix][niz]:
                        # Use periodic boundary in x (horizontal)
                        dx = abs(x - cx)
                        if dx > 0.5 * W_nm:
                            dx = W_nm - dx
                        dz = abs(z - cz)
                        if (dx * dx + dz * dz) < (r_outer + cr_out) ** 2:
                            overlaps = True
                            break
                    if overlaps:
                        break
                if overlaps:
                    break
            
            if not overlaps:
                circles.append((x, z, r_outer, r_core))
                grid[ix][iz].append((x, z, r_outer, r_core))
                placed = True
                consecutive_failures = 0
                break
        
        if not placed:
            consecutive_failures += 1
            if consecutive_failures >= 20:
                print(f"  Stopping early: placed {len(circles)}/{N_target} circles")
                sys.stdout.flush()
                break
        
        if len(circles) % 500 == 0 and len(circles) > 0:
            print(f"  Placed {len(circles)}/{N_target} circles...")
            sys.stdout.flush()
    
    # Estimate actual phi from placed circles
    # This is approximate since circles have varying radii
    total_circle_area = sum(np.pi * r_out**2 for _, _, r_out, _ in circles)
    slice_area = W_nm * L_nm
    area_fraction = total_circle_area / slice_area
    # Convert back to approximate volume fraction
    # area_fraction ≈ (average cross-section) * n_3D * 2R
    # For uniform h: <pi*r²> = <pi*(R²-h²)> = pi*R² * (2/3) 
    # So area_fraction ≈ (2/3)*pi*R² * n_3D * 2R = (4/3)*pi*R³ * n_3D = phi
    actual_phi = area_fraction  # approximate
    
    if len(circles) >= N_target:
        print(f"  Placed all {len(circles)} circles")
    print(f"  Area fraction in slice: {area_fraction:.4f}")
    sys.stdout.flush()
    
    return circles, actual_phi


def trace_ray_2d(x_ray, circles, W_nm):
    """
    Trace a vertical ray (along z, through slab) in the 2D slice and count shell crossings.
    
    The ray travels vertically at horizontal position x_ray, passing through
    the entire slab thickness. This matches the 3D geometry where particles
    traverse the slab along its thickness.
    
    Parameters:
        x_ray: horizontal position of ray in nm
        circles: list of (x, z, r_outer, r_core) tuples
        W_nm: slice width for periodic boundary
    
    Returns:
        n_crossings: number of shell crossings
        n_circles: number of circles intersected
    """
    n_crossings = 0
    n_circles = 0
    
    for cx, cz, r_outer, r_core in circles:
        # Horizontal distance from ray to circle center (with periodic boundary)
        dx = abs(x_ray - cx)
        dx = min(dx, W_nm - dx)  # periodic boundary in x
        
        if dx < r_outer:
            # Ray intersects this circle
            n_circles += 1
            
            if r_core > 0 and dx < r_core:
                # Ray goes through core: 2 shell crossings
                n_crossings += 2
            else:
                # Ray only in shell region (or no core visible): 1 crossing
                n_crossings += 1
    
    return n_crossings, n_circles


def trace_rays_2d(circles, W_nm, n_rays, seed, verbose=False):
    """
    Trace multiple rays and return crossings and circle-hit arrays.
    """
    np.random.seed(seed)
    crossings = np.zeros(n_rays, dtype=int)
    circles_hit = np.zeros(n_rays, dtype=int)

    for i in range(n_rays):
        x_ray = np.random.uniform(0, W_nm)
        n_cross, n_circ = trace_ray_2d(x_ray, circles, W_nm)
        crossings[i] = n_cross
        circles_hit[i] = n_circ
        if verbose and (i + 1) % 2000 == 0:
            print(f"  {i + 1}/{n_rays} rays traced...")
            sys.stdout.flush()

    return crossings, circles_hit


def plot_slice_configuration(circles, W_nm, L_nm, out_base, title, params_text,
                             max_circles=None, seed=123):
    """
    Plot an example 2D slice configuration (outer + inner circles).
    """
    import ROOT
    from ROOT import TCanvas, TH2D, TEllipse, TLatex, TPaveText, gStyle, gPad

    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)

    n_total = len(circles)
    n_show = n_total
    indices = np.arange(n_total)
    if max_circles is not None and max_circles > 0 and n_total > max_circles:
        rng = np.random.default_rng(seed)
        indices = rng.choice(indices, size=max_circles, replace=False)
        n_show = len(indices)

    W_um = W_nm * 1.0e-3
    L_um = L_nm * 1.0e-3
    # Choose canvas size so the drawable area (after margins) matches the data aspect.
    aspect = L_um / W_um if W_um > 0 else 1.0
    left_margin = 0.12
    right_margin = 0.06
    bottom_margin = 0.12
    top_margin = 0.28

    width = 900
    drawable_w = width * (1.0 - left_margin - right_margin)
    drawable_h = drawable_w * aspect
    height = int(drawable_h / (1.0 - top_margin - bottom_margin))

    if height > 1200:
        height = 1200
        drawable_h = height * (1.0 - top_margin - bottom_margin)
        drawable_w = drawable_h / aspect if aspect > 0 else drawable_h
        width = int(drawable_w / (1.0 - left_margin - right_margin))
    if height < 500:
        height = 500
        drawable_h = height * (1.0 - top_margin - bottom_margin)
        drawable_w = drawable_h / aspect if aspect > 0 else drawable_h
        width = int(drawable_w / (1.0 - left_margin - right_margin))
    width = max(500, min(1200, width))

    c = TCanvas(f"c_slice_{format_param(seed)}", title, width, height)
    c.SetGrid()

    gPad.SetLeftMargin(left_margin)
    gPad.SetBottomMargin(bottom_margin)
    gPad.SetTopMargin(top_margin)
    gPad.SetRightMargin(right_margin)

    frame = TH2D(f"frame_slice_{format_param(seed)}", "", 10, 0.0, W_um, 10, 0.0, L_um)
    frame.SetXTitle("x (#mum)")
    frame.SetYTitle("z (#mum)")
    frame.Draw("AXIS")
    c.Update()
    gPad.SetFixedAspectRatio(1)
    gPad.Modified()
    gPad.Update()

    ellipses = []
    for idx in indices:
        x, z, r_out, r_core = circles[int(idx)]
        x_um = x * 1.0e-3
        z_um = z * 1.0e-3
        r_out_um = r_out * 1.0e-3
        outer = TEllipse(x_um, z_um, r_out_um)
        outer.SetLineColor(ROOT.kBlue + 2)
        outer.SetLineWidth(1)
        outer.SetFillStyle(0)
        outer.Draw()
        ellipses.append(outer)

        if r_core > 0.0:
            r_core_um = r_core * 1.0e-3
            inner = TEllipse(x_um, z_um, r_core_um)
            inner.SetLineColor(ROOT.kRed + 1)
            inner.SetLineWidth(1)
            inner.SetFillStyle(0)
            inner.Draw()
            ellipses.append(inner)

    text = TPaveText(0.15, 0.80, 0.85, 0.98, "NDC")
    text.SetFillColor(0)
    text.SetBorderSize(1)
    text.SetTextAlign(12)
    text.SetTextSize(0.028)
    text.SetTextFont(42)
    text.AddText(title)
    text.AddText(params_text)
    if n_show < n_total:
        text.AddText(f"Showing {n_show} of {n_total} circles")
    else:
        text.AddText(f"Showing {n_show} circles")
    text.Draw()

    c.Update()
    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".root")


def plot_scan_summary(scan_name, values, mc_mean, mc_sem, ana_target, ana_achieved,
                      plots_dir, annotation_lines):
    import ROOT
    from ROOT import TCanvas, TGraphErrors, TGraph, TLatex, TLegend, gStyle

    gStyle.SetOptStat(0)

    x_label_map = {
        "phi": "Filling factor #phi",
        "d_nm": "Shell thickness d (nm)",
        "a_nm": "Core radius a (nm)",
        "L_um": "Slab thickness L (#mum)",
    }
    x_label = x_label_map.get(scan_name, scan_name)

    n = len(values)
    if n == 0:
        return

    gr_mc = TGraphErrors(n)
    gr_ana_target = TGraph(n)
    gr_ana_ach = TGraph(n)

    for i, x in enumerate(values):
        gr_mc.SetPoint(i, float(x), float(mc_mean[i]))
        gr_mc.SetPointError(i, 0.0, float(mc_sem[i]))
        gr_ana_target.SetPoint(i, float(x), float(ana_target[i]))
        gr_ana_ach.SetPoint(i, float(x), float(ana_achieved[i]))

    c = TCanvas(f"c_scan_{scan_name}", f"2D MC scan: {scan_name}", 900, 650)
    c.SetGrid()

    gr_mc.SetMarkerStyle(20)
    gr_mc.SetMarkerSize(1.0)
    gr_mc.SetLineWidth(2)
    gr_mc.SetTitle(f"2D MC scan vs analytical;{x_label};Mean shell crossings")
    gr_mc.Draw("AP")

    gr_ana_target.SetLineColor(ROOT.kRed + 1)
    gr_ana_target.SetLineStyle(2)
    gr_ana_target.SetLineWidth(2)
    gr_ana_target.Draw("L SAME")

    gr_ana_ach.SetLineColor(ROOT.kBlue + 2)
    gr_ana_ach.SetLineStyle(3)
    gr_ana_ach.SetLineWidth(2)
    gr_ana_ach.Draw("L SAME")

    leg = TLegend(0.60, 0.75, 0.88, 0.90)
    leg.SetBorderSize(0)
    leg.AddEntry(gr_mc, "2D MC mean (SEM)", "lp")
    leg.AddEntry(gr_ana_target, "Analytical (target #phi)", "l")
    leg.AddEntry(gr_ana_ach, "Analytical (achieved #phi)", "l")
    leg.Draw()

    lat = TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.03)
    y = 0.70
    for line in annotation_lines:
        lat.DrawLatex(0.15, y, line)
        y -= 0.04

    out_base = os.path.join(plots_dir, f"mc_2d_scan_{scan_name}")
    c.Update()
    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".root")


def plot_core_radius_mc_vs_ana(values, mc_mean, mc_sem, ana_achieved, plots_dir, annotation_lines):
    import ROOT
    from ROOT import TCanvas, TGraphErrors, TGraph, TLatex, TLegend, gStyle

    gStyle.SetOptStat(0)
    n = len(values)
    if n == 0:
        return

    gr_mc = TGraphErrors(n)
    gr_ana = TGraph(n)
    for i, x in enumerate(values):
        gr_mc.SetPoint(i, float(x), float(mc_mean[i]))
        gr_mc.SetPointError(i, 0.0, float(mc_sem[i]))
        gr_ana.SetPoint(i, float(x), float(ana_achieved[i]))

    c = TCanvas("c_core_radius_mc_vs_ana_2d", "Core radius MC vs analytical (2D)", 800, 600)
    c.SetGrid()
    gr_mc.SetTitle("Core radius scan (2D MC vs analytical);Core radius a (nm);Shell crossings")
    gr_mc.SetMarkerStyle(20)
    gr_mc.SetMarkerSize(1.1)
    gr_mc.SetLineWidth(2)
    gr_mc.Draw("AP")

    gr_ana.SetLineColor(ROOT.kRed + 1)
    gr_ana.SetLineStyle(2)
    gr_ana.SetLineWidth(2)
    gr_ana.Draw("L SAME")

    leg = TLegend(0.15, 0.75, 0.45, 0.88)
    leg.SetBorderSize(1)
    leg.SetFillStyle(0)
    leg.AddEntry(gr_mc, "2D MC mean (SEM)", "lpe")
    leg.AddEntry(gr_ana, "Analytical (achieved #phi)", "l")
    leg.Draw()

    lat = TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.03)
    y = 0.70
    for line in annotation_lines:
        lat.DrawLatex(0.55, y, line)
        y -= 0.04

    out_base = os.path.join(plots_dir, "crossings_vs_core_radius_mc_vs_ana_2d")
    c.Update()
    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".root")

def run_montecarlo_2d_scan(cfg, ref, mc_cfg, plots_dir, output_dir, quick=False):
    scan_ranges = cfg.get("scan_ranges", {})
    if not scan_ranges:
        print("No scan_ranges found in config; skipping scan.")
        return None

    W_um = mc_cfg.get("W_um", 10.0)
    n_rays = mc_cfg.get("n_rays", 2000)
    seed_base = mc_cfg.get("seed", 42)
    plot_max_circles = mc_cfg.get("plot_max_circles", None)
    if plot_max_circles is not None and plot_max_circles <= 0:
        plot_max_circles = None
    scan_plot_mode = mc_cfg.get("scan_slice_plots", "first_last")
    quick_max_values = mc_cfg.get("quick_scan_max_values", 2)
    quick_n_rays = mc_cfg.get("quick_scan_n_rays", max(200, int(n_rays / 4)))
    quick_plot_max_circles = mc_cfg.get("quick_scan_plot_max_circles", 400)
    slice_dir = os.path.join(plots_dir, "slice_configs_2d")
    os.makedirs(slice_dir, exist_ok=True)

    if quick:
        print("[Quick scan] limiting values and rays for faster debug timing.")
        n_rays = quick_n_rays
        if plot_max_circles is None:
            plot_max_circles = quick_plot_max_circles

    summary = {
        "meta": {
            "W_um": W_um,
            "n_rays": n_rays,
            "seed_base": seed_base,
            "plot_max_circles": plot_max_circles,
            "scan_slice_plots": scan_plot_mode,
            "quick_scan": quick,
            "quick_scan_max_values": quick_max_values if quick else None,
        },
        "scans": {},
    }

    scan_keys = ["phi", "d_nm", "a_nm", "L_um"]
    scan_index = 0
    t_scan_start = time.perf_counter()

    for scan_name in scan_keys:
        values = scan_ranges.get(scan_name, [])
        if not values:
            continue
        if quick and len(values) > quick_max_values:
            values = values[:quick_max_values]

        plot_indices = select_plot_indices(len(values), scan_plot_mode)

        print(f"\n[Scan] {scan_name} ({len(values)} values)")
        sys.stdout.flush()

        results = {
            "values": [],
            "mc_mean": [],
            "mc_std": [],
            "mc_sem": [],
            "ana_target": [],
            "ana_achieved": [],
            "area_fraction": [],
            "n_circles": [],
            "mean_circles": [],
        }

        for idx, val in enumerate(values):
            # Base parameters
            L_um = mc_cfg.get("L_um", ref["L_um"])
            a_nm = ref["a_nm"]
            d_nm = ref["d_nm"]
            phi = ref["phi"]

            if scan_name == "phi":
                phi = float(val)
            elif scan_name == "d_nm":
                d_nm = float(val)
            elif scan_name == "a_nm":
                a_nm = float(val)
            elif scan_name == "L_um":
                L_um = float(val)

            seed_case = seed_base + scan_index * 1000 + idx * 7
            print(f"  {scan_name}={val}: L={L_um} um, W={W_um} um, a={a_nm} nm, d={d_nm} nm, phi={phi}")
            sys.stdout.flush()

            t0 = time.perf_counter()
            circles, area_fraction = generate_circles_2d(
                L_um, W_um, a_nm, d_nm, phi, seed=seed_case
            )
            t1 = time.perf_counter()
            crossings_arr, circles_hit_arr = trace_rays_2d(
                circles, W_um * 1000.0, n_rays, seed_case + 1
            )
            t2 = time.perf_counter()

            mean_cross = float(np.mean(crossings_arr))
            std_cross = float(np.std(crossings_arr))
            sem_cross = std_cross / np.sqrt(n_rays) if n_rays > 0 else 0.0
            mean_circles = float(np.mean(circles_hit_arr)) if len(circles_hit_arr) else 0.0

            ana_target = expected_crossings(L_um, a_nm, d_nm, phi)
            ana_achieved = expected_crossings(L_um, a_nm, d_nm, area_fraction)

            results["values"].append(float(val))
            results["mc_mean"].append(mean_cross)
            results["mc_std"].append(std_cross)
            results["mc_sem"].append(sem_cross)
            results["ana_target"].append(float(ana_target))
            results["ana_achieved"].append(float(ana_achieved))
            results["area_fraction"].append(float(area_fraction))
            results["n_circles"].append(int(len(circles)))
            results["mean_circles"].append(mean_circles)

            t3 = None
            if idx in plot_indices:
                # Plot example slice configuration for this parameter set
                tag = f"{scan_name}_{format_param(val)}"
                out_base = os.path.join(slice_dir, f"slice_config_{tag}")
                title = f"2D slice: {scan_name} = {val}"
                params_text = f"L={L_um} #mum, W={W_um} #mum, a={a_nm} nm, d={d_nm} nm, #phi={phi}"
                plot_slice_configuration(
                    circles,
                    W_um * 1000.0,
                    L_um * 1000.0,
                    out_base,
                    title,
                    params_text,
                    max_circles=plot_max_circles,
                    seed=seed_case + 2,
                )
                t3 = time.perf_counter()
            else:
                t3 = time.perf_counter()
            print(
                "    Timings:"
                f" circles={format_seconds(t1 - t0)},"
                f" rays={format_seconds(t2 - t1)},"
                f" plot={format_seconds(t3 - t2)},"
                f" total={format_seconds(t3 - t0)}"
            )
            sys.stdout.flush()

        summary["scans"][scan_name] = results

        annotation = [
            f"W={W_um} #mum, N_rays={n_rays}",
            f"a={ref['a_nm']} nm, d={ref['d_nm']} nm, #phi={ref['phi']}",
        ]
        if scan_name == "a_nm":
            annotation[1] = f"d={ref['d_nm']} nm, #phi={ref['phi']}"
        elif scan_name == "d_nm":
            annotation[1] = f"a={ref['a_nm']} nm, #phi={ref['phi']}"
        elif scan_name == "phi":
            annotation[1] = f"a={ref['a_nm']} nm, d={ref['d_nm']} nm"
        elif scan_name == "L_um":
            annotation[1] = f"a={ref['a_nm']} nm, d={ref['d_nm']} nm, #phi={ref['phi']}"

        plot_scan_summary(
            scan_name,
            results["values"],
            results["mc_mean"],
            results["mc_sem"],
            results["ana_target"],
            results["ana_achieved"],
            plots_dir,
            annotation,
        )

        if scan_name == "a_nm":
            core_annotation = [
                f"L={mc_cfg.get('L_um', ref['L_um'])} #mum, W={W_um} #mum",
                f"d={ref['d_nm']} nm, #phi={ref['phi']}",
            ]
            plot_core_radius_mc_vs_ana(
                results["values"],
                results["mc_mean"],
                results["mc_sem"],
                results["ana_achieved"],
                plots_dir,
                core_annotation,
            )

        scan_index += 1

    t_scan_end = time.perf_counter()
    print(f"\n[Scan] Total time: {format_seconds(t_scan_end - t_scan_start)}")
    sys.stdout.flush()

    summary_path = os.path.join(output_dir, "geometry_mc_2d_scan_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nScan summary written to: {summary_path}")
    sys.stdout.flush()
    return summary

def run_montecarlo_2d(config_path=None):
    """
    Run 2D slice Monte Carlo simulation.
    """
    import ROOT
    from ROOT import TCanvas, TH1D, TGraph, TGraphErrors, TLatex, TLegend, gStyle
    
    if config_path is None:
        config_path = "config/geometry/geometry_config.json"
    
    with open(config_path, 'r') as f:
        cfg = json.load(f)
    
    ref = cfg["reference"]
    mc_cfg = cfg.get("montecarlo_2d", cfg.get("montecarlo", {}))
    plots_dir = cfg.get("plots_dir", "plots/geometry")
    output_dir = cfg.get("output_dir", "results/geometry")
    
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # Parameters
    L_um = mc_cfg.get("L_um", ref["L_um"])
    W_um = mc_cfg.get("W_um", 10.0)
    a_nm = ref["a_nm"]
    d_nm = ref["d_nm"]
    phi = ref["phi"]
    n_rays = mc_cfg.get("n_rays", 10000)
    seed = mc_cfg.get("seed", 42)
    plot_slice = mc_cfg.get("plot_slice", True)
    plot_max_circles = mc_cfg.get("plot_max_circles", None)
    if plot_max_circles is not None and plot_max_circles <= 0:
        plot_max_circles = None
    
    L_nm = L_um * 1000.0
    W_nm = W_um * 1000.0
    R_nm = a_nm + d_nm
    
    print("Shell Crossing Geometry - 2D Slice Monte Carlo")
    print("=" * 55)
    print(f"Parameters: L={L_um} um, W={W_um} um, a={a_nm} nm, d={d_nm} nm, phi={phi}")
    print(f"N_rays: {n_rays}")
    print()
    sys.stdout.flush()
    
    # Generate circles
    print("Generating 2D slice circles...")
    sys.stdout.flush()
    circles, area_fraction = generate_circles_2d(L_um, W_um, a_nm, d_nm, phi, seed=seed)
    print(f"Generated {len(circles)} circles")
    print()
    sys.stdout.flush()

    if plot_slice:
        slice_dir = os.path.join(plots_dir, "slice_configs_2d")
        os.makedirs(slice_dir, exist_ok=True)
        out_base = os.path.join(slice_dir, "slice_config_reference")
        title = "2D slice: reference configuration"
        params_text = f"L={L_um} #mum, W={W_um} #mum, a={a_nm} nm, d={d_nm} nm, #phi={phi}"
        plot_slice_configuration(
            circles,
            W_nm,
            L_nm,
            out_base,
            title,
            params_text,
            max_circles=plot_max_circles,
            seed=seed + 2,
        )
    
    # Trace rays (vertical rays through slab, at random x positions)
    print(f"Tracing {n_rays} vertical rays...")
    sys.stdout.flush()
    crossings_arr, circles_hit_arr = trace_rays_2d(
        circles, W_nm, n_rays, seed + 1, verbose=True
    )
    
    # Statistics
    mean_crossings = np.mean(crossings_arr)
    std_crossings = np.std(crossings_arr)
    mean_circles = np.mean(circles_hit_arr)
    
    # Analytical prediction (3D formula with target phi)
    N_analytical_target = expected_crossings(L_um, a_nm, d_nm, phi)
    # Analytical prediction with achieved area fraction (≈ effective phi)
    N_analytical_achieved = expected_crossings(L_um, a_nm, d_nm, area_fraction)
    
    print()
    print("Results:")
    print(f"  Mean shell crossings (2D MC): {mean_crossings:.2f} +/- {std_crossings:.2f}")
    print(f"  Analytical (target phi={phi}): {N_analytical_target:.2f}")
    print(f"  Analytical (achieved phi={area_fraction:.3f}): {N_analytical_achieved:.2f}")
    print(f"  Ratio 2D_MC / Ana(achieved): {mean_crossings / N_analytical_achieved:.3f}")
    print()
    print(f"  Mean circles hit per ray: {mean_circles:.2f}")
    if mean_circles > 0:
        print(f"  Crossings per circle: {mean_crossings / mean_circles:.3f}")
    sys.stdout.flush()
    
    gStyle.SetOptStat(110)
    
    # Plot: Distribution of shell crossings
    max_cross = int(crossings_arr.max()) if crossings_arr.max() > 0 else 10
    h_cross = TH1D("h_crossings_2d", "Shell crossings distribution (2D slice);Number of shell crossings;Number of rays",
                   max_cross + 1, -0.5, max_cross + 0.5)
    for c in crossings_arr:
        h_cross.Fill(c)
    
    c1 = TCanvas("c_mc_dist_2d", "2D MC crossings distribution", 800, 600)
    c1.SetGrid()
    h_cross.SetFillColor(38)
    h_cross.SetFillStyle(1001)
    h_cross.Draw()
    
    lat = TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.03)
    lat.DrawLatex(0.55, 0.85, f"2D MC mean = {mean_crossings:.1f}")
    lat.DrawLatex(0.55, 0.80, f"Ana (achieved #phi) = {N_analytical_achieved:.1f}")
    lat.DrawLatex(0.55, 0.75, f"Ratio = {mean_crossings / N_analytical_achieved:.3f}")
    lat.DrawLatex(0.15, 0.85, f"L={L_um} #mum, W={W_um} #mum")
    lat.DrawLatex(0.15, 0.80, f"#phi={phi}, {len(circles)} circles")
    c1.Update()
    c1.SaveAs(os.path.join(plots_dir, "mc_crossings_distribution_2d.pdf"))
    c1.SaveAs(os.path.join(plots_dir, "mc_crossings_distribution_2d.root"))
    print(f"\nSaved: {plots_dir}/mc_crossings_distribution_2d.pdf")
    sys.stdout.flush()
    
    # Write summary
    summary_path = os.path.join(output_dir, "geometry_mc_2d_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("Shell Crossing Geometry - 2D Slice Monte Carlo Results\n")
        f.write("=" * 55 + "\n\n")
        f.write("Parameters:\n")
        f.write(f"  L = {L_um} um (slab thickness)\n")
        f.write(f"  W = {W_um} um (slice width)\n")
        f.write(f"  a = {a_nm} nm (core radius)\n")
        f.write(f"  d = {d_nm} nm (shell thickness)\n")
        f.write(f"  R = {R_nm} nm (total sphere radius)\n")
        f.write(f"  phi = {phi} (3D volume fraction)\n")
        f.write(f"  N_rays = {n_rays}\n")
        f.write(f"  N_circles = {len(circles)}\n")
        f.write(f"  Area fraction in slice = {area_fraction:.4f}\n\n")
        f.write("Results:\n")
        f.write(f"  Mean shell crossings (2D MC) = {mean_crossings:.2f} +/- {std_crossings:.2f}\n")
        f.write(f"  Analytical (target phi={phi}) = {N_analytical_target:.2f}\n")
        f.write(f"  Analytical (achieved phi={area_fraction:.4f}) = {N_analytical_achieved:.2f}\n")
        f.write(f"  Ratio 2D_MC / Ana(achieved) = {mean_crossings / N_analytical_achieved:.4f}\n\n")
        f.write(f"  Mean circles hit per ray = {mean_circles:.2f}\n")
        if mean_circles > 0:
            f.write(f"  Crossings per circle = {mean_crossings / mean_circles:.4f}\n")
    
    print(f"Summary written to: {summary_path}")
    sys.stdout.flush()
    
    return {
        "mean_crossings": mean_crossings,
        "std_crossings": std_crossings,
        "analytical_target": N_analytical_target,
        "analytical_achieved": N_analytical_achieved,
        "mean_circles": mean_circles,
        "n_circles": len(circles),
        "area_fraction": area_fraction,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="2D Slice Monte Carlo for shell crossings")
    parser.add_argument("--config", "-c", default="config/geometry/geometry_config.json",
                        help="Path to config JSON")
    parser.add_argument("--scan", action="store_true",
                        help="Run parameter scans and summary plots")
    parser.add_argument("--scan-quick", action="store_true",
                        help="Run a reduced scan with fewer values/rays (debug timing)")
    args = parser.parse_args()
    
    if args.scan:
        with open(args.config, "r") as f:
            cfg = json.load(f)
        ref = cfg["reference"]
        mc_cfg = cfg.get("montecarlo_2d", cfg.get("montecarlo", {}))
        plots_dir = cfg.get("plots_dir", "plots/geometry")
        output_dir = cfg.get("output_dir", "results/geometry")
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        run_montecarlo_2d_scan(cfg, ref, mc_cfg, plots_dir, output_dir, quick=args.scan_quick)
    else:
        run_montecarlo_2d(args.config)


if __name__ == "__main__":
    main()
