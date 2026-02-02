#!/usr/bin/env python3
"""
Monte Carlo simulation of shell crossings for a particle traversing a slab
filled with core-shell spheres.

Algorithm:
1. Generate random sphere positions in slab (non-overlapping)
2. Trace rays (particles) along z-axis through the slab
3. For each ray-sphere intersection, determine shell crossings
4. Collect statistics and compare with analytical prediction
"""

import json
import os
import array
import numpy as np
from geometry_analytical import expected_crossings, expected_sphere_intersections, crossings_per_sphere


def generate_sphere_positions(L_um, a_nm, d_nm, phi, xy_size_um, seed=None, max_attempts_per_sphere=10000):
    """
    Generate non-overlapping sphere positions in a slab using random sequential addition.
    
    Parameters:
        L_um: slab thickness in um (z direction)
        a_nm: core radius in nm
        d_nm: shell thickness in nm
        phi: target filling factor
        xy_size_um: size of slab in x and y directions (um)
        seed: random seed
        max_attempts_per_sphere: max attempts per sphere before giving up
    
    Returns:
        positions: array of sphere centers (x, y, z) in nm
        actual_phi: achieved filling factor
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Convert to nm
    L_nm = L_um * 1000.0
    xy_size_nm = xy_size_um * 1000.0
    R_nm = a_nm + d_nm  # total sphere radius
    
    # Volume of slab and single sphere
    V_slab = xy_size_nm * xy_size_nm * L_nm
    V_sphere = (4.0 / 3.0) * np.pi * R_nm**3
    
    # Target number of spheres
    N_target = int(phi * V_slab / V_sphere)
    print(f"  Target: {N_target} spheres for phi={phi}")
    
    positions = []
    consecutive_failures = 0
    
    while len(positions) < N_target:
        placed = False
        for _ in range(max_attempts_per_sphere):
            # Random position (sphere center must be at least R from boundaries in z,
            # use periodic boundaries in x, y)
            x = np.random.uniform(0, xy_size_nm)
            y = np.random.uniform(0, xy_size_nm)
            z = np.random.uniform(R_nm, L_nm - R_nm)
            
            # Check overlap with existing spheres (including periodic images in x, y)
            overlaps = False
            for px, py, pz in positions:
                # Check distance considering periodic boundaries in x, y
                dx = min(abs(x - px), xy_size_nm - abs(x - px))
                dy = min(abs(y - py), xy_size_nm - abs(y - py))
                dz = abs(z - pz)
                dist_sq = dx**2 + dy**2 + dz**2
                if dist_sq < (2 * R_nm)**2:  # spheres overlap if centers closer than 2R
                    overlaps = True
                    break
            
            if not overlaps:
                positions.append((x, y, z))
                placed = True
                consecutive_failures = 0
                break
        
        if not placed:
            consecutive_failures += 1
            if consecutive_failures >= 10:
                print(f"  Stopping early: placed {len(positions)}/{N_target} spheres")
                break
        
        if len(positions) % 500 == 0:
            print(f"  Placed {len(positions)}/{N_target} spheres...")
    
    positions = np.array(positions) if positions else np.array([]).reshape(0, 3)
    actual_phi = len(positions) * V_sphere / V_slab
    
    if len(positions) < N_target:
        print(f"  Final: {len(positions)}/{N_target} spheres, phi = {actual_phi:.4f} (target: {phi})")
    else:
        print(f"  Placed all {len(positions)} spheres, phi = {actual_phi:.4f}")
    
    return positions, actual_phi


def trace_ray(x_ray, y_ray, positions, L_nm, a_nm, d_nm, xy_size_nm):
    """
    Trace a ray along z-axis and count shell crossings.
    
    Parameters:
        x_ray, y_ray: ray position in x, y (nm)
        positions: sphere centers array
        L_nm: slab thickness in nm
        a_nm: core radius in nm
        d_nm: shell thickness in nm
        xy_size_nm: slab xy size (for periodic boundaries)
    
    Returns:
        n_crossings: number of shell crossings
        n_spheres: number of spheres intersected
    """
    R_nm = a_nm + d_nm
    n_crossings = 0
    n_spheres = 0
    
    for sx, sy, sz in positions:
        # Distance from ray to sphere center (considering periodic boundaries)
        dx = min(abs(x_ray - sx), xy_size_nm - abs(x_ray - sx))
        dy = min(abs(y_ray - sy), xy_size_nm - abs(y_ray - sy))
        b = np.sqrt(dx**2 + dy**2)  # impact parameter
        
        if b < R_nm:
            # Ray intersects this sphere
            n_spheres += 1
            
            if b < a_nm:
                # Ray goes through core: 2 shell crossings
                n_crossings += 2
            else:
                # Ray grazes shell only: 1 shell crossing
                n_crossings += 1
    
    return n_crossings, n_spheres


def run_montecarlo(config_path=None):
    """
    Run Monte Carlo simulation and compare with analytical prediction.
    """
    import ROOT
    from ROOT import TCanvas, TH1D, TGraph, TGraphErrors, TLatex, TLegend, gStyle
    
    if config_path is None:
        config_path = "config/geometry/geometry_config.json"
    
    with open(config_path, 'r') as f:
        cfg = json.load(f)
    
    ref = cfg["reference"]
    mc_cfg = cfg["montecarlo"]
    plots_dir = cfg.get("plots_dir", "plots/geometry")
    output_dir = cfg.get("output_dir", "results/geometry")
    
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # Parameters (MC can override L for faster validation)
    L_um = mc_cfg.get("L_um_override", ref["L_um"])
    a_nm = ref["a_nm"]
    d_nm = ref["d_nm"]
    phi = ref["phi"]
    xy_size_um = mc_cfg["slab_xy_um"]
    n_rays = mc_cfg["n_rays"]
    seed = mc_cfg["seed"]
    
    L_nm = L_um * 1000.0
    xy_size_nm = xy_size_um * 1000.0
    R_nm = a_nm + d_nm
    
    import sys
    print("Shell Crossing Geometry - Monte Carlo Simulation")
    print("=" * 50)
    print(f"Parameters: L={L_um} um, a={a_nm} nm, d={d_nm} nm, phi={phi}")
    print(f"Slab xy size: {xy_size_um} um, N_rays: {n_rays}")
    print()
    sys.stdout.flush()
    
    # Generate sphere positions
    print("Generating sphere positions...")
    positions, actual_phi = generate_sphere_positions(
        L_um, a_nm, d_nm, phi, xy_size_um, seed=seed
    )
    print(f"Placed {len(positions)} spheres, actual phi = {actual_phi:.4f}")
    print()
    
    # Trace rays
    print(f"Tracing {n_rays} rays...")
    np.random.seed(seed + 1)  # different seed for ray positions
    
    crossings_list = []
    spheres_list = []
    
    for i in range(n_rays):
        x_ray = np.random.uniform(0, xy_size_nm)
        y_ray = np.random.uniform(0, xy_size_nm)
        n_cross, n_sph = trace_ray(x_ray, y_ray, positions, L_nm, a_nm, d_nm, xy_size_nm)
        crossings_list.append(n_cross)
        spheres_list.append(n_sph)
        
        if (i + 1) % 2000 == 0:
            print(f"  {i + 1}/{n_rays} rays traced...")
    
    crossings_arr = np.array(crossings_list)
    spheres_arr = np.array(spheres_list)
    
    # Statistics
    mean_crossings = np.mean(crossings_arr)
    std_crossings = np.std(crossings_arr)
    mean_spheres = np.mean(spheres_arr)
    std_spheres = np.std(spheres_arr)
    
    # Analytical predictions (using actual phi)
    N_analytical = expected_crossings(L_um, a_nm, d_nm, actual_phi)
    N_spheres_analytical = expected_sphere_intersections(L_um, a_nm, d_nm, actual_phi)
    crossings_per_sph = crossings_per_sphere(a_nm, d_nm)
    
    print()
    print("Results:")
    print(f"  Mean shell crossings: {mean_crossings:.2f} +/- {std_crossings:.2f}")
    print(f"  Analytical prediction: {N_analytical:.2f}")
    print(f"  Ratio MC/Analytical: {mean_crossings / N_analytical:.3f}")
    print()
    print(f"  Mean sphere intersections: {mean_spheres:.2f} +/- {std_spheres:.2f}")
    print(f"  Analytical prediction: {N_spheres_analytical:.2f}")
    print(f"  Crossings per sphere (MC): {mean_crossings / mean_spheres:.3f}")
    print(f"  Crossings per sphere (analytical): {crossings_per_sph:.3f}")
    
    gStyle.SetOptStat(110)
    
    # Plot 1: Distribution of shell crossings
    max_cross = int(crossings_arr.max())
    h_cross = TH1D("h_crossings", "Shell crossings distribution;Number of shell crossings;Number of rays",
                   max_cross + 1, -0.5, max_cross + 0.5)
    for c in crossings_list:
        h_cross.Fill(c)
    
    c1 = TCanvas("c_mc_dist", "MC crossings distribution", 800, 600)
    c1.SetGrid()
    h_cross.SetFillColor(38)
    h_cross.SetFillStyle(1001)
    h_cross.Draw()
    
    lat = TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.03)
    lat.DrawLatex(0.55, 0.85, f"MC mean = {mean_crossings:.1f}")
    lat.DrawLatex(0.55, 0.80, f"Analytical = {N_analytical:.1f}")
    lat.DrawLatex(0.55, 0.75, f"Ratio = {mean_crossings / N_analytical:.3f}")
    lat.DrawLatex(0.15, 0.85, f"L={L_um} #mum, #phi={actual_phi:.3f}")
    c1.Update()
    c1.SaveAs(os.path.join(plots_dir, "mc_crossings_distribution.pdf"))
    c1.SaveAs(os.path.join(plots_dir, "mc_crossings_distribution.root"))
    print(f"\nSaved: {plots_dir}/mc_crossings_distribution.pdf")
    
    # Plot 2: Analytical vs MC for varying phi
    print("\nRunning phi scan for comparison plot...")
    sys.stdout.flush()
    phi_scan = cfg["scan_ranges"]["phi"]
    mc_means = []
    mc_stds = []
    analytical_values = []
    actual_phis = []
    
    for phi_val in phi_scan:
        print(f"  phi = {phi_val}...")
        sys.stdout.flush()
        # Generate new sphere configuration
        pos, actual_p = generate_sphere_positions(
            L_um, a_nm, d_nm, phi_val, xy_size_um, seed=seed + int(phi_val * 1000)
        )
        actual_phis.append(actual_p)
        
        # Trace fewer rays for scan (faster)
        n_scan_rays = min(500, n_rays)
        crossings = []
        for _ in range(n_scan_rays):
            x_ray = np.random.uniform(0, xy_size_nm)
            y_ray = np.random.uniform(0, xy_size_nm)
            n_cross, _ = trace_ray(x_ray, y_ray, pos, L_nm, a_nm, d_nm, xy_size_nm)
            crossings.append(n_cross)
        
        mc_means.append(np.mean(crossings))
        mc_stds.append(np.std(crossings) / np.sqrt(n_scan_rays))  # standard error
        analytical_values.append(expected_crossings(L_um, a_nm, d_nm, actual_p))
        print(f"    MC={mc_means[-1]:.1f}, Ana={analytical_values[-1]:.1f}")
        sys.stdout.flush()
    
    # Use actual achieved phi values for x-axis
    c2 = TCanvas("c_comparison", "Analytical vs MC", 800, 600)
    c2.SetGrid()
    
    n_pts = len(phi_scan)
    gr_mc = TGraphErrors(n_pts,
                         array.array('d', actual_phis),
                         array.array('d', mc_means),
                         array.array('d', [0.0] * n_pts),
                         array.array('d', mc_stds))
    gr_mc.SetMarkerStyle(20)
    gr_mc.SetMarkerSize(1.2)
    gr_mc.SetMarkerColor(4)
    gr_mc.SetLineColor(4)
    gr_mc.SetLineWidth(2)
    
    gr_ana = TGraph(n_pts,
                    array.array('d', actual_phis),
                    array.array('d', analytical_values))
    gr_ana.SetMarkerStyle(22)
    gr_ana.SetMarkerSize(1.2)
    gr_ana.SetMarkerColor(2)
    gr_ana.SetLineColor(2)
    gr_ana.SetLineWidth(2)
    gr_ana.SetLineStyle(2)
    
    # Draw
    gr_mc.SetTitle("Analytical vs Monte Carlo;Filling factor (#phi);Shell crossings")
    gr_mc.Draw("APL")
    gr_ana.Draw("PL SAME")
    
    leg = TLegend(0.15, 0.75, 0.45, 0.88)
    leg.SetBorderSize(1)
    leg.SetFillStyle(0)
    leg.AddEntry(gr_mc, "Monte Carlo", "lpe")
    leg.AddEntry(gr_ana, "Analytical", "lp")
    leg.Draw()
    
    lat.DrawLatex(0.55, 0.85, f"L={L_um} #mum, a={a_nm} nm, d={d_nm} nm")
    c2.Update()
    c2.SaveAs(os.path.join(plots_dir, "analytical_vs_mc.pdf"))
    c2.SaveAs(os.path.join(plots_dir, "analytical_vs_mc.root"))
    print(f"Saved: {plots_dir}/analytical_vs_mc.pdf")
    
    # Write summary
    summary_path = os.path.join(output_dir, "geometry_mc_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("Shell Crossing Geometry - Monte Carlo Results\n")
        f.write("=" * 50 + "\n\n")
        f.write("Parameters:\n")
        f.write(f"  L = {L_um} um\n")
        f.write(f"  a = {a_nm} nm\n")
        f.write(f"  d = {d_nm} nm\n")
        f.write(f"  Target phi = {phi}\n")
        f.write(f"  Actual phi = {actual_phi:.4f}\n")
        f.write(f"  Slab xy size = {xy_size_um} um\n")
        f.write(f"  N_rays = {n_rays}\n")
        f.write(f"  N_spheres placed = {len(positions)}\n\n")
        f.write("Results:\n")
        f.write(f"  Mean shell crossings (MC) = {mean_crossings:.2f} +/- {std_crossings:.2f}\n")
        f.write(f"  Analytical prediction = {N_analytical:.2f}\n")
        f.write(f"  Ratio MC/Analytical = {mean_crossings / N_analytical:.4f}\n\n")
        f.write(f"  Mean sphere intersections (MC) = {mean_spheres:.2f} +/- {std_spheres:.2f}\n")
        f.write(f"  Analytical prediction = {N_spheres_analytical:.2f}\n\n")
        f.write(f"  Crossings per sphere (MC) = {mean_crossings / mean_spheres:.4f}\n")
        f.write(f"  Crossings per sphere (analytical) = {crossings_per_sph:.4f}\n")
    
    print(f"\nSummary written to: {summary_path}")
    
    return {
        "mean_crossings": mean_crossings,
        "std_crossings": std_crossings,
        "analytical": N_analytical,
        "mean_spheres": mean_spheres,
        "actual_phi": actual_phi,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Monte Carlo shell crossing simulation")
    parser.add_argument("--config", "-c", default="config/geometry/geometry_config.json",
                        help="Path to config JSON")
    args = parser.parse_args()
    
    run_montecarlo(args.config)


if __name__ == "__main__":
    main()
