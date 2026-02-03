#!/usr/bin/env python3
"""
Analytical calculation of shell crossings for a particle traversing a slab
filled with core-shell spheres.

Physical model:
- Slab of thickness L filled with spheres at volume fraction phi
- Each sphere has core radius a and shell thickness d (total radius R = a + d)
- A particle traveling through the slab intersects some spheres
- For each intersection:
  - If impact parameter b < a: particle crosses shell twice (through core)
  - If a < b < R: particle grazes shell once (shell-only)

Expected number of shell crossings:
  N_crossings = (3 * L * phi / (4 * R)) * (1 + (a/R)^2)
"""

import json
import os
import numpy as np


def expected_crossings(L_um, a_nm, d_nm, phi):
    """
    Calculate expected number of shell crossings.
    
    Parameters:
        L_um: slab thickness in micrometers
        a_nm: core radius in nanometers
        d_nm: shell thickness in nanometers
        phi: filling factor (volume fraction of spheres)
    
    Returns:
        Expected number of shell crossings
    """
    # Convert to consistent units (nm)
    L_nm = L_um * 1000.0
    R_nm = a_nm + d_nm  # total sphere radius
    
    # Expected number of sphere intersections
    N_spheres = (3.0 * L_nm * phi) / (4.0 * R_nm)
    
    # Expected crossings per sphere
    # P(through core) = (a/R)^2, gives 2 crossings
    # P(graze shell) = 1 - (a/R)^2, gives 1 crossing
    # E[crossings] = 2*(a/R)^2 + 1*(1 - (a/R)^2) = 1 + (a/R)^2
    ratio = a_nm / R_nm
    crossings_per_sphere = 1.0 + ratio**2
    
    return N_spheres * crossings_per_sphere


def expected_sphere_intersections(L_um, a_nm, d_nm, phi):
    """Calculate expected number of sphere intersections."""
    L_nm = L_um * 1000.0
    R_nm = a_nm + d_nm
    return (3.0 * L_nm * phi) / (4.0 * R_nm)


def crossings_per_sphere(a_nm, d_nm):
    """Calculate expected shell crossings per sphere intersection."""
    R_nm = a_nm + d_nm
    ratio = a_nm / R_nm
    return 1.0 + ratio**2


def analytical_formula_lines():
    """Compact analytical formula lines for plot annotations."""
    return [
        "#LTN_{cross}#GT = #frac{3 L #phi}{4 R} (1 + (#frac{a}{R})^{2})",
    ]


def load_config(config_path):
    """Load configuration from JSON file."""
    with open(config_path, 'r') as f:
        return json.load(f)


def run_analytical_scans(config_path=None):
    """
    Run analytical parameter scans and generate plots using ROOT.
    """
    import ROOT
    from ROOT import TCanvas, TGraph, TLatex, TPaveText, gStyle
    
    if config_path is None:
        config_path = "config/geometry/geometry_config.json"
    
    cfg = load_config(config_path)
    ref = cfg["reference"]
    scans = cfg["scan_ranges"]
    plots_dir = cfg.get("plots_dir", "plots/geometry")
    output_dir = cfg.get("output_dir", "results/geometry")
    
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # Reference values
    L_ref = ref["L_um"]
    a_ref = ref["a_nm"]
    d_ref = ref["d_nm"]
    phi_ref = ref["phi"]
    
    # Calculate reference value
    N_ref = expected_crossings(L_ref, a_ref, d_ref, phi_ref)
    print(f"Reference parameters: L={L_ref} um, a={a_ref} nm, d={d_ref} nm, phi={phi_ref}")
    print(f"Expected crossings at reference: {N_ref:.1f}")
    print(f"Expected sphere intersections: {expected_sphere_intersections(L_ref, a_ref, d_ref, phi_ref):.1f}")
    print(f"Crossings per sphere: {crossings_per_sphere(a_ref, d_ref):.3f}")
    print()
    
    gStyle.SetOptStat(0)
    
    import array
    
    # 1. N_crossings vs phi (filling factor)
    phi_values = scans["phi"]
    N_vs_phi = [expected_crossings(L_ref, a_ref, d_ref, p) for p in phi_values]
    
    c1 = TCanvas("c_phi", "Crossings vs phi", 800, 600)
    c1.SetGrid()
    c1.SetRightMargin(0.32)
    gr_phi = TGraph(len(phi_values), 
                    array.array('d', phi_values), 
                    array.array('d', N_vs_phi))
    gr_phi.SetTitle("Shell crossings vs filling factor;Filling factor (#phi);Expected shell crossings")
    gr_phi.SetMarkerStyle(20)
    gr_phi.SetMarkerSize(1.2)
    gr_phi.SetMarkerColor(4)
    gr_phi.SetLineColor(4)
    gr_phi.SetLineWidth(2)
    gr_phi.Draw("APL")
    lat = TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.03)
    text1 = TPaveText(0.70, 0.64, 0.98, 0.88, "NDC")
    text1.SetFillColor(0)
    text1.SetBorderSize(1)
    text1.SetTextAlign(12)
    text1.SetTextSize(0.027)
    text1.AddText(f"L={L_ref} #mum")
    text1.AddText(f"a={a_ref} nm")
    text1.AddText(f"d={d_ref} nm")
    for line in analytical_formula_lines():
        text1.AddText(line)
    text1.Draw()
    c1.Update()
    c1.SaveAs(os.path.join(plots_dir, "crossings_vs_phi.pdf"))
    c1.SaveAs(os.path.join(plots_dir, "crossings_vs_phi.root"))
    print(f"Saved: {plots_dir}/crossings_vs_phi.pdf")
    
    # 2. N_crossings vs d (shell thickness)
    d_values = scans["d_nm"]
    N_vs_d = [expected_crossings(L_ref, a_ref, d, phi_ref) for d in d_values]
    
    c2 = TCanvas("c_d", "Crossings vs shell thickness", 800, 600)
    c2.SetGrid()
    c2.SetRightMargin(0.32)
    gr_d = TGraph(len(d_values),
                  array.array('d', d_values),
                  array.array('d', N_vs_d))
    gr_d.SetTitle("Shell crossings vs shell thickness;Shell thickness d (nm);Expected shell crossings")
    gr_d.SetMarkerStyle(20)
    gr_d.SetMarkerSize(1.2)
    gr_d.SetMarkerColor(3)
    gr_d.SetLineColor(3)
    gr_d.SetLineWidth(2)
    gr_d.Draw("APL")
    text2 = TPaveText(0.70, 0.64, 0.98, 0.88, "NDC")
    text2.SetFillColor(0)
    text2.SetBorderSize(1)
    text2.SetTextAlign(12)
    text2.SetTextSize(0.027)
    text2.AddText(f"L={L_ref} #mum")
    text2.AddText(f"a={a_ref} nm")
    text2.AddText(f"#phi={phi_ref}")
    for line in analytical_formula_lines():
        text2.AddText(line)
    text2.Draw()
    c2.Update()
    c2.SaveAs(os.path.join(plots_dir, "crossings_vs_shell_thickness.pdf"))
    c2.SaveAs(os.path.join(plots_dir, "crossings_vs_shell_thickness.root"))
    print(f"Saved: {plots_dir}/crossings_vs_shell_thickness.pdf")
    
    # 3. N_crossings vs a (core radius)
    a_values = scans["a_nm"]
    N_vs_a = [expected_crossings(L_ref, a, d_ref, phi_ref) for a in a_values]
    
    c3 = TCanvas("c_a", "Crossings vs core radius", 800, 600)
    c3.SetGrid()
    c3.SetRightMargin(0.32)
    gr_a = TGraph(len(a_values),
                  array.array('d', a_values),
                  array.array('d', N_vs_a))
    gr_a.SetTitle("Shell crossings vs core radius;Core radius a (nm);Expected shell crossings")
    gr_a.SetMarkerStyle(20)
    gr_a.SetMarkerSize(1.2)
    gr_a.SetMarkerColor(6)
    gr_a.SetLineColor(6)
    gr_a.SetLineWidth(2)
    gr_a.Draw("APL")
    text3 = TPaveText(0.70, 0.64, 0.98, 0.88, "NDC")
    text3.SetFillColor(0)
    text3.SetBorderSize(1)
    text3.SetTextAlign(12)
    text3.SetTextSize(0.027)
    text3.AddText(f"L={L_ref} #mum")
    text3.AddText(f"d={d_ref} nm")
    text3.AddText(f"#phi={phi_ref}")
    for line in analytical_formula_lines():
        text3.AddText(line)
    text3.Draw()
    c3.Update()
    c3.SaveAs(os.path.join(plots_dir, "crossings_vs_core_radius.pdf"))
    c3.SaveAs(os.path.join(plots_dir, "crossings_vs_core_radius.root"))
    print(f"Saved: {plots_dir}/crossings_vs_core_radius.pdf")
    
    # 4. N_crossings vs L (slab thickness)
    L_values = scans["L_um"]
    N_vs_L = [expected_crossings(L, a_ref, d_ref, phi_ref) for L in L_values]
    
    c4 = TCanvas("c_L", "Crossings vs slab thickness", 800, 600)
    c4.SetGrid()
    c4.SetRightMargin(0.32)
    gr_L = TGraph(len(L_values),
                  array.array('d', L_values),
                  array.array('d', N_vs_L))
    gr_L.SetTitle("Shell crossings vs slab thickness;Slab thickness L (#mum);Expected shell crossings")
    gr_L.SetMarkerStyle(20)
    gr_L.SetMarkerSize(1.2)
    gr_L.SetMarkerColor(7)
    gr_L.SetLineColor(7)
    gr_L.SetLineWidth(2)
    gr_L.Draw("APL")
    text = TPaveText(0.70, 0.64, 0.98, 0.88, "NDC")
    text.SetFillColor(0)
    text.SetBorderSize(1)
    text.SetTextAlign(12)
    text.SetTextSize(0.027)
    text.AddText(f"a={a_ref} nm")
    text.AddText(f"d={d_ref} nm")
    text.AddText(f"#phi={phi_ref}")
    for line in analytical_formula_lines():
        text.AddText(line)
    text.Draw()
    c4.Update()
    c4.SaveAs(os.path.join(plots_dir, "crossings_vs_slab_thickness.pdf"))
    c4.SaveAs(os.path.join(plots_dir, "crossings_vs_slab_thickness.root"))
    print(f"Saved: {plots_dir}/crossings_vs_slab_thickness.pdf")
    
    # Write summary
    summary_path = os.path.join(output_dir, "geometry_analytical_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("Shell Crossing Geometry - Analytical Results\n")
        f.write("=" * 50 + "\n\n")
        f.write("Reference parameters:\n")
        f.write(f"  L = {L_ref} um (slab thickness)\n")
        f.write(f"  a = {a_ref} nm (core radius)\n")
        f.write(f"  d = {d_ref} nm (shell thickness)\n")
        f.write(f"  phi = {phi_ref} (filling factor)\n\n")
        f.write(f"At reference:\n")
        f.write(f"  Total sphere radius R = {a_ref + d_ref} nm\n")
        f.write(f"  Expected sphere intersections = {expected_sphere_intersections(L_ref, a_ref, d_ref, phi_ref):.1f}\n")
        f.write(f"  Crossings per sphere = {crossings_per_sphere(a_ref, d_ref):.3f}\n")
        f.write(f"  Expected shell crossings = {N_ref:.1f}\n\n")
        f.write("Formula:\n")
        f.write("  N_crossings = (3 * L * phi / (4 * R)) * (1 + (a/R)^2)\n\n")
        f.write("Parameter scans:\n")
        f.write(f"  phi: {scans['phi']}\n")
        f.write(f"  d_nm: {scans['d_nm']}\n")
        f.write(f"  a_nm: {scans['a_nm']}\n")
        f.write(f"  L_um: {scans['L_um']}\n")
    
    print(f"\nSummary written to: {summary_path}")
    
    return {
        "reference": {"L": L_ref, "a": a_ref, "d": d_ref, "phi": phi_ref, "N": N_ref},
        "scans": {
            "phi": (phi_values, N_vs_phi),
            "d": (d_values, N_vs_d),
            "a": (a_values, N_vs_a),
            "L": (L_values, N_vs_L),
        }
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Analytical shell crossing calculation")
    parser.add_argument("--config", "-c", default="config/geometry/geometry_config.json",
                        help="Path to config JSON")
    args = parser.parse_args()
    
    run_analytical_scans(args.config)


if __name__ == "__main__":
    main()
