#!/usr/bin/env python
import argparse
import math
import os
from dataclasses import dataclass


@dataclass
class DionneParams:
    B: float
    A: float
    n: float
    alpha: float
    Ea: float
    label: str


PRESETS = {
    "al2o3": DionneParams(B=0.46, A=37.0, n=1.61, alpha=0.0075, Ea=1.0, label="Al2O3 (Table II)"),
    "si": DionneParams(B=0.26, A=20.0, n=1.38, alpha=0.040, Ea=1.0, label="Si (Table II)"),
    "al2o3_1nm": DionneParams(B=0.382, A=25.0, n=1.45, alpha=0.030, Ea=1.0, label="Al2O3 1 nm (Table III)"),
    "al2o3_3nm": DionneParams(B=0.425, A=35.0, n=1.57, alpha=0.013, Ea=1.0, label="Al2O3 3 nm (Table III)"),
    "al2o3_5nm": DionneParams(B=0.450, A=36.0, n=1.60, alpha=0.0090, Ea=1.0, label="Al2O3 5 nm (Table III)"),
}


def dionne_delta(Ep_eV: float, params: DionneParams) -> float:
    if Ep_eV <= 0:
        return 0.0
    d = (Ep_eV ** params.n) / (params.A ** params.n)
    term1 = (params.B / params.Ea)
    term2 = ((params.A * params.n) / params.alpha) ** (1.0 / params.n)
    term3 = (params.alpha * d) ** (1.0 / params.n - 1.0)
    term4 = 1.0 - math.exp(-params.alpha * d)
    return term1 * term2 * term3 * term4


def main():
    parser = argparse.ArgumentParser(description="Validate MC SEY against Dionne analytic model (Fig. 9 paper).")
    parser.add_argument("--results-dir", required=True, help="Results directory with ROOT files from scan.")
    parser.add_argument("--material", default="al2o3",
                        help="Preset name: " + ", ".join(PRESETS.keys()))
    parser.add_argument("--B", type=float, default=None, help="Override B parameter")
    parser.add_argument("--A", type=float, default=None, help="Override A parameter")
    parser.add_argument("--n", type=float, default=None, help="Override n parameter")
    parser.add_argument("--alpha", type=float, default=None, help="Override alpha parameter (1/nm)")
    parser.add_argument("--Ea", type=float, default=None, help="Override Ea parameter (eV)")
    parser.add_argument("--output-dir", default="plots/validation", help="Output directory for plots")
    parser.add_argument("--title", default="SEY validation vs Dionne model", help="Plot title")
    args = parser.parse_args()

    if args.material not in PRESETS:
        raise SystemExit(f"Unknown material preset '{args.material}'. Valid: {', '.join(PRESETS.keys())}")

    params = PRESETS[args.material]
    params = DionneParams(
        B=params.B if args.B is None else args.B,
        A=params.A if args.A is None else args.A,
        n=params.n if args.n is None else args.n,
        alpha=params.alpha if args.alpha is None else args.alpha,
        Ea=params.Ea if args.Ea is None else args.Ea,
        label=params.label,
    )

    import ROOT
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    # Collect MC points
    energies_eV = []
    sey_values = []
    sey_errors = []
    thickness_nm = None
    em_model = None

    root_files = [
        os.path.join(args.results_dir, f)
        for f in os.listdir(args.results_dir)
        if f.endswith(".root") and not f.startswith("summary")
    ]
    if not root_files:
        raise SystemExit(f"No ROOT files found in {args.results_dir}")

    for path in sorted(root_files):
        f = ROOT.TFile.Open(path)
        if not f or not f.IsOpen():
            continue
        meta = f.Get("RunMeta")
        if not meta:
            f.Close()
            continue
        meta.GetEntry(0)
        try:
            energy_mev = float(meta.primaryEnergyMeV)
            energy_eV = energy_mev * 1.0e6
            sey = float(meta.sey)
            n_primary = int(meta.primaryElectrons)
            try:
                n_emitted = int(meta.emittedElectrons)
            except Exception:
                n_emitted = int(meta.secondaryElectrons)
            if n_primary > 0:
                err = math.sqrt(max(n_emitted, 0)) / n_primary
            else:
                err = 0.0
            energies_eV.append(energy_eV)
            sey_values.append(sey)
            sey_errors.append(err)
            if thickness_nm is None:
                thickness_nm = float(meta.sampleThicknessNm)
            if em_model is None:
                em_model = str(meta.emModel)
            if str(meta.primaryParticle) not in ("e-", "e+"):
                raise SystemExit(
                    f"RunMeta primaryParticle={meta.primaryParticle} is not an electron. "
                    "Dionne validation is defined only for electron primaries."
                )
        except Exception:
            pass
        f.Close()

    if not energies_eV:
        raise SystemExit("No RunMeta entries with SEY found. Re-run simulations with updated RunAction.")

    # Sort by energy
    order = sorted(range(len(energies_eV)), key=lambda i: energies_eV[i])
    energies_eV = [energies_eV[i] for i in order]
    sey_values = [sey_values[i] for i in order]
    sey_errors = [sey_errors[i] for i in order]

    # Analytic curve
    min_e = min(energies_eV)
    max_e = max(energies_eV)
    n_curve = 400
    step = (max_e - min_e) / (n_curve - 1)
    curve_x = []
    curve_y = []
    for i in range(n_curve):
        e = min_e + i * step
        curve_x.append(e)
        curve_y.append(dionne_delta(e, params))

    # Plot
    c = ROOT.TCanvas("c_sey_dionne", "SEY vs Dionne", 1000, 700)
    frame = ROOT.TH2D("frame", "", 10, min_e * 0.95, max_e * 1.05, 10, 0.0, max(max(sey_values), max(curve_y)) * 1.25)
    frame.SetXTitle("Primary electron energy (eV)")
    frame.SetYTitle("Secondary electron yield (SEY)")
    frame.Draw("AXIS")

    gr_mc = ROOT.TGraphErrors(len(energies_eV))
    for i, (x, y, e) in enumerate(zip(energies_eV, sey_values, sey_errors)):
        gr_mc.SetPoint(i, x, y)
        gr_mc.SetPointError(i, 0.0, e)
    gr_mc.SetMarkerStyle(20)
    gr_mc.SetMarkerSize(1.4)
    gr_mc.SetLineWidth(2)
    gr_mc.SetLineColor(ROOT.kBlue + 2)
    gr_mc.SetMarkerColor(ROOT.kBlue + 2)
    gr_mc.Draw("P")

    gr_ana = ROOT.TGraph(len(curve_x))
    for i, (x, y) in enumerate(zip(curve_x, curve_y)):
        gr_ana.SetPoint(i, x, y)
    gr_ana.SetLineWidth(3)
    gr_ana.SetLineColor(ROOT.kRed + 1)
    gr_ana.Draw("L")

    leg = ROOT.TLegend(0.58, 0.72, 0.92, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(gr_mc, "MC SEY", "p")
    leg.AddEntry(gr_ana, "Dionne model", "l")
    leg.Draw()

    text = ROOT.TPaveText(0.12, 0.72, 0.54, 0.88, "NDC")
    text.SetFillColor(0)
    text.SetBorderSize(1)
    text.SetTextFont(42)
    text.SetTextSize(0.028)
    if thickness_nm is not None:
        text.AddText(f"Thickness: {thickness_nm:.0f} nm")
    if em_model:
        text.AddText(f"EM model: {em_model}")
    text.AddText(params.label)
    text.AddText(f"B={params.B:.3g}, A={params.A:.3g}, n={params.n:.3g}, #alpha={params.alpha:.3g}")
    text.AddText(f"E_a={params.Ea:.3g} eV")
    text.Draw()

    formula = ROOT.TLatex()
    formula.SetNDC(True)
    formula.SetTextFont(42)
    formula.SetTextSize(0.028)
    formula.DrawLatex(
        0.12,
        0.66,
        "#delta = (B/E_{a}) (A n/#alpha)^{1/n} (#alpha d)^{1/n-1} (1-e^{-#alpha d})",
    )

    os.makedirs(args.output_dir, exist_ok=True)
    thick_tag = f"{int(round(thickness_nm))}nm" if thickness_nm is not None else "unknown"
    out_base = os.path.join(args.output_dir, f"sey_dionne_vs_mc_{args.material}_{thick_tag}")
    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".root")


if __name__ == "__main__":
    main()
