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


def _as_str(value) -> str:
    try:
        if value.__class__.__name__ == "LowLevelView":
            return bytes(value).decode("utf-8", errors="ignore").rstrip("\x00")
    except Exception:
        pass
    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8", errors="ignore").rstrip("\x00")
    try:
        return str(value)
    except Exception:
        return ""


def _norm_model(value: str) -> str:
    if not value:
        return ""
    lower = value.strip().lower()
    if lower in ("g4emlivermorephysics", "livermore", "livermorephysics"):
        return "livermore"
    if lower in ("g4empenelopephysics", "penelope", "penelopephysics"):
        return "penelope"
    if lower in ("pai", "g4pai", "g4empai"):
        return "pai"
    return lower


def dionne_delta(Ep_eV: float, params: DionneParams, d_mode: str = "A_pow") -> float:
    if Ep_eV <= 0:
        return 0.0
    if d_mode == "A_linear":
        # alternative: d = E^n / A
        d = (Ep_eV ** params.n) / params.A
    else:
        # paper Eq. 5: d = E^n / A^n
        d = (Ep_eV ** params.n) / (params.A ** params.n)
    term1 = (params.B / params.Ea)
    term2 = ((params.A * params.n) / params.alpha) ** (1.0 / params.n)
    term3 = (params.alpha * d) ** (1.0 / params.n - 1.0)
    term4 = 1.0 - math.exp(-params.alpha * d)
    return term1 * term2 * term3 * term4


def fit_ea_for_peak(params: DionneParams, e_peak: float, delta_peak: float, d_mode: str) -> float:
    if delta_peak <= 0 or e_peak <= 0:
        return params.Ea
    unit_params = DionneParams(
        B=params.B, A=params.A, n=params.n, alpha=params.alpha, Ea=1.0, label=params.label
    )
    delta_unit = dionne_delta(e_peak, unit_params, d_mode=d_mode)
    if delta_unit <= 0:
        return params.Ea
    return delta_unit / delta_peak


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
    parser.add_argument("--em-model", default=None,
                        help="Require a specific EM model (PAI/Livermore/Penelope)")
    parser.add_argument("--mc-source", choices=["geant4", "toy"], default="geant4",
                        help="Source for MC points: geant4 (RunMeta.sey) or toy (SEY_MonteCarlo histogram)")
    parser.add_argument("--toy-hist", default="SEY_MonteCarlo",
                        help="Histogram name in toy MC output (default: SEY_MonteCarlo)")
    parser.add_argument("--fit-fig9", action="store_true",
                        help="Fit Ea to Fig. 9 peak (default: 375 eV, δ=4.08) and plot curve")
    parser.add_argument("--fig9-ep", type=float, default=375.0,
                        help="Fig. 9 peak energy in eV (default: 375)")
    parser.add_argument("--fig9-delta", type=float, default=4.08,
                        help="Fig. 9 peak SEY (default: 4.08)")
    parser.add_argument("--fig9-material", default="al2o3",
                        help="Preset for Fig. 9 fit parameters (default: al2o3)")
    parser.add_argument("--fit-mc-peak", action="store_true",
                        help="Fit Ea to MC peak and plot curve (uses current --material parameters)")
    parser.add_argument("--fig9-ea-on-material", action="store_true",
                        help="When --fit-fig9 is set, also apply the Fig. 9-fitted Ea to the "
                             "current material parameters (useful for Fig. 10 thickness curves).")
    parser.add_argument("--d-mode", choices=["A_linear", "A_pow"], default="A_pow",
                        help="Penetration depth formula: A_pow uses d=E^n/A^n (paper Eq. 5), A_linear uses d=E^n/A")
    parser.add_argument("--alpha-unit", choices=["nm", "angstrom"], default="angstrom",
                        help="Unit used for alpha in the paper table (default: angstrom). "
                             "angstrom values are converted to 1/nm by multiplying by 10.")
    parser.add_argument("--alpha-scale", type=float, default=None,
                        help="Override alpha scaling factor (applied after --alpha overrides). "
                             "Example: 10 for Angstrom->nm conversion.")
    parser.add_argument("--output-dir", default="plots/MC_electrons_on_shell_dionne-model/mc_vs_dionne",
                        help="Output directory for plots")
    parser.add_argument("--title", default="SEY validation vs Dionne model", help="Plot title")
    args = parser.parse_args()

    if args.material not in PRESETS:
        raise SystemExit(f"Unknown material preset '{args.material}'. Valid: {', '.join(PRESETS.keys())}")
    if args.fig9_material not in PRESETS:
        raise SystemExit(f"Unknown fig9 material preset '{args.fig9_material}'. Valid: {', '.join(PRESETS.keys())}")

    params = PRESETS[args.material]
    params = DionneParams(
        B=params.B if args.B is None else args.B,
        A=params.A if args.A is None else args.A,
        n=params.n if args.n is None else args.n,
        alpha=params.alpha if args.alpha is None else args.alpha,
        Ea=params.Ea if args.Ea is None else args.Ea,
        label=params.label,
    )
    alpha_raw = params.alpha
    if args.alpha_scale is not None:
        alpha_scale = args.alpha_scale
    else:
        alpha_scale = 10.0 if args.alpha_unit == "angstrom" else 1.0
    params = DionneParams(
        B=params.B,
        A=params.A,
        n=params.n,
        alpha=params.alpha * alpha_scale,
        Ea=params.Ea,
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
    expected_model = _norm_model(args.em_model) if args.em_model else ""

    if args.mc_source == "toy":
        root_files = [
            os.path.join(args.results_dir, f)
            for f in os.listdir(args.results_dir)
            if f.endswith("_SEY_MonteCarlo.root")
        ]
    else:
        root_files = [
            os.path.join(args.results_dir, f)
            for f in os.listdir(args.results_dir)
            if f.endswith(".root") and not f.startswith("summary") and "_SEY_MonteCarlo" not in f
        ]
    if not root_files:
        raise SystemExit(f"No ROOT files found in {args.results_dir}")

    for path in sorted(root_files):
        f = ROOT.TFile.Open(path)
        if not f or not f.IsOpen():
            continue
        try:
            sample_thickness_val = None
            em_model_val = None
            particle_val = None
            if args.mc_source == "toy":
                toy_hist = f.Get(args.toy_hist)
                if not toy_hist:
                    f.Close()
                    continue
                sey = float(toy_hist.GetMean())
                n_events = float(toy_hist.GetEntries())
                rms = float(toy_hist.GetRMS())
                err = (rms / math.sqrt(n_events)) if n_events > 0 else 0.0

                base_path = path.replace("_SEY_MonteCarlo.root", ".root")
                if not os.path.isfile(base_path):
                    f.Close()
                    continue
                base = ROOT.TFile.Open(base_path)
                meta = base.Get("RunMeta")
                if not meta:
                    base.Close()
                    f.Close()
                    continue
                meta.GetEntry(0)
                energy_mev = float(meta.primaryEnergyMeV)
                energy_eV = energy_mev * 1.0e6
                n_primary = int(meta.primaryElectrons)
                sample_thickness_val = float(meta.sampleThicknessNm)
                em_model_val = _as_str(meta.emModel)
                particle_val = _as_str(meta.primaryParticle)
                base.Close()
            else:
                meta = f.Get("RunMeta")
                if not meta:
                    f.Close()
                    continue
                meta.GetEntry(0)
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
                sample_thickness_val = float(meta.sampleThicknessNm)
                em_model_val = _as_str(meta.emModel)
                particle_val = _as_str(meta.primaryParticle)
            energies_eV.append(energy_eV)
            sey_values.append(sey)
            sey_errors.append(err)
            if thickness_nm is None and sample_thickness_val is not None:
                thickness_nm = sample_thickness_val
            if em_model is None and em_model_val is not None:
                em_model = em_model_val
            if expected_model:
                if _norm_model(em_model) != expected_model:
                    raise SystemExit(
                        f"RunMeta emModel={em_model} does not match requested {args.em_model}."
                    )
            if particle_val not in ("e-", "e+"):
                raise SystemExit(
                    f"RunMeta primaryParticle={particle_val} is not an electron. "
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

    # Analytic curve(s)
    min_e = min(energies_eV)
    max_e = max(energies_eV)
    n_curve = 400
    step = (max_e - min_e) / (n_curve - 1)
    curve_x = [min_e + i * step for i in range(n_curve)]

    fit_curves = []
    fig9_ea = None
    if args.fit_fig9:
        fig_params = PRESETS[args.fig9_material]
        fig_params_scaled = DionneParams(
            B=fig_params.B, A=fig_params.A, n=fig_params.n,
            alpha=fig_params.alpha * alpha_scale, Ea=fig_params.Ea, label=fig_params.label
        )
        fig9_ea = fit_ea_for_peak(fig_params_scaled, args.fig9_ep, args.fig9_delta, d_mode=args.d_mode)
        fig_params = DionneParams(
            B=fig_params.B, A=fig_params.A, n=fig_params.n,
            alpha=fig_params.alpha * alpha_scale, Ea=fig9_ea, label=fig_params.label
        )
        curve_y = [dionne_delta(e, fig_params, d_mode=args.d_mode) for e in curve_x]
        fit_curves.append(("Fig. 9 fit", fig_params, curve_y))
        if args.fig9_ea_on_material:
            fig9_on_mat = DionneParams(
                B=params.B, A=params.A, n=params.n,
                alpha=params.alpha, Ea=fig9_ea, label=params.label
            )
            curve_y = [dionne_delta(e, fig9_on_mat, d_mode=args.d_mode) for e in curve_x]
            fit_curves.append(("Fig. 9 Ea on material", fig9_on_mat, curve_y))

    if args.fit_mc_peak:
        # Use current MC peak as target
        idx_peak = max(range(len(sey_values)), key=lambda i: sey_values[i])
        e_peak_mc = energies_eV[idx_peak]
        delta_peak_mc = sey_values[idx_peak]
        ea_mc = fit_ea_for_peak(params, e_peak_mc, delta_peak_mc, d_mode=args.d_mode)
        mc_params = DionneParams(
            B=params.B, A=params.A, n=params.n,
            alpha=params.alpha, Ea=ea_mc, label=params.label
        )
        curve_y = [dionne_delta(e, mc_params, d_mode=args.d_mode) for e in curve_x]
        fit_curves.append(("MC peak fit", mc_params, curve_y))

    if not fit_curves:
        curve_y = [dionne_delta(e, params, d_mode=args.d_mode) for e in curve_x]
        fit_curves.append(("Dionne model", params, curve_y))

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

    colors = [ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1]
    styles = [1, 2, 7]
    analytic_graphs = []
    for idx, (label, pfit, yvals) in enumerate(fit_curves):
        gr = ROOT.TGraph(len(curve_x))
        for i, (x, y) in enumerate(zip(curve_x, yvals)):
            gr.SetPoint(i, x, y)
        gr.SetLineWidth(3)
        gr.SetLineColor(colors[idx % len(colors)])
        gr.SetLineStyle(styles[idx % len(styles)])
        gr.Draw("L")
        analytic_graphs.append((label, pfit, gr))

    leg = ROOT.TLegend(0.58, 0.70, 0.92, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    mc_label = "Toy MC SEY" if args.mc_source == "toy" else "MC SEY"
    leg.AddEntry(gr_mc, mc_label, "p")
    for label, pfit, gr in analytic_graphs:
        leg.AddEntry(gr, f"{label} (E_a={pfit.Ea:.2g} eV)", "l")
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
    text.AddText(f"B={params.B:.3g}, A={params.A:.3g}, n={params.n:.3g}, #alpha={params.alpha:.3g} 1/nm")
    if alpha_scale != 1.0:
        text.AddText(f"#alpha (table)={alpha_raw:.3g} 1/#AA  (x{alpha_scale:.0f})")
    text.AddText(f"E_a={params.Ea:.3g} eV")
    text.Draw()

    if args.fit_fig9 or args.fit_mc_peak:
        fit_box = ROOT.TPaveText(0.58, 0.52, 0.92, 0.70, "NDC")
        fit_box.SetFillColor(0)
        fit_box.SetBorderSize(1)
        fit_box.SetTextFont(42)
        fit_box.SetTextSize(0.024)
        fit_box.AddText("Dionne fits")
        if args.fit_fig9:
            fit_box.AddText(f"Fig. 9: E_p={args.fig9_ep:.0f} eV, #delta={args.fig9_delta:.2f}")
        if args.fit_mc_peak:
            idx_peak = max(range(len(sey_values)), key=lambda i: sey_values[i])
            fit_box.AddText(f"MC peak: E_p={energies_eV[idx_peak]:.0f} eV, #delta={sey_values[idx_peak]:.2f}")
        fit_box.Draw()

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
    model_tag = _norm_model(em_model) if em_model else ""
    model_suffix = f"_{model_tag}" if model_tag else ""
    source_suffix = "_toy" if args.mc_source == "toy" else ""
    fit_suffix = ""
    if args.fit_fig9 or args.fit_mc_peak:
        fit_suffix = "_fit"
    out_base = os.path.join(
        args.output_dir,
        f"sey_dionne_vs_mc_{args.material}_{thick_tag}{model_suffix}{source_suffix}{fit_suffix}",
    )
    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".root")


if __name__ == "__main__":
    main()
