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
    "al2o3": DionneParams(B=0.46, A=37.0, n=1.61, alpha=0.0075, Ea=10.0, label="Al2O3 (Table II)"),
    "si": DionneParams(B=0.26, A=20.0, n=1.38, alpha=0.040, Ea=10.0, label="Si (Table II)"),
    "al2o3_1nm": DionneParams(B=0.382, A=25.0, n=1.45, alpha=0.030, Ea=10.0, label="Al2O3 1 nm (Table III)"),
    "al2o3_3nm": DionneParams(B=0.425, A=35.0, n=1.57, alpha=0.013, Ea=10.0, label="Al2O3 3 nm (Table III)"),
    "al2o3_5nm": DionneParams(B=0.450, A=36.0, n=1.60, alpha=0.0090, Ea=10.0, label="Al2O3 5 nm (Table III)"),
}


def _apply_alpha_scale(params: DionneParams, alpha_scale: float) -> DionneParams:
    return DionneParams(
        B=params.B,
        A=params.A,
        n=params.n,
        alpha=params.alpha * alpha_scale,
        Ea=params.Ea,
        label=params.label,
    )


def _mix_params(
    thickness_nm: float,
    energy_eV: float,
    al2o3_params: DionneParams,
    si_params: DionneParams,
    max_iter: int = 3,
) -> DionneParams:
    """
    Mixed-layer parameters for Al2O3 on Si based on penetration depth.

    Assumption (from paper text): proportions are given by the fraction of the
    penetration depth within the Al2O3 layer. We approximate:
        f_al2o3 = min(1, thickness_nm / d)
        f_si = 1 - f_al2o3
    and linearly mix B, A, n, alpha, Ea.

    We iterate a few times because d depends on A and n.
    """
    if energy_eV <= 0 or thickness_nm <= 0:
        return al2o3_params
    # Start from Al2O3 parameters for d
    mixed = al2o3_params
    for _ in range(max_iter):
        # Dionne Eq. 5: d = E^n / A^n
        d = (energy_eV ** mixed.n) / (mixed.A ** mixed.n)
        if d <= 0:
            f_al = 1.0
        else:
            f_al = min(1.0, thickness_nm / d)
        f_si = 1.0 - f_al
        mixed = DionneParams(
            B=f_al * al2o3_params.B + f_si * si_params.B,
            A=f_al * al2o3_params.A + f_si * si_params.A,
            n=f_al * al2o3_params.n + f_si * si_params.n,
            alpha=f_al * al2o3_params.alpha + f_si * si_params.alpha,
            Ea=f_al * al2o3_params.Ea + f_si * si_params.Ea,
            label="Al2O3/Si mixed",
        )
    return mixed


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
    parser.add_argument("--fit-mc-peak", action="store_true",
                        help="Fit Ea to MC peak and plot curve (uses current --material parameters)")
    parser.add_argument("--d-mode", choices=["A_linear", "A_pow"], default="A_pow",
                        help="Penetration depth formula: A_pow uses d=E^n/A^n (paper Eq. 5), A_linear uses d=E^n/A")
    parser.add_argument("--dionne-model", choices=["intrinsic", "mixed"], default="intrinsic",
                        help="Use intrinsic (Table II/III) parameters or mixed-layer model (Al2O3/Si).")
    parser.add_argument("--alpha-unit", choices=["nm", "angstrom"], default="nm",
                        help="Unit used for alpha in the paper table (default: nm). "
                             "angstrom values are converted to 1/nm by multiplying by 10.")
    parser.add_argument("--alpha-scale", type=float, default=None,
                        help="Override alpha scaling factor (applied after --alpha overrides). "
                             "Example: 10 for Angstrom->nm conversion.")
    parser.add_argument("--paper-curve", default=None,
                        help="Optional CSV with digitized paper curve (energy_eV,sey).")
    parser.add_argument("--paper-curve-label", default="Paper curve (digitized)",
                        help="Legend label for --paper-curve.")
    parser.add_argument("--fit-dionne", action="store_true",
                        help="Fit MC points to Dionne model and report fitted parameters.")
    parser.add_argument("--fit-ea", action="store_true",
                        help="Include E_a as a free parameter in the Dionne fit.")
    parser.add_argument("--fit-n-min", type=float, default=None,
                        help="Allow n to vary with lower bound (e.g. 1.3).")
    parser.add_argument("--fit-n-max", type=float, default=None,
                        help="Allow n to vary with upper bound (e.g. 2.0).")
    parser.add_argument("--fit-output", default=None,
                        help="Optional path to save fitted parameters as a text report.")
    parser.add_argument("--output-dir", default="plots/MC_electrons_on_shell_dionne-model/mc_vs_dionne",
                        help="Output directory for plots")
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
    alpha_raw = params.alpha
    if args.alpha_scale is not None:
        alpha_scale = args.alpha_scale
    else:
        alpha_scale = 10.0 if args.alpha_unit == "angstrom" else 1.0
    params = _apply_alpha_scale(params, alpha_scale)

    # Intrinsic parameters for mixed-layer model
    al2o3_intrinsic = _apply_alpha_scale(PRESETS["al2o3"], alpha_scale)
    si_intrinsic = _apply_alpha_scale(PRESETS["si"], alpha_scale)

    if args.dionne_model == "mixed" and args.fit_mc_peak:
        raise SystemExit("--fit-mc-peak is not supported for the mixed-layer model.")
    if args.dionne_model == "mixed" and args.fit_dionne:
        raise SystemExit("--fit-dionne is currently supported only for intrinsic parameters.")

    import ROOT
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    # Collect MC points
    energies_eV = []
    sey_values = []
    sey_errors = []
    edep_primary_means = []
    edep_primary_errors = []
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
            edep_primary_mean = float("nan")
            edep_primary_err = 0.0
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
                h_edep = base.Get("EdepPrimary")
                if h_edep:
                    edep_primary_mean = float(h_edep.GetMean())
                    n_edep = float(h_edep.GetEntries())
                    if n_edep > 0:
                        edep_primary_err = float(h_edep.GetRMS()) / math.sqrt(n_edep)
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
                h_edep = f.Get("EdepPrimary")
                if h_edep:
                    edep_primary_mean = float(h_edep.GetMean())
                    n_edep = float(h_edep.GetEntries())
                    if n_edep > 0:
                        edep_primary_err = float(h_edep.GetRMS()) / math.sqrt(n_edep)
            energies_eV.append(energy_eV)
            sey_values.append(sey)
            sey_errors.append(err)
            edep_primary_means.append(edep_primary_mean)
            edep_primary_errors.append(edep_primary_err)
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
    edep_primary_means = [edep_primary_means[i] for i in order]
    edep_primary_errors = [edep_primary_errors[i] for i in order]

    # Analytic curve(s)
    min_e = min(energies_eV)
    max_e = max(energies_eV)
    n_curve = 400
    step = (max_e - min_e) / (n_curve - 1)
    curve_x = [min_e + i * step for i in range(n_curve)]

    fit_curves = []
    if args.dionne_model == "mixed":
        if thickness_nm is None:
            raise SystemExit("Mixed-layer model requires thickness from RunMeta.")
        base_curve_y = [
            dionne_delta(
                e,
                _mix_params(thickness_nm, e, al2o3_intrinsic, si_intrinsic),
                d_mode=args.d_mode,
            )
            for e in curve_x
        ]
        fit_curves.append(("Dionne mixed-layer model", None, base_curve_y))
    else:
        # Base Dionne model (Table parameters) only when not fitting MC
        if not args.fit_dionne:
            base_curve_y = [dionne_delta(e, params, d_mode=args.d_mode) for e in curve_x]
            fit_curves.append(("Dionne model", params, base_curve_y))

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

    # Optional paper curve
    paper_graph = None
    paper_max = None
    if args.paper_curve:
        try:
            xs = []
            ys = []
            with open(args.paper_curve, "r") as f:
                for line in f:
                    if line.strip().startswith("energy_eV"):
                        continue
                    parts = line.strip().split(",")
                    if len(parts) < 2:
                        continue
                    xs.append(float(parts[0]))
                    ys.append(float(parts[1]))
            if xs:
                paper_graph = ROOT.TGraph(len(xs))
                for i, (x, y) in enumerate(zip(xs, ys)):
                    paper_graph.SetPoint(i, x, y)
                paper_graph.SetLineWidth(2)
                paper_graph.SetLineColor(ROOT.kGray + 2)
                paper_graph.SetLineStyle(7)
                paper_max = max(ys)
        except Exception:
            paper_graph = None

    # Plot
    if fit_curves:
        max_curve = max(max(yvals) for _, _, yvals in fit_curves)
    else:
        max_curve = max(sey_values)
    if paper_max is not None:
        max_curve = max(max_curve, paper_max)
    c = ROOT.TCanvas("c_sey_dionne", "SEY vs Dionne", 1000, 700)
    c.SetTopMargin(0.28)
    frame = ROOT.TH2D("frame", "", 10, min_e * 0.95, max_e * 1.05, 10, 0.0, max(max(sey_values), max_curve) * 1.25)
    frame.SetXTitle("Primary electron energy (eV)")
    frame.SetYTitle("Secondary electron yield (SEY)")
    frame.Draw("AXIS")
    y1_min = frame.GetYaxis().GetXmin()
    y1_max = frame.GetYaxis().GetXmax()

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

    # Overlay mean EdepPrimary on a right-hand axis.
    # Values are mapped to the left SEY axis for drawing.
    gr_edep = None
    edep_axis = None
    valid_edep = [
        i for i, v in enumerate(edep_primary_means)
        if isinstance(v, (int, float)) and math.isfinite(v)
    ]
    if valid_edep:
        edep_max = max(edep_primary_means[i] for i in valid_edep)
        edep_axis_min = 0.0
        edep_axis_max = edep_max * 1.15 if edep_max > 0 else 1.0
        edep_span = edep_axis_max - edep_axis_min
        y1_span = y1_max - y1_min

        gr_edep = ROOT.TGraphErrors(len(valid_edep))
        for ip, i in enumerate(valid_edep):
            x = energies_eV[i]
            v = edep_primary_means[i]
            ve = edep_primary_errors[i] if i < len(edep_primary_errors) else 0.0
            y = y1_min + ((v - edep_axis_min) / edep_span) * y1_span
            ye = (ve / edep_span) * y1_span
            gr_edep.SetPoint(ip, x, y)
            gr_edep.SetPointError(ip, 0.0, ye)
        gr_edep.SetLineColor(ROOT.kGray + 2)
        gr_edep.SetMarkerColor(ROOT.kGray + 2)
        gr_edep.SetLineWidth(2)
        gr_edep.SetLineStyle(2)
        gr_edep.SetMarkerStyle(24)
        gr_edep.SetMarkerSize(1.0)
        gr_edep.Draw("PL")

        edep_axis = ROOT.TGaxis(
            max_e * 1.05, y1_min, max_e * 1.05, y1_max,
            edep_axis_min, edep_axis_max, 510, "+L"
        )
        edep_axis.SetTitle("Mean EdepPrimary per event (eV)")
        edep_axis.SetTitleOffset(1.15)
        edep_axis.SetLabelFont(42)
        edep_axis.SetTitleFont(42)
        edep_axis.SetLineColor(ROOT.kGray + 2)
        edep_axis.SetLabelColor(ROOT.kGray + 2)
        edep_axis.SetTitleColor(ROOT.kGray + 2)
        edep_axis.Draw()

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

    fit_result = None
    fit_params = None
    if args.fit_dionne:
        def _dionne_tf1(x, p):
            E = x[0]
            if E <= 0:
                return 0.0
            B, A, n, alpha, Ea = p[0], p[1], p[2], p[3], p[4]
            if A <= 0 or n <= 0 or alpha <= 0 or Ea <= 0:
                return 0.0
            if args.d_mode == "A_linear":
                d = (E ** n) / A
            else:
                d = (E ** n) / (A ** n)
            term1 = (B / Ea)
            term2 = ((A * n) / alpha) ** (1.0 / n)
            term3 = (alpha * d) ** (1.0 / n - 1.0)
            term4 = 1.0 - math.exp(-alpha * d)
            return term1 * term2 * term3 * term4

        tf1 = ROOT.TF1("dionne_fit", _dionne_tf1, min_e, max_e, 5)
        tf1.SetParNames("B", "A", "n", "alpha", "Ea")
        tf1.SetParameters(params.B, params.A, params.n, params.alpha, params.Ea)
        tf1.SetParLimits(0, 0.0, 5.0)
        tf1.SetParLimits(1, 1.0, 200.0)
        tf1.SetParLimits(2, 0.5, 3.0)
        tf1.SetParLimits(3, 0.001, 1.0)
        tf1.SetParLimits(4, 0.1, 50.0)
        if args.fit_n_min is not None or args.fit_n_max is not None:
            n_min = args.fit_n_min if args.fit_n_min is not None else 0.5
            n_max = args.fit_n_max if args.fit_n_max is not None else 3.0
            tf1.SetParLimits(2, n_min, n_max)
        else:
            # If fitting MC, keep n fixed to the preset value unless user loosens it.
            # This stabilizes the fit against over-flexible n.
            tf1.FixParameter(2, params.n)
        if not args.fit_ea:
            tf1.FixParameter(4, params.Ea)
        fit_result = gr_mc.Fit(tf1, "QS")
        fit_params = DionneParams(
            B=tf1.GetParameter(0),
            A=tf1.GetParameter(1),
            n=tf1.GetParameter(2),
            alpha=tf1.GetParameter(3),
            Ea=tf1.GetParameter(4),
            label="Dionne fit",
        )
        fit_graph = ROOT.TGraph(len(curve_x))
        for i, x in enumerate(curve_x):
            fit_graph.SetPoint(i, x, tf1.Eval(x))
        fit_graph.SetLineWidth(3)
        fit_graph.SetLineColor(ROOT.kOrange + 7)
        fit_graph.SetLineStyle(1)
        fit_graph.Draw("L")
        analytic_graphs.append(("Dionne fit", fit_params, fit_graph))

    if paper_graph:
        paper_graph.Draw("L")

    leg = ROOT.TLegend(0.58, 0.58, 0.92, 0.72)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    mc_label = "Toy MC SEY" if args.mc_source == "toy" else "MC SEY"
    leg.AddEntry(gr_mc, mc_label, "p")
    for label, pfit, gr in analytic_graphs:
        if pfit is None:
            leg.AddEntry(gr, label, "l")
        else:
            leg.AddEntry(gr, f"{label} (E_a={pfit.Ea:.2g} eV)", "l")
    if gr_edep:
        leg.AddEntry(gr_edep, "Mean EdepPrimary (right axis)", "pl")
    if paper_graph:
        leg.AddEntry(paper_graph, args.paper_curve_label, "l")
    leg.Draw()

    text = ROOT.TPaveText(0.12, 0.78, 0.54, 0.96, "NDC")
    text.SetFillColor(0)
    text.SetBorderSize(1)
    text.SetTextFont(42)
    text.SetTextSize(0.028)
    if thickness_nm is not None:
        text.AddText(f"Thickness: {thickness_nm:.0f} nm")
    if em_model:
        text.AddText(f"EM model: {em_model}")
    if args.dionne_model == "mixed":
        text.AddText("Mixed-layer model (Al2O3 + Si)")
        text.AddText("f_{Al2O3} = min(1, t/d), d=E^{n}/A^{n}")
        text.AddText(
            f"Al2O3: B={al2o3_intrinsic.B:.3g}, A={al2o3_intrinsic.A:.3g}, "
            f"n={al2o3_intrinsic.n:.3g}, #alpha={al2o3_intrinsic.alpha:.3g} 1/nm"
        )
        text.AddText(
            f"Si: B={si_intrinsic.B:.3g}, A={si_intrinsic.A:.3g}, "
            f"n={si_intrinsic.n:.3g}, #alpha={si_intrinsic.alpha:.3g} 1/nm"
        )
        text.AddText(f"E_a={al2o3_intrinsic.Ea:.3g} eV")
    else:
        if not args.fit_dionne:
            text.AddText(params.label)
            text.AddText(
                f"B={params.B:.3g}, A={params.A:.3g}, n={params.n:.3g}, #alpha={params.alpha:.3g} 1/nm"
            )
            if alpha_scale != 1.0:
                text.AddText(f"#alpha (table)={alpha_raw:.3g} 1/#AA  (x{alpha_scale:.0f})")
            text.AddText(f"E_a={params.Ea:.3g} eV")
    if fit_params is not None:
        text.AddText(
            f"Fit: B={fit_params.B:.3g}, A={fit_params.A:.3g}, n={fit_params.n:.3g}, "
            f"#alpha={fit_params.alpha:.3g} 1/nm"
        )
        text.AddText(f"Fit: E_a={fit_params.Ea:.3g} eV")
    text.Draw()

    if args.fit_mc_peak:
        fit_box = ROOT.TPaveText(0.58, 0.78, 0.92, 0.96, "NDC")
        fit_box.SetFillColor(0)
        fit_box.SetBorderSize(1)
        fit_box.SetTextFont(42)
        fit_box.SetTextSize(0.024)
        fit_box.AddText("Dionne fits")
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
        0.50,
        "#delta = (B/E_{a}) (A n/#alpha)^{1/n} (#alpha d)^{1/n-1} (1-e^{-#alpha d})",
    )

    os.makedirs(args.output_dir, exist_ok=True)
    thick_tag = f"{int(round(thickness_nm))}nm" if thickness_nm is not None else "unknown"
    model_tag = _norm_model(em_model) if em_model else ""
    model_suffix = f"_{model_tag}" if model_tag else ""
    source_suffix = "_toy" if args.mc_source == "toy" else ""
    fit_suffix = ""
    if args.fit_mc_peak:
        fit_suffix = "_fit"
    out_base = os.path.join(
        args.output_dir,
        f"sey_dionne_vs_mc_{args.material}_{thick_tag}{model_suffix}{source_suffix}{fit_suffix}",
    )
    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".root")

    if args.fit_dionne:
        report_lines = []
        report_lines.append("Dionne fit results")
        report_lines.append(f"Material: {args.material}")
        report_lines.append(f"Thickness: {thickness_nm} nm")
        if em_model:
            report_lines.append(f"EM model: {em_model}")
        if fit_result:
            report_lines.append(f"chi2/ndf = {fit_result.Chi2():.3g}/{int(fit_result.Ndf())}")
        if fit_params:
            report_lines.append(f"B = {fit_params.B:.6g}")
            report_lines.append(f"A = {fit_params.A:.6g}")
            report_lines.append(f"n = {fit_params.n:.6g}")
            report_lines.append(f"alpha = {fit_params.alpha:.6g} 1/nm")
            report_lines.append(f"E_a = {fit_params.Ea:.6g} eV")
        report_text = "\n".join(report_lines) + "\n"
        if args.fit_output:
            with open(args.fit_output, "w") as f:
                f.write(report_text)
        else:
            print(report_text)


if __name__ == "__main__":
    main()
