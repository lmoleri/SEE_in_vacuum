#!/usr/bin/env python3
import argparse
import os


def _format_energy_ev(value):
    if abs(value - round(value)) < 1e-6:
        return f"{int(round(value))}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def _meta_string(value):
    if value is None:
        return ""
    try:
        if isinstance(value, (bytes, bytearray)):
            raw = bytes(value)
        elif value.__class__.__name__ == "LowLevelView":
            raw = bytes(value)
        else:
            return str(value).split("\x00", 1)[0].strip()
        return raw.split(b"\x00", 1)[0].decode("utf-8", errors="ignore").strip()
    except Exception:
        return str(value).split("\x00", 1)[0].strip()


def _load_curves(results_dir, hist_name):
    import ROOT

    root_files = [
        f for f in os.listdir(results_dir)
        if f.endswith(".root") and "_SEY_MonteCarlo" not in f and not f.startswith("summary")
    ]
    root_files.sort()

    curves = []
    thickness_nm = None
    em_model = None
    n_events = None

    for name in root_files:
        path = os.path.join(results_dir, name)
        f = ROOT.TFile.Open(path)
        if not f or not f.IsOpen():
            continue
        meta = f.Get("RunMeta")
        if not meta:
            f.Close()
            continue
        meta.GetEntry(0)
        energy_ev = float(meta.primaryEnergyMeV) * 1.0e6
        n_events = int(meta.primaryElectrons)
        thickness_nm = float(meta.sampleThicknessNm)
        em_model = _meta_string(meta.emModel)

        h = f.Get(hist_name)
        if not h:
            f.Close()
            raise SystemExit(
                f"Missing {hist_name} in {path}. Re-run the Geant4 scan."
            )
        hc = h.Clone()
        hc.SetDirectory(0)
        curves.append((energy_ev, hc))
        f.Close()

    curves.sort(key=lambda x: x[0])
    return curves, thickness_nm, em_model, n_events


def _draw_overlay(
    curves,
    title,
    x_title,
    y_title,
    out_base,
    thickness_nm,
    em_model,
    n_events,
    root_file=None,
    canvas_key=None,
):
    import ROOT

    colors = [
        ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1,
        ROOT.kOrange + 7, ROOT.kCyan + 1, ROOT.kViolet + 1, ROOT.kAzure + 2,
        ROOT.kTeal + 1, ROOT.kPink + 2, ROOT.kGray + 2
    ]

    c = ROOT.TCanvas("c_overlay", "Depth overlay", 1100, 750)
    c.SetRightMargin(0.28)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)

    max_y = 0.0
    x_min = float("inf")
    x_max = float("-inf")
    for _, h in curves:
        max_y = max(max_y, h.GetMaximum())
        x_min = min(x_min, h.GetXaxis().GetXmin())
        x_max = max(x_max, h.GetXaxis().GetXmax())
    if max_y <= 0:
        max_y = 1.0
    if x_min >= x_max:
        x_min, x_max = 0.0, 1.0

    frame = ROOT.TH2D(
        f"frame_{abs(hash(out_base)) % 100000000}",
        "",
        10,
        x_min,
        x_max,
        10,
        0.0,
        max_y * 1.15,
    )
    frame.SetTitle(title)
    frame.SetXTitle(x_title)
    frame.SetYTitle(y_title)
    frame.Draw("AXIS")

    leg = ROOT.TLegend(0.74, 0.14, 0.98, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(1)

    for idx, (energy_ev, h) in enumerate(curves):
        color = colors[idx % len(colors)]
        h.SetLineWidth(2)
        h.SetLineColor(color)
        h.SetMarkerColor(color)
        h.Draw("HIST SAME")
        leg.AddEntry(h, f"{_format_energy_ev(energy_ev)} eV", "l")

    leg.Draw()

    text = ROOT.TPaveText(0.12, 0.84, 0.58, 0.96, "NDC")
    text.SetFillColor(0)
    text.SetBorderSize(1)
    text.SetTextFont(42)
    text.SetTextSize(0.028)
    if thickness_nm is not None:
        text.AddText(f"Thickness: {thickness_nm:.0f} nm")
    if em_model:
        text.AddText(f"EM model: {em_model}")
    if n_events is not None:
        text.AddText(f"MC events: {n_events}")
    text.Draw()

    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".png")
    if root_file:
        root_file.cd()
        key = canvas_key if canvas_key else os.path.basename(out_base)
        c.Write(key)
    c.Close()


def main():
    parser = argparse.ArgumentParser(
        description="Overlay depth distributions across primary energies."
    )
    parser.add_argument("--results-dir", required=True, help="Results directory with ROOT files.")
    parser.add_argument(
        "--output-dir",
        default="plots/MC_electrons_on_shell_dionne-model/diagnostics_edep_depth/overlays",
        help="Output directory for plots.",
    )
    parser.add_argument("--label", default=None, help="Label for output filenames.")
    parser.add_argument(
        "--save-results",
        action="store_true",
        help="Also save plots under <results-dir>/plots.",
    )
    args = parser.parse_args()

    import ROOT
    ROOT.TH1.AddDirectory(False)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    os.makedirs(args.output_dir, exist_ok=True)
    label = args.label or os.path.basename(os.path.normpath(args.results_dir))
    output_root_path = os.path.join(args.output_dir, f"overlays_{label}.root")
    output_root = ROOT.TFile.Open(output_root_path, "RECREATE")
    if not output_root or output_root.IsZombie():
        raise SystemExit(f"Cannot create ROOT output file: {output_root_path}")

    results_root = None
    results_plot_dir = None
    if args.save_results:
        results_plot_dir = os.path.join(args.results_dir, "plots")
        os.makedirs(results_plot_dir, exist_ok=True)
        results_root_path = os.path.join(results_plot_dir, f"overlays_{label}.root")
        results_root = ROOT.TFile.Open(results_root_path, "RECREATE")
        if not results_root or results_root.IsZombie():
            raise SystemExit(f"Cannot create ROOT output file: {results_root_path}")

    specs = [
        (
            "EdepDepthPrimary",
            "Primary e- energy deposition vs depth in Al_{2}O_{3}",
            "Depth from entrance (nm)",
            "Energy deposition (eV)",
            f"edep_depth_overlay_{label}",
        ),
        (
            "EdepDepthPrimaryWeighted",
            "Depth-weighted primary e- energy deposition vs depth in Al_{2}O_{3}",
            "Depth from entrance (nm)",
            "Weighted energy deposition (eV)",
            f"edep_depth_weighted_overlay_{label}",
        ),
        (
            "EdepDepthPrimaryCounts",
            "Primary e- energy-depositing step counts vs depth in Al_{2}O_{3}",
            "Depth from entrance (nm)",
            "Energy-depositing steps per event",
            f"edep_depth_counts_overlay_{label}",
        ),
        (
            "EdepPrimary",
            "Primary e- total energy deposition in Al_{2}O_{3}",
            "Primary energy deposition (eV)",
            "Number of events",
            f"edep_primary_overlay_{label}",
        ),
    ]

    for hist_name, title, x_title, y_title, out_name in specs:
        curves, thickness_nm, em_model, n_events = _load_curves(args.results_dir, hist_name)
        if not curves:
            raise SystemExit(f"No curves loaded for {hist_name}")
        out_base = os.path.join(args.output_dir, out_name)
        _draw_overlay(
            curves,
            title,
            x_title,
            y_title,
            out_base,
            thickness_nm,
            em_model,
            n_events,
            root_file=output_root,
            canvas_key=out_name,
        )

        if args.save_results:
            out_base_results = os.path.join(results_plot_dir, out_name)
            _draw_overlay(
                curves,
                title,
                x_title,
                y_title,
                out_base_results,
                thickness_nm,
                em_model,
                n_events,
                root_file=results_root,
                canvas_key=out_name,
            )

    output_root.Close()
    if results_root:
        results_root.Close()

    print(f"Saved overlay plots in {args.output_dir}")
    print(f"Saved overlay canvases ROOT file: {output_root_path}")
    if args.save_results:
        print(f"Saved results overlay canvases ROOT file: {results_root_path}")


if __name__ == "__main__":
    main()
