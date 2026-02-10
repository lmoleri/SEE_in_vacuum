#!/usr/bin/env python3
import argparse
import os

def _format_energy_ev(value):
    if abs(value - round(value)) < 1e-6:
        return f"{int(round(value))}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def main():
    parser = argparse.ArgumentParser(description="Overlay EdepDepthPrimaryCounts curves by energy.")
    parser.add_argument("--results-dir", required=True, help="Results directory with ROOT files.")
    parser.add_argument(
        "--output-dir",
        default="plots/MC_electrons_on_shell_dionne-model/diagnostics_edep_depth/per_energy_counts",
        help="Output directory for plots.",
    )
    parser.add_argument("--label", default=None, help="Label for output filenames.")
    args = parser.parse_args()

    import ROOT
    ROOT.TH1.AddDirectory(False)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    def _meta_string(value):
        try:
            raw = bytes(value)
            if raw:
                return raw.decode(errors="ignore").rstrip("\x00").strip()
        except Exception:
            pass
        return str(value).strip()

    root_files = [
        f for f in os.listdir(args.results_dir)
        if f.endswith(".root") and "_SEY_MonteCarlo" not in f and not f.startswith("summary")
    ]
    root_files.sort()

    if not root_files:
        raise SystemExit(f"No ROOT files found in {args.results_dir}")

    os.makedirs(args.output_dir, exist_ok=True)
    label = args.label or os.path.basename(os.path.normpath(args.results_dir))

    curves = []
    thickness_nm = None
    em_model = None
    n_events = None

    for name in root_files:
        path = os.path.join(args.results_dir, name)
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
        h_counts = f.Get("EdepDepthPrimaryCounts")
        if not h_counts:
            f.Close()
            raise SystemExit(
                f"Missing EdepDepthPrimaryCounts in {path}. Re-run the Geant4 scan."
            )
        h = h_counts.Clone()
        h.SetDirectory(0)
        # Histograms are stored per-event in results.
        curves.append((energy_ev, h))
        f.Close()

    if not curves:
        raise SystemExit("No curves loaded.")

    # Sort by energy
    curves.sort(key=lambda x: x[0])

    colors = [
        ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1,
        ROOT.kOrange + 7, ROOT.kCyan + 1, ROOT.kViolet + 1, ROOT.kAzure + 2,
        ROOT.kTeal + 1, ROOT.kPink + 2, ROOT.kGray + 2
    ]

    c = ROOT.TCanvas("c_overlay", "Depth counts overlay", 1000, 700)
    c.SetRightMargin(0.3)

    max_y = 0.0
    for _, h in curves:
        max_y = max(max_y, h.GetMaximum())

    first = True
    leg = ROOT.TLegend(0.72, 0.12, 0.98, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(1)

    for idx, (energy_ev, h) in enumerate(curves):
        color = colors[idx % len(colors)]
        h.SetLineWidth(2)
        h.SetLineColor(color)
        if first:
            h.SetTitle("")
            h.SetXTitle("Depth from entrance (nm)")
            h.SetYTitle("Energy-depositing steps per event")
            h.SetMaximum(max_y * 1.15)
            h.Draw("HIST")
            first = False
        else:
            h.Draw("HIST SAME")
        leg.AddEntry(h, f"{_format_energy_ev(energy_ev)} eV", "l")

    leg.Draw()

    text = ROOT.TPaveText(0.12, 0.82, 0.5, 0.95, "NDC")
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

    out_base = os.path.join(args.output_dir, f"edep_depth_counts_overlay_{label}")
    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".png")
    c.Close()

    print(f"Saved overlay plot: {out_base}.pdf/.png")


if __name__ == "__main__":
    main()
