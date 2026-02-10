#!/usr/bin/env python3
import argparse
import os

def _format_energy_ev(value):
    if abs(value - round(value)) < 1e-6:
        return f"{int(round(value))}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def main():
    parser = argparse.ArgumentParser(description="Plot EdepDepthPrimaryCounts per primary energy.")
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
    ROOT.gStyle.SetOptStat(0)

    root_files = [
        f for f in os.listdir(args.results_dir)
        if f.endswith(".root") and "_SEY_MonteCarlo" not in f and not f.startswith("summary")
    ]
    root_files.sort()

    if not root_files:
        raise SystemExit(f"No ROOT files found in {args.results_dir}")

    os.makedirs(args.output_dir, exist_ok=True)
    label = args.label or os.path.basename(os.path.normpath(args.results_dir))

    for idx, name in enumerate(root_files):
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
        if n_events > 0:
            h.Scale(1.0 / n_events)

        c = ROOT.TCanvas(f"c_{idx}", "Edep depth counts", 900, 650)
        h.SetTitle("")
        h.SetLineWidth(2)
        h.SetLineColor(ROOT.kBlue + 1)
        h.SetXTitle("Depth from entrance (nm)")
        h.SetYTitle("Energy-depositing steps per event")
        h.Draw("HIST")

        text = ROOT.TPaveText(0.12, 0.82, 0.5, 0.95, "NDC")
        text.SetFillColor(0)
        text.SetBorderSize(1)
        text.SetTextFont(42)
        text.SetTextSize(0.028)
        text.AddText(f"Energy: {_format_energy_ev(energy_ev)} eV")
        text.AddText(f"Thickness: {thickness_nm:.0f} nm")
        if em_model:
            text.AddText(f"EM model: {em_model}")
        text.AddText(f"MC events: {n_events}")
        text.Draw()

        e_tag = _format_energy_ev(energy_ev)
        out_base = os.path.join(
            args.output_dir,
            f"edep_depth_counts_e{e_tag}eV_{label}"
        )
        c.SaveAs(out_base + ".pdf")
        c.SaveAs(out_base + ".png")
        c.Close()
        f.Close()

    print(f"Saved per-energy count plots in {args.output_dir}")


if __name__ == "__main__":
    main()
