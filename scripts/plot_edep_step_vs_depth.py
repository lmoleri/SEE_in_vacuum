#!/usr/bin/env python3
import argparse
import os
def _format_energy_ev(value):
    if abs(value - round(value)) < 1e-6:
        return f"{int(round(value))}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def _meta_string(value):
    try:
        raw = bytes(value)
        if raw:
            return raw.decode(errors="ignore").rstrip("\x00").strip()
    except Exception:
        pass
    return str(value).strip()


def _energy_edges(energies):
    if len(energies) == 1:
        e = energies[0]
        left = max(0.0, e * 0.5)
        return [left, e * 1.5]
    edges = [energies[0] - (energies[1] - energies[0]) / 2.0]
    for i in range(1, len(energies)):
        edges.append((energies[i - 1] + energies[i]) / 2.0)
    edges.append(energies[-1] + (energies[-1] - energies[-2]) / 2.0)
    if edges[0] < 0.0:
        edges[0] = 0.0
    return edges


def _max_nonzero_y(h2):
    ny = h2.GetNbinsY()
    nx = h2.GetNbinsX()
    for iy in range(ny, 0, -1):
        for ix in range(1, nx + 1):
            if h2.GetBinContent(ix, iy) > 0:
                return h2.GetYaxis().GetBinUpEdge(iy)
    return None


def main():
    parser = argparse.ArgumentParser(description="Plot EdepStepVsDepthPrimary per energy and stacked.")
    parser.add_argument("--results-dir", required=True, help="Results directory with ROOT files.")
    parser.add_argument(
        "--output-dir",
        default="plots/MC_electrons_on_shell_dionne-model/diagnostics_edep_depth/step_vs_depth",
        help="Output directory for plots.",
    )
    parser.add_argument("--label", default=None, help="Label for output filenames.")
    parser.add_argument("--per-energy", action="store_true", help="Produce per-energy plots.")
    parser.add_argument(
        "--save-results",
        action="store_true",
        help="Also save plots under <results-dir>/plots.",
    )
    args = parser.parse_args()

    if not args.per_energy:
        args.per_energy = True

    import ROOT
    ROOT.TH1.AddDirectory(False)
    ROOT.gROOT.SetBatch(True)
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

    records = []
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
        max_step_nm = float(getattr(meta, "maxStepNm", 0.0))
        em_model = _meta_string(meta.emModel)
        h2 = f.Get("EdepStepVsDepthPrimary")
        if not h2:
            f.Close()
            raise SystemExit(
                f"Missing EdepStepVsDepthPrimary in {path}. Re-run the Geant4 scan."
            )
        h2c = h2.Clone()
        h2c.SetDirectory(0)
        records.append({
            "energy_eV": energy_ev,
            "n_events": n_events,
            "thickness_nm": thickness_nm,
            "max_step_nm": max_step_nm,
            "em_model": em_model,
            "h2": h2c,
        })
        f.Close()

    records.sort(key=lambda r: r["energy_eV"])
    if not records:
        raise SystemExit("No records loaded.")

    thickness_nm = records[0]["thickness_nm"]
    em_model = records[0]["em_model"]
    n_events = records[0]["n_events"]

    if args.per_energy:
        for rec in records:
            h2 = rec["h2"].Clone()
            h2.SetDirectory(0)
            if rec["n_events"] > 0:
                h2.Scale(1.0 / rec["n_events"])
            if rec["max_step_nm"] > 0:
                bin_width = h2.GetXaxis().GetBinWidth(1)
                if bin_width > 0:
                    factor = int(round(rec["max_step_nm"] / bin_width))
                    if factor > 1:
                        h2.RebinX(factor)

            c = ROOT.TCanvas("c", "Step edep vs depth", 1000, 700)
            c.SetRightMargin(0.16)
            h2.SetTitle("")
            h2.GetXaxis().SetTitle("Depth from entrance (nm)")
            h2.GetYaxis().SetTitle("Energy deposition per step (eV)")
            y_max = _max_nonzero_y(h2)
            if y_max is None:
                y_max = rec["energy_eV"]
            else:
                # Zoom to the occupied range; keep a minimum span so the band is visible.
                y_max = max(5.0, min(rec["energy_eV"], y_max * 1.2))
            h2.GetYaxis().SetRangeUser(0.0, y_max)
            h2.GetZaxis().SetTitle("Steps per event")
            h2.Draw("COLZ")

            text = ROOT.TPaveText(0.12, 0.82, 0.5, 0.95, "NDC")
            text.SetFillColor(0)
            text.SetBorderSize(1)
            text.SetTextFont(42)
            text.SetTextSize(0.028)
            text.AddText(f"Energy: {_format_energy_ev(rec['energy_eV'])} eV")
            text.AddText(f"Thickness: {thickness_nm:.0f} nm")
            if em_model:
                text.AddText(f"EM model: {em_model}")
            text.AddText(f"MC events: {rec['n_events']}")
            text.Draw()

            e_tag = _format_energy_ev(rec["energy_eV"])
            out_base = os.path.join(args.output_dir, f"edep_step_vs_depth_e{e_tag}eV_{label}")
            c.SaveAs(out_base + ".pdf")
            c.SaveAs(out_base + ".png")
            if args.save_results:
                results_plot_dir = os.path.join(args.results_dir, "plots")
                os.makedirs(results_plot_dir, exist_ok=True)
                out_base_results = os.path.join(
                    results_plot_dir, f"edep_step_vs_depth_e{e_tag}eV_{label}"
                )
                c.SaveAs(out_base_results + ".pdf")
                c.SaveAs(out_base_results + ".png")
            c.Close()

    print(f"Saved plots to {args.output_dir}")


if __name__ == "__main__":
    main()
