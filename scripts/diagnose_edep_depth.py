#!/usr/bin/env python3
import argparse
import os
from array import array


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


def _load_records(results_dir):
    import ROOT
    ROOT.TH1.AddDirectory(False)

    root_files = [
        f for f in os.listdir(results_dir)
        if f.endswith(".root") and "_SEY_MonteCarlo" not in f and not f.startswith("summary")
    ]
    root_files.sort()
    records = []
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
        energy_eV = float(meta.primaryEnergyMeV) * 1.0e6
        n_events = int(meta.primaryElectrons)
        thickness_nm = float(meta.sampleThicknessNm)
        em_model = str(meta.emModel).strip()
        h_depth = f.Get("EdepDepthPrimary")
        h_depth_w = f.Get("EdepDepthPrimaryWeighted")
        h_counts = f.Get("EdepDepthPrimaryCounts")
        if not h_depth or not h_depth_w:
            f.Close()
            raise SystemExit(
                "Missing EdepDepthPrimary/EdepDepthPrimaryWeighted in {}. "
                "Rebuild and re-run the Geant4 scan.".format(path)
            )
        h_depth_clone = h_depth.Clone()
        h_depth_clone.SetDirectory(0)
        h_depth_w_clone = h_depth_w.Clone()
        h_depth_w_clone.SetDirectory(0)
        h_counts_clone = None
        if h_counts:
            h_counts_clone = h_counts.Clone()
            h_counts_clone.SetDirectory(0)
        records.append({
            "energy_eV": energy_eV,
            "n_events": n_events,
            "thickness_nm": thickness_nm,
            "em_model": em_model,
            "h_depth": h_depth_clone,
            "h_depth_w": h_depth_w_clone,
            "h_counts": h_counts_clone,
        })
        f.Close()
    records.sort(key=lambda r: r["energy_eV"])
    return records


def main():
    parser = argparse.ArgumentParser(description="Diagnose depth-resolved energy deposition vs energy.")
    parser.add_argument("--results-dir", required=True, help="Results directory containing ROOT files.")
    parser.add_argument("--output-dir", default="plots/MC_electrons_on_shell_dionne-model/diagnostics_edep_depth",
                        help="Output directory for plots.")
    parser.add_argument("--label", default=None, help="Label for output filenames.")
    parser.add_argument("--mode", choices=["edep", "counts"], default="edep",
                        help="Plot mode: energy deposition (edep) or step counts (counts).")
    args = parser.parse_args()

    import ROOT

    ROOT.gStyle.SetOptStat(0)
    records = _load_records(args.results_dir)
    if not records:
        raise SystemExit("No ROOT files found in {}".format(args.results_dir))

    energies = [r["energy_eV"] for r in records]
    n_events = records[0]["n_events"]
    thickness_nm = records[0]["thickness_nm"]
    em_model = records[0]["em_model"]

    use_counts = args.mode == "counts"
    h_ref = records[0]["h_counts"] if use_counts else records[0]["h_depth"]
    if use_counts and h_ref is None:
        raise SystemExit("Counts mode requires EdepDepthPrimaryCounts. Re-run the Geant4 scan.")
    n_depth_bins = h_ref.GetNbinsX()
    depth_edges = [h_ref.GetXaxis().GetBinLowEdge(1)]
    for i in range(1, n_depth_bins + 1):
        depth_edges.append(h_ref.GetXaxis().GetBinUpEdge(i))

    energy_edges = _energy_edges(energies)

    x_edges = array("d", energy_edges)
    y_edges = array("d", depth_edges)

    h2 = ROOT.TH2D(
        "h_edep_depth",
        "",
        len(energies),
        x_edges,
        n_depth_bins,
        y_edges,
    )
    h2w = ROOT.TH2D(
        "h_edep_depth_weighted",
        "",
        len(energies),
        x_edges,
        n_depth_bins,
        y_edges,
    )

    total_edep = []
    total_edep_w = []

    for i, rec in enumerate(records, start=1):
        if use_counts:
            h = rec["h_counts"]
            hw = None
        else:
            h = rec["h_depth"]
            hw = rec["h_depth_w"]
        events = max(1, rec["n_events"])
        scale = 1.0 / events
        s = 0.0
        sw = 0.0
        for b in range(1, n_depth_bins + 1):
            val = h.GetBinContent(b) * scale
            valw = hw.GetBinContent(b) * scale if hw else 0.0
            h2.SetBinContent(i, b, val)
            if hw:
                h2w.SetBinContent(i, b, valw)
            s += val
            sw += valw
        total_edep.append(s)
        total_edep_w.append(sw)

    os.makedirs(args.output_dir, exist_ok=True)
    label = args.label or os.path.basename(os.path.normpath(args.results_dir))

    def _add_text(pad):
        text = ROOT.TPaveText(0.12, 0.83, 0.5, 0.95, "NDC")
        text.SetFillColor(0)
        text.SetBorderSize(1)
        text.SetTextFont(42)
        text.SetTextSize(0.028)
        text.AddText(f"Thickness: {thickness_nm:.0f} nm")
        if em_model:
            text.AddText(f"EM model: {em_model}")
        text.AddText(f"MC events: {n_events}")
        text.Draw()

    if args.mode == "edep":
        # Unweighted heatmap
        c1 = ROOT.TCanvas("c_edep_depth", "Edep depth vs energy", 1000, 700)
        c1.SetRightMargin(0.16)
        h2.SetTitle("")
        h2.SetXTitle("Primary electron energy (eV)")
        h2.SetYTitle("Depth from entrance (nm)")
        h2.SetZTitle("Mean energy deposition per event (eV)")
        h2.Draw("COLZ")
        _add_text(c1)
        out1 = os.path.join(args.output_dir, f"edep_depth_vs_energy_{label}.pdf")
        c1.SaveAs(out1)
        c1.SaveAs(out1.replace(".pdf", ".png"))

        # Weighted heatmap
        c2 = ROOT.TCanvas("c_edep_depth_weighted", "Weighted edep depth vs energy", 1000, 700)
        c2.SetRightMargin(0.16)
        h2w.SetTitle("")
        h2w.SetXTitle("Primary electron energy (eV)")
        h2w.SetYTitle("Depth from entrance (nm)")
        h2w.SetZTitle("Mean weighted energy deposition per event (eV)")
        h2w.Draw("COLZ")
        _add_text(c2)
        out2 = os.path.join(args.output_dir, f"edep_depth_weighted_vs_energy_{label}.pdf")
        c2.SaveAs(out2)
        c2.SaveAs(out2.replace(".pdf", ".png"))

        # Integrals vs energy
        gr = ROOT.TGraph(len(energies))
        grw = ROOT.TGraph(len(energies))
        for i, (e, s, sw) in enumerate(zip(energies, total_edep, total_edep_w)):
            gr.SetPoint(i, e, s)
            grw.SetPoint(i, e, sw)
        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(1.2)
        gr.SetLineWidth(2)
        gr.SetLineColor(ROOT.kBlue + 1)
        gr.SetMarkerColor(ROOT.kBlue + 1)
        grw.SetMarkerStyle(21)
        grw.SetMarkerSize(1.2)
        grw.SetLineWidth(2)
        grw.SetLineColor(ROOT.kRed + 1)
        grw.SetMarkerColor(ROOT.kRed + 1)

        c3 = ROOT.TCanvas("c_edep_integrals", "Edep integrals vs energy", 900, 650)
        frame = c3.DrawFrame(min(energies) * 0.9, 0.0, max(energies) * 1.05, max(max(total_edep), max(total_edep_w)) * 1.2)
        frame.SetTitle("")
        frame.SetXTitle("Primary electron energy (eV)")
        frame.SetYTitle("Mean energy deposition per event (eV)")
        gr.Draw("LP")
        grw.Draw("LP")
        leg = ROOT.TLegend(0.58, 0.72, 0.88, 0.86)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(gr, "Total Edep", "lp")
        leg.AddEntry(grw, "Weighted Edep", "lp")
        leg.Draw()
        _add_text(c3)
        out3 = os.path.join(args.output_dir, f"edep_integrals_vs_energy_{label}.pdf")
        c3.SaveAs(out3)
        c3.SaveAs(out3.replace(".pdf", ".png"))

        print("Saved:")
        print(out1)
        print(out2)
        print(out3)
    else:
        # Counts heatmap
        c1 = ROOT.TCanvas("c_edep_depth_counts", "Step counts vs depth", 1000, 700)
        c1.SetRightMargin(0.16)
        h2.SetTitle("")
        h2.SetXTitle("Primary electron energy (eV)")
        h2.SetYTitle("Depth from entrance (nm)")
        h2.SetZTitle("Mean energy-depositing steps per event")
        h2.Draw("COLZ")
        _add_text(c1)
        out1 = os.path.join(args.output_dir, f"edep_depth_counts_vs_energy_{label}.pdf")
        c1.SaveAs(out1)
        c1.SaveAs(out1.replace(".pdf", ".png"))
        print("Saved:")
        print(out1)


if __name__ == "__main__":
    main()
