#!/usr/bin/env python3
import argparse
import os


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
        h_edep = f.Get("EdepPrimary")
        h_len = f.Get("PrimaryTrackLengthAl2O3")
        if not h_edep or not h_len:
            f.Close()
            raise SystemExit(
                "Missing EdepPrimary/PrimaryTrackLengthAl2O3 in {}. "
                "Rebuild and re-run the Geant4 scan.".format(path)
            )
        h_edep_clone = h_edep.Clone()
        h_edep_clone.SetDirectory(0)
        h_len_clone = h_len.Clone()
        h_len_clone.SetDirectory(0)
        records.append({
            "energy_eV": energy_eV,
            "n_events": n_events,
            "thickness_nm": thickness_nm,
            "em_model": em_model,
            "h_edep": h_edep_clone,
            "h_len": h_len_clone,
        })
        f.Close()
    records.sort(key=lambda r: r["energy_eV"])
    return records


def main():
    parser = argparse.ArgumentParser(description="Diagnose mean dE/dx vs energy from Geant4 output.")
    parser.add_argument("--results-dir", required=True, help="Results directory containing ROOT files.")
    parser.add_argument("--output-dir", default="plots/MC_electrons_on_shell_dionne-model/diagnostics_dedx",
                        help="Output directory for plots.")
    parser.add_argument("--label", default=None, help="Label for output filenames.")
    args = parser.parse_args()

    import ROOT
    ROOT.gStyle.SetOptStat(0)

    records = _load_records(args.results_dir)
    if not records:
        raise SystemExit("No ROOT files found in {}".format(args.results_dir))

    energies = []
    dedx_vals = []
    length_vals = []
    n_events = records[0]["n_events"]
    thickness_nm = records[0]["thickness_nm"]
    em_model = records[0]["em_model"]

    for rec in records:
        mean_edep = rec["h_edep"].GetMean()  # eV
        mean_len = rec["h_len"].GetMean()    # nm
        if mean_len <= 0:
            continue
        energies.append(rec["energy_eV"])
        dedx_vals.append(mean_edep / mean_len)
        length_vals.append(mean_len)

    if not energies:
        raise SystemExit("No valid dE/dx points (mean track length was zero).")

    os.makedirs(args.output_dir, exist_ok=True)
    label = args.label or os.path.basename(os.path.normpath(args.results_dir))

    gr = ROOT.TGraph(len(energies))
    for i, (e, v) in enumerate(zip(energies, dedx_vals)):
        gr.SetPoint(i, e, v)
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(1.2)
    gr.SetLineWidth(2)
    gr.SetLineColor(ROOT.kBlue + 1)
    gr.SetMarkerColor(ROOT.kBlue + 1)

    c1 = ROOT.TCanvas("c_dedx", "dE/dx vs energy", 900, 650)
    frame = c1.DrawFrame(min(energies) * 0.9, 0.0, max(energies) * 1.05, max(dedx_vals) * 1.2)
    frame.SetTitle("")
    frame.SetXTitle("Primary electron energy (eV)")
    frame.SetYTitle("Mean dE/dx (eV/nm)")
    gr.Draw("LP")

    text = ROOT.TPaveText(0.12, 0.80, 0.55, 0.94, "NDC")
    text.SetFillColor(0)
    text.SetBorderSize(1)
    text.SetTextFont(42)
    text.SetTextSize(0.028)
    text.AddText(f"Thickness: {thickness_nm:.0f} nm")
    if em_model:
        text.AddText(f"EM model: {em_model}")
    text.AddText(f"MC events: {n_events}")
    text.AddText("dE/dx ≈ <Edep>/<L> (mean per event)")
    text.Draw()

    out = os.path.join(args.output_dir, f"dedx_vs_energy_{label}.pdf")
    c1.SaveAs(out)

    # Print peak for quick diagnosis
    peak_idx = max(range(len(dedx_vals)), key=lambda i: dedx_vals[i])
    print(f"Peak dE/dx at {energies[peak_idx]:.0f} eV (value={dedx_vals[peak_idx]:.4g} eV/nm)")
    print("Saved:")
    print(out)


if __name__ == "__main__":
    main()
