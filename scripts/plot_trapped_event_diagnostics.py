#!/usr/bin/env python3
import argparse
import os

import ROOT


def _find_root_files(results_dir):
    files = []
    for root, _, names in os.walk(results_dir):
        for name in names:
            if not name.endswith(".root"):
                continue
            if "_SEY_MonteCarlo" in name or name.startswith("summary"):
                continue
            files.append(os.path.join(root, name))
    return sorted(files)


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


def _load_runs(results_dir):
    runs = []
    for path in _find_root_files(results_dir):
        f = ROOT.TFile.Open(path)
        if not f or not f.IsOpen():
            continue
        meta = f.Get("RunMeta")
        tree = f.Get("EventDiagnostics")
        if not meta or not tree:
            f.Close()
            continue
        meta.GetEntry(0)
        runs.append(
            {
                "path": path,
                "file": f,
                "tree": tree,
                "step_nm": float(meta.maxStepNm),
                "e0_eV": float(meta.primaryEnergyMeV) * 1.0e6,
                "events": int(meta.primaryElectrons),
                "model": _meta_string(meta.emModel),
                "thickness_nm": float(meta.sampleThicknessNm),
            }
        )
    runs.sort(key=lambda r: r["step_nm"])
    return runs


def _draw_overlay_pad(pad, runs, var, title, x_title, bins, x_min, x_max, use_logy):
    pad.cd()
    pad.SetLeftMargin(0.12)
    pad.SetRightMargin(0.03)
    pad.SetBottomMargin(0.14)
    pad.SetTopMargin(0.08)
    pad.SetLogy(use_logy)

    colors = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7]
    hists = []
    y_max = 0.0

    leg = ROOT.TLegend(0.46, 0.58, 0.94, 0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)

    sel = "primaryExitClass==4"
    for i, run in enumerate(runs):
        h = ROOT.TH1D(
            f"h_{var}_{i}_{ROOT.gRandom.Integer(10_000_000)}",
            "",
            bins,
            x_min,
            x_max,
        )
        n = run["tree"].Draw(f"{var}>>{h.GetName()}", sel, "goff")
        h.SetDirectory(0)
        if n <= 0:
            hists.append(h)
            continue
        y_max = max(y_max, h.GetMaximum())
        color = colors[i % len(colors)]
        h.SetLineColor(color)
        h.SetLineWidth(3)
        hists.append(h)
        leg.AddEntry(h, f"step={run['step_nm']:.3g} nm (N={n})", "l")

    if y_max <= 0:
        frame = ROOT.TH2D(
            f"frame_{var}_{ROOT.gRandom.Integer(10_000_000)}",
            "",
            10,
            x_min,
            x_max,
            10,
            0.0,
            1.0,
        )
        frame.SetTitle(f"{title};{x_title};Entries")
        frame.Draw("AXIS")
        txt = ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(42)
        txt.SetTextSize(0.045)
        txt.DrawLatex(0.35, 0.50, "No entries")
        return [frame, txt, leg] + hists

    y_min = 1e-1 if use_logy else 0.0
    y_top = max(1.3 * y_max, 1.0)
    if use_logy:
        y_top = max(y_top, 10.0)

    frame = ROOT.TH2D(
        f"frame_{var}_{ROOT.gRandom.Integer(10_000_000)}",
        "",
        10,
        x_min,
        x_max,
        10,
        y_min if use_logy else 0.0,
        y_top,
    )
    frame.SetTitle(f"{title};{x_title};Entries")
    frame.Draw("AXIS")

    drawn = False
    for h in hists:
        if h.GetEntries() > 0:
            if use_logy:
                h.SetMinimum(y_min)
            if not drawn:
                h.Draw("HIST SAME")
                drawn = True
            else:
                h.Draw("HIST SAME")
    leg.Draw()
    return [frame, leg] + hists


def main():
    parser = argparse.ArgumentParser(
        description="Plot trapped-event (primaryExitClass==4) diagnostics vs max step."
    )
    parser.add_argument("--results-dir", required=True, help="Directory with ROOT outputs.")
    parser.add_argument("--output-dir", required=True, help="Directory for output plots.")
    parser.add_argument("--label", default="trapped", help="Output filename suffix.")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    runs = _load_runs(args.results_dir)
    if not runs:
        raise SystemExit("No valid ROOT files with EventDiagnostics found.")

    e0 = runs[0]["e0_eV"]
    model = runs[0]["model"]
    events = runs[0]["events"]
    thickness = runs[0]["thickness_nm"]

    os.makedirs(args.output_dir, exist_ok=True)

    variables = [
        ("nDirectionReversalsZ", "z-direction reversals (trapped)", "Reversal count", 30, 0.0, 30.0, True),
        ("nBoundaryCrossings", "Boundary crossings (trapped)", "Crossing count", 8, -0.5, 7.5, False),
        ("nEdepSteps", "Energy-depositing steps (trapped)", "Step count", 150, 0.0, 300.0, True),
        ("primaryTrackLengthNm", "Track length in Al2O3 (trapped)", "Track length (nm)", 120, 0.0, 30.0, True),
        ("maxDepthNm", "Max depth reached (trapped)", "Depth from entrance (nm)", 100, 0.0, max(5.0, thickness), True),
        ("primaryResidualEv", "Residual energy at end (trapped)", "Residual energy (eV)", 120, 0.0, max(1.0, 1.02 * e0), True),
        ("edepPrimaryEv", "Total deposited energy (trapped)", "EdepPrimary (eV)", 120, 0.0, max(1.0, 1.02 * e0), True),
        ("edepByEIoniEv", "Deposited by eIoni (trapped)", "Energy (eV)", 120, 0.0, max(50.0, 0.2 * e0), True),
        ("edepByMscEv", "Deposited by msc (trapped)", "Energy (eV)", 120, 0.0, max(1.0, 1.02 * e0), True),
        ("edepByOtherEv", "Deposited by other processes (trapped)", "Energy (eV)", 120, 0.0, max(1.0, 1.02 * e0), True),
        ("firstDirectionReversalStep", "First reversal step (trapped)", "Step index", 140, 0.0, 280.0, True),
        ("firstDirectionReversalDepthNm", "First reversal depth (trapped)", "Depth from entrance (nm)", 100, 0.0, max(5.0, thickness), True),
        ("firstDirectionReversalEnergyEv", "First reversal energy (trapped)", "Kinetic energy (eV)", 120, 0.0, max(1.0, 1.02 * e0), True),
    ]

    c = ROOT.TCanvas("c_trapped_diag", "Trapped diagnostics", 2100, 1800)
    c.Divide(4, 4)
    keep = [c]

    for i, (var, title, x_title, bins, x_min, x_max, logy) in enumerate(variables, start=1):
        keep.extend(_draw_overlay_pad(c.cd(i), runs, var, title, x_title, bins, x_min, x_max, logy))

    c.cd(16)
    ROOT.gPad.SetLeftMargin(0.08)
    ROOT.gPad.SetRightMargin(0.02)
    ROOT.gPad.SetBottomMargin(0.08)
    ROOT.gPad.SetTopMargin(0.08)
    info = ROOT.TPaveText(0.05, 0.05, 0.95, 0.95, "NDC")
    info.SetFillColor(0)
    info.SetBorderSize(1)
    info.SetTextFont(42)
    info.SetTextSize(0.045)
    info.AddText("Trapped-event diagnostics")
    info.AddText("Selection: primaryExitClass == 4")
    info.AddText(f"E0 = {e0:.0f} eV")
    info.AddText(f"thickness = {thickness:.1f} nm")
    info.AddText(f"EM model = {model}")
    info.AddText(f"MC events per run = {events}")
    info.Draw()
    keep.append(info)

    title_box = ROOT.TPaveText(0.08, 0.965, 0.92, 0.995, "NDC")
    title_box.SetFillStyle(0)
    title_box.SetBorderSize(0)
    title_box.SetTextFont(42)
    title_box.SetTextSize(0.03)
    title_box.AddText("Class-4 (trapped) diagnostics vs max step")
    title_box.Draw()
    keep.append(title_box)

    out_base = os.path.join(args.output_dir, f"trapped_event_diagnostics_{args.label}")
    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".png")

    for run in runs:
        run["file"].Close()

    print("Saved:")
    print(out_base + ".pdf")
    print(out_base + ".png")


if __name__ == "__main__":
    main()
