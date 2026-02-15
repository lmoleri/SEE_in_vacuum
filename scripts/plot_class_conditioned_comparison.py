#!/usr/bin/env python3
import argparse
import os


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


def _find_root_files(results_dir):
    files = []
    for root, _, names in os.walk(results_dir):
        for name in names:
            if not name.endswith(".root"):
                continue
            if "_SEY_MonteCarlo" in name or name.startswith("summary"):
                continue
            files.append(os.path.join(root, name))
    files.sort()
    return files


def _load_runs(results_dir):
    import ROOT

    ROOT.TH1.AddDirectory(False)

    hnames = [
        "EdepPrimaryStop",
        "EdepPrimaryExitEntrance",
        "EdepPrimaryExitOpposite",
        "EdepPrimaryExitLateral",
    ]

    runs = []
    for path in _find_root_files(results_dir):
        f = ROOT.TFile.Open(path)
        if not f or not f.IsOpen():
            continue

        meta = f.Get("RunMeta")
        if not meta:
            f.Close()
            continue
        meta.GetEntry(0)

        step_nm = float(meta.maxStepNm)
        energy_eV = float(meta.primaryEnergyMeV) * 1.0e6
        n_events = int(meta.primaryElectrons)
        em_model = _meta_string(meta.emModel)

        hists = {}
        missing = False
        for hname in hnames:
            h = f.Get(hname)
            if not h:
                missing = True
                break
            hc = h.Clone(f"{hname}_step{step_nm}")
            hc.SetDirectory(0)
            hists[hname] = hc

        if missing:
            f.Close()
            continue

        runs.append(
            {
                "path": path,
                "step_nm": step_nm,
                "energy_eV": energy_eV,
                "n_events": n_events,
                "em_model": em_model,
                "hists": hists,
            }
        )
        f.Close()

    runs.sort(key=lambda r: r["step_nm"])
    return runs


def _normalize(h):
    integral = h.Integral()
    if integral > 0:
        h.Scale(1.0 / integral)


def main():
    parser = argparse.ArgumentParser(
        description="Class-conditioned EdepPrimary comparison plots from EventDiagnostics outputs."
    )
    parser.add_argument("--results-dir", required=True, help="Directory with ROOT outputs.")
    parser.add_argument("--output-dir", required=True, help="Output directory for plots.")
    parser.add_argument("--label", default="class_conditioned", help="Output label suffix.")
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize class histograms to unit area before overlay (recommended).",
    )
    args = parser.parse_args()

    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    runs = _load_runs(args.results_dir)
    if not runs:
        raise SystemExit("No usable ROOT files found with class-conditioned histograms.")

    energy_eV = runs[0]["energy_eV"]
    n_events = runs[0]["n_events"]
    em_model = runs[0]["em_model"]

    for r in runs[1:]:
        if abs(r["energy_eV"] - energy_eV) > 1e-6:
            raise SystemExit("Mixed energies found. Use one energy campaign per call.")

    os.makedirs(args.output_dir, exist_ok=True)

    colors = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1]
    class_defs = [
        ("EdepPrimaryStop", "Stop / No Valid Exit"),
        ("EdepPrimaryExitEntrance", "Entrance Exit"),
        ("EdepPrimaryExitOpposite", "Opposite Exit"),
        ("EdepPrimaryExitLateral", "Lateral Exit"),
    ]

    keep = []

    # Plot 1: 2x2 class-conditioned overlays.
    c1 = ROOT.TCanvas("c_class_overlay", "Class overlays", 1500, 1100)
    keep.append(c1)
    c1.Divide(2, 2)

    x_max = max(energy_eV * 1.02, 1.0)

    for ipad, (hname, title) in enumerate(class_defs, start=1):
        c1.cd(ipad)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.03)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTopMargin(0.08)

        clones = []
        y_max = 0.0
        for i, run in enumerate(runs):
            h = run["hists"][hname].Clone(f"{hname}_clone_{i}")
            h.SetDirectory(0)
            if args.normalize:
                _normalize(h)
            y_max = max(y_max, h.GetMaximum())
            clones.append((run, h))
            keep.append(h)

        if y_max <= 0:
            y_max = 1.0

        frame = ROOT.TH2D(
            f"frame_{hname}",
            "",
            10,
            0.0,
            x_max,
            10,
            0.0,
            y_max * 1.2,
        )
        keep.append(frame)
        frame.SetTitle(f"{title};EdepPrimary (eV);{'Arbitrary units (area=1)' if args.normalize else 'Events'}")
        frame.Draw("AXIS")

        leg = ROOT.TLegend(0.56, 0.60, 0.94, 0.90)
        keep.append(leg)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.034)

        for i, (run, h) in enumerate(clones):
            color = colors[i % len(colors)]
            h.SetLineColor(color)
            h.SetLineWidth(3)
            h.Draw("HIST SAME")
            frac = run["hists"][hname].Integral() / float(run["n_events"]) if run["n_events"] > 0 else 0.0
            leg.AddEntry(h, f"step={run['step_nm']:.3g} nm  (f={frac:.3f})", "l")
        leg.Draw()

    out1 = os.path.join(args.output_dir, f"class_conditioned_edep_overlays_{args.label}")
    c1.SaveAs(out1 + ".pdf")
    c1.SaveAs(out1 + ".png")

    # Plot 2: class fractions vs step.
    c2 = ROOT.TCanvas("c_fractions", "Class fractions vs step", 1000, 750)
    keep.append(c2)
    c2.SetLeftMargin(0.12)
    c2.SetRightMargin(0.04)
    c2.SetBottomMargin(0.12)
    c2.SetTopMargin(0.08)

    step_min = min(r["step_nm"] for r in runs)
    step_max = max(r["step_nm"] for r in runs)
    frame2 = ROOT.TH2D(
        "frame_fractions",
        "",
        10,
        step_min * 0.9,
        step_max * 1.1,
        10,
        0.0,
        1.05,
    )
    keep.append(frame2)
    frame2.SetTitle("Primary event class fractions vs max step;max step in Al_{2}O_{3} (nm);Fraction of events")
    frame2.Draw("AXIS")

    frac_specs = [
        ("EdepPrimaryStop", ROOT.kRed + 1, 20, "Stop / no-valid-exit"),
        ("EdepPrimaryExitEntrance", ROOT.kOrange + 7, 21, "Entrance exit"),
        ("EdepPrimaryExitOpposite", ROOT.kBlue + 1, 22, "Opposite exit"),
        ("EdepPrimaryExitLateral", ROOT.kGreen + 2, 33, "Lateral exit"),
    ]

    leg2 = ROOT.TLegend(0.56, 0.64, 0.93, 0.90)
    keep.append(leg2)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.03)

    for hname, color, marker, label in frac_specs:
        g = ROOT.TGraph(len(runs))
        keep.append(g)
        for i, run in enumerate(runs):
            frac = run["hists"][hname].Integral() / float(run["n_events"]) if run["n_events"] > 0 else 0.0
            g.SetPoint(i, run["step_nm"], frac)
        g.SetLineColor(color)
        g.SetMarkerColor(color)
        g.SetMarkerStyle(marker)
        g.SetMarkerSize(1.3)
        g.SetLineWidth(2)
        g.Draw("LP SAME")
        leg2.AddEntry(g, label, "lp")
    leg2.Draw()

    info = ROOT.TPaveText(0.14, 0.82, 0.46, 0.92, "NDC")
    keep.append(info)
    info.SetFillColor(0)
    info.SetBorderSize(1)
    info.SetTextFont(42)
    info.SetTextSize(0.028)
    info.AddText(f"E0 = {energy_eV:.0f} eV")
    info.AddText(f"EM model: {em_model}")
    info.AddText(f"MC events: {n_events}")
    info.Draw()

    out2 = os.path.join(args.output_dir, f"class_fractions_vs_step_{args.label}")
    c2.SaveAs(out2 + ".pdf")
    c2.SaveAs(out2 + ".png")

    print("Saved:")
    print(out1 + ".pdf")
    print(out1 + ".png")
    print(out2 + ".pdf")
    print(out2 + ".png")


if __name__ == "__main__":
    main()
