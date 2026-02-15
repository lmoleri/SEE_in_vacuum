#!/usr/bin/env python3
import argparse
import os

import ROOT


CLASS_DEFS = [
    (1, "Entrance exit"),
    (2, "Opposite exit"),
    (3, "Lateral exit"),
    (4, "Stop / no-valid-exit"),
]


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
        needed = [
            "primaryExitClass",
            "firstDirectionReversalStep",
            "firstDirectionReversalDepthNm",
            "firstDirectionReversalEnergyEv",
            "firstBoundaryStep",
            "firstBoundaryDepthNm",
            "firstBoundaryEnergyEv",
            "firstBoundaryType",
        ]
        ok = all(bool(tree.GetListOfBranches().FindObject(b)) for b in needed)
        if not ok:
            f.Close()
            continue
        meta.GetEntry(0)
        runs.append(
            {
                "path": path,
                "step_nm": float(meta.maxStepNm),
                "e0_eV": float(meta.primaryEnergyMeV) * 1.0e6,
                "events": int(meta.primaryElectrons),
                "model": _meta_string(meta.emModel),
                "thickness_nm": float(meta.sampleThicknessNm),
                "file": f,
                "tree": tree,
            }
        )
    runs.sort(key=lambda r: r["step_nm"])
    return runs


def _selection(metric, cls):
    if metric in ("firstDirectionReversalStep", "firstBoundaryStep"):
        return f"primaryExitClass=={cls} && {metric}>0"
    if metric in ("firstDirectionReversalDepthNm", "firstBoundaryDepthNm"):
        return f"primaryExitClass=={cls} && {metric}>=0"
    if metric in ("firstDirectionReversalEnergyEv", "firstBoundaryEnergyEv"):
        return f"primaryExitClass=={cls} && {metric}>0"
    if metric == "firstBoundaryType":
        return f"primaryExitClass=={cls} && {metric}>0"
    return f"primaryExitClass=={cls}"


def _draw_metric_overlays(runs, metric, title, x_title, bins, x_min, x_max, out_base):
    c = ROOT.TCanvas(
        f"c_{metric}_{ROOT.gRandom.Integer(10_000_000)}",
        f"{title} step overlays",
        1600,
        1100,
    )
    c.Divide(2, 2)
    keep = [c]

    colors = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7]

    for ipad, (cls, cls_label) in enumerate(CLASS_DEFS, start=1):
        c.cd(ipad)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.03)
        ROOT.gPad.SetBottomMargin(0.13)
        ROOT.gPad.SetTopMargin(0.08)
        is_logy = metric != "firstBoundaryType"
        ROOT.gPad.SetLogy(is_logy)

        leg = ROOT.TLegend(0.52, 0.60, 0.94, 0.92)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.032)
        keep.append(leg)

        hists = []
        y_max = 0.0
        for i, run in enumerate(runs):
            h = ROOT.TH1D(
                f"h_{metric}_cls{cls}_step{i}_{ROOT.gRandom.Integer(10_000_000)}",
                "",
                bins,
                x_min,
                x_max,
            )
            n = run["tree"].Draw(f"{metric}>>{h.GetName()}", _selection(metric, cls), "goff")
            h.SetDirectory(0)
            if n > 0:
                y_max = max(y_max, h.GetMaximum())
                color = colors[i % len(colors)]
                h.SetLineColor(color)
                h.SetLineWidth(4)
                leg.AddEntry(h, f"step={run['step_nm']:.3g} nm (N={n})", "l")
                hists.append(h)
            keep.append(h)

        if hists:
            y_min = 1e-1 if is_logy else 0.0
            y_top = max(1.3 * y_max, 1.0)
            if is_logy:
                y_top = max(y_top, 10.0)
            h0 = hists[0]
            h0.SetTitle(f"{cls_label};{x_title};Entries")
            h0.GetYaxis().SetRangeUser(y_min, y_top)
            h0.Draw("HIST")
            for h in hists[1:]:
                h.Draw("HIST SAME")
            leg.Draw()
        else:
            if is_logy:
                ROOT.gPad.SetLogy(False)
            txt = ROOT.TLatex()
            txt.SetNDC(True)
            txt.SetTextFont(42)
            txt.SetTextSize(0.045)
            txt.DrawLatex(0.30, 0.50, "No entries")
            keep.append(txt)

    title_box = ROOT.TPaveText(0.10, 0.955, 0.90, 0.995, "NDC")
    title_box.SetFillStyle(0)
    title_box.SetBorderSize(0)
    title_box.SetTextFont(42)
    title_box.SetTextSize(0.03)
    title_box.AddText(title)
    title_box.Draw()
    keep.append(title_box)

    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".png")


def main():
    parser = argparse.ArgumentParser(description="Overlay transition diagnostics across max-step values.")
    parser.add_argument("--results-dir", required=True, help="Folder with ROOT files.")
    parser.add_argument("--output-dir", required=True, help="Output folder.")
    parser.add_argument("--label", default="diag", help="Output label.")
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

    metric_defs = [
        (
            "firstDirectionReversalStep",
            f"First z-reversal step vs max step (E0={e0:.0f} eV, {model}, MC events={events})",
            "Step index",
            140,
            0.0,
            280.0,
        ),
        (
            "firstDirectionReversalDepthNm",
            f"First z-reversal depth vs max step (E0={e0:.0f} eV, {model}, MC events={events})",
            "Depth from entrance (nm)",
            100,
            0.0,
            max(5.0, thickness),
        ),
        (
            "firstDirectionReversalEnergyEv",
            f"First z-reversal energy vs max step (E0={e0:.0f} eV, {model}, MC events={events})",
            "Kinetic energy (eV)",
            140,
            0.0,
            max(1.0, 1.02 * e0),
        ),
        (
            "firstBoundaryStep",
            f"First outward boundary step vs max step (E0={e0:.0f} eV, {model}, MC events={events})",
            "Step index",
            140,
            0.0,
            280.0,
        ),
        (
            "firstBoundaryDepthNm",
            f"First outward boundary depth vs max step (E0={e0:.0f} eV, {model}, MC events={events})",
            "Depth from entrance (nm)",
            100,
            0.0,
            max(5.0, thickness),
        ),
        (
            "firstBoundaryEnergyEv",
            f"First outward boundary energy vs max step (E0={e0:.0f} eV, {model}, MC events={events})",
            "Kinetic energy (eV)",
            140,
            0.0,
            max(1.0, 1.02 * e0),
        ),
        (
            "firstBoundaryType",
            f"First outward boundary type vs max step (E0={e0:.0f} eV, {model}, MC events={events})",
            "Boundary type code",
            4,
            0.5,
            4.5,
        ),
    ]

    outs = []
    for metric, title, xtitle, bins, xmin, xmax in metric_defs:
        out_base = os.path.join(args.output_dir, f"{metric}_step_overlays_{args.label}")
        _draw_metric_overlays(runs, metric, title, xtitle, bins, xmin, xmax, out_base)
        outs.append(out_base + ".pdf")
        outs.append(out_base + ".png")

    for run in runs:
        run["file"].Close()

    print("Saved:")
    for out in outs:
        print(out)


if __name__ == "__main__":
    main()
