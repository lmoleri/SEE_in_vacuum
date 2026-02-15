#!/usr/bin/env python3
import argparse
import os
from array import array

import ROOT


CLASS_DEFS = [
    (2, "Opposite exit"),
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


def _format_step_tag(step_nm):
    text = f"{step_nm:.3g}".replace(".", "p")
    return text


def _collect_runs(results_dir):
    runs = []
    for path in _find_root_files(results_dir):
        f = ROOT.TFile.Open(path)
        if not f or not f.IsOpen():
            continue
        meta = f.Get("RunMeta")
        tree = f.Get("PrimaryTrajectoryDiagnostics")
        if not meta or not tree:
            f.Close()
            continue
        meta.GetEntry(0)
        runs.append(
            {
                "file_path": path,
                "file": f,
                "tree": tree,
                "step_nm": float(meta.maxStepNm),
                "e0_eV": float(meta.primaryEnergyMeV) * 1.0e6,
                "thickness_nm": float(meta.sampleThicknessNm),
                "events": int(meta.primaryElectrons),
                "model": _meta_string(meta.emModel),
            }
        )
    runs.sort(key=lambda r: r["step_nm"])
    return runs


def _collect_class_traces(tree, exit_class, max_events):
    selected = {}
    order = []
    n_entries = tree.GetEntries()
    for i in range(n_entries):
        tree.GetEntry(i)
        if int(tree.primaryExitClass) != exit_class:
            continue
        event_id = int(tree.eventId)
        if event_id not in selected:
            if len(order) >= max_events:
                continue
            selected[event_id] = []
            order.append(event_id)
        selected[event_id].append(
            {
                "step": int(tree.stepNumber),
                "pre_depth": float(tree.preDepthNm),
                "post_depth": float(tree.postDepthNm),
                "pre_e": float(tree.preEnergyEv),
                "post_e": float(tree.postEnergyEv),
                "process": _meta_string(tree.process),
                "reversal": int(tree.reversalOnStep),
                "first_reversal": int(tree.isFirstReversalStep),
            }
        )

    traces = []
    for event_id in order:
        steps = sorted(selected[event_id], key=lambda s: s["step"])
        if not steps:
            continue
        x_vals = [steps[0]["pre_depth"]]
        y_vals = [steps[0]["pre_e"]]
        for s in steps:
            x_vals.append(s["post_depth"])
            y_vals.append(s["post_e"])

        first_msc = None
        first_rev = None
        for s in steps:
            if first_rev is None and s["first_reversal"] > 0:
                first_rev = s
            if first_msc is None and s["process"] == "msc" and s["reversal"] > 0:
                first_msc = s
            if first_msc is not None and first_rev is not None:
                break

        traces.append(
            {
                "event_id": event_id,
                "x": x_vals,
                "y": y_vals,
                "first_msc": first_msc,
                "first_rev": first_rev,
            }
        )
    return traces


def _draw_class_panel(pad, class_label, traces, thickness_nm, e0_eV):
    pad.cd()
    pad.SetLeftMargin(0.12)
    pad.SetRightMargin(0.03)
    pad.SetBottomMargin(0.12)
    pad.SetTopMargin(0.10)
    pad.SetGrid(0, 0)

    frame = ROOT.TH2D(
        f"frame_{ROOT.gRandom.Integer(10000000)}",
        f";Depth from entrance (nm);Kinetic energy (eV)",
        100,
        0.0,
        max(5.0, thickness_nm),
        100,
        0.0,
        max(1.0, 1.02 * e0_eV),
    )
    frame.SetStats(0)
    frame.Draw("AXIS")

    keep = [frame]
    class_title = ROOT.TLatex()
    class_title.SetNDC(True)
    class_title.SetTextFont(42)
    class_title.SetTextSize(0.05)
    class_title.DrawLatex(0.13, 0.94, class_label)
    keep.append(class_title)

    colors = [
        ROOT.kBlue + 1,
        ROOT.kGreen + 2,
        ROOT.kAzure + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
        ROOT.kCyan + 2,
        ROOT.kViolet + 5,
        ROOT.kTeal + 3,
    ]

    n_with_msc = 0
    n_with_rev = 0
    first_msc_depth_sum = 0.0
    first_msc_energy_sum = 0.0
    first_rev_depth_sum = 0.0
    first_rev_energy_sum = 0.0

    for i, trace in enumerate(traces):
        g = ROOT.TGraph(len(trace["x"]), array("d", trace["x"]), array("d", trace["y"]))
        g.SetLineColor(colors[i % len(colors)])
        g.SetLineWidth(2)
        g.Draw("L SAME")
        keep.append(g)

        first_rev = trace["first_rev"]
        if first_rev is not None:
            n_with_rev += 1
            first_rev_depth_sum += first_rev["post_depth"]
            first_rev_energy_sum += first_rev["post_e"]
            m_rev = ROOT.TMarker(first_rev["post_depth"], first_rev["post_e"], 29)
            m_rev.SetMarkerColor(ROOT.kMagenta + 2)
            m_rev.SetMarkerSize(0.7)
            m_rev.Draw("SAME")
            keep.append(m_rev)

        first_msc = trace["first_msc"]
        if first_msc is not None:
            n_with_msc += 1
            first_msc_depth_sum += first_msc["post_depth"]
            first_msc_energy_sum += first_msc["post_e"]
            m_msc = ROOT.TMarker(first_msc["post_depth"], first_msc["post_e"], 20)
            m_msc.SetMarkerColor(ROOT.kRed + 1)
            m_msc.SetMarkerSize(0.55)
            m_msc.Draw("SAME")
            keep.append(m_msc)

    stats = {
        "n_samples": len(traces),
        "n_with_rev": n_with_rev,
        "n_with_msc": n_with_msc,
        "mean_msc_depth": (first_msc_depth_sum / n_with_msc) if n_with_msc > 0 else None,
        "mean_msc_energy": (first_msc_energy_sum / n_with_msc) if n_with_msc > 0 else None,
        "mean_rev_depth": (first_rev_depth_sum / n_with_rev) if n_with_rev > 0 else None,
        "mean_rev_energy": (first_rev_energy_sum / n_with_rev) if n_with_rev > 0 else None,
    }
    return keep, stats


def _draw_run_canvas(run, output_dir, label, max_events_per_class):
    traces_by_class = {}
    for cls, _ in CLASS_DEFS:
        traces_by_class[cls] = _collect_class_traces(run["tree"], cls, max_events_per_class)

    step_tag = _format_step_tag(run["step_nm"])
    out_base = os.path.join(
        output_dir,
        f"sampled_primary_trajectories_depth_energy_{label}_step{step_tag}nm",
    )

    c = ROOT.TCanvas(
        f"c_traj_{step_tag}_{ROOT.gRandom.Integer(10000000)}",
        "Sampled primary trajectories",
        1900,
        950,
    )

    top = ROOT.TPad(f"top_{step_tag}", "", 0.00, 0.93, 1.00, 1.00)
    left = ROOT.TPad(f"left_{step_tag}", "", 0.00, 0.22, 0.50, 0.93)
    right = ROOT.TPad(f"right_{step_tag}", "", 0.50, 0.22, 1.00, 0.93)
    bottom = ROOT.TPad(f"bottom_{step_tag}", "", 0.00, 0.00, 1.00, 0.22)
    top.Draw()
    left.Draw()
    right.Draw()
    bottom.Draw()

    keep = [c]
    keep.extend([top, left, right, bottom])

    class_stats = {}
    for ipad, (cls, cls_label) in enumerate(CLASS_DEFS, start=1):
        panel = left if ipad == 1 else right
        panel_keep, stats = _draw_class_panel(
            panel,
            f"{cls_label} (class {cls})",
            traces_by_class[cls],
            run["thickness_nm"],
            run["e0_eV"],
        )
        class_stats[cls] = stats
        keep.extend(panel_keep)

    top.cd()
    top.SetLeftMargin(0.03)
    top.SetRightMargin(0.03)
    top.SetBottomMargin(0.02)
    top.SetTopMargin(0.02)
    title = ROOT.TPaveText(0.01, 0.05, 0.99, 0.95, "NDC")
    title.SetFillStyle(0)
    title.SetBorderSize(1)
    title.SetTextFont(42)
    title.SetTextSize(0.26)
    title.AddText(
        f"Sampled full trajectories: E0={run['e0_eV']:.0f} eV, step={run['step_nm']:.3g} nm, "
        f"EM={run['model']}, MC events={run['events']}"
    )
    title.Draw()
    keep.append(title)

    bottom.cd()
    bottom.SetLeftMargin(0.02)
    bottom.SetRightMargin(0.02)
    bottom.SetTopMargin(0.06)
    bottom.SetBottomMargin(0.10)

    info2 = ROOT.TPaveText(0.01, 0.08, 0.37, 0.95, "NDC")
    info2.SetFillStyle(0)
    info2.SetBorderSize(1)
    info2.SetTextFont(42)
    info2.SetTextAlign(12)
    info2.SetTextSize(0.19)
    s2 = class_stats[2]
    info2.AddText("Class 2 (Opposite exit)")
    info2.AddText(f"Sampled events: {s2['n_samples']}")
    info2.AddText(f"With first z-reversal: {s2['n_with_rev']}")
    info2.AddText(f"With first msc reversal: {s2['n_with_msc']}")
    if s2["mean_msc_depth"] is not None:
        info2.AddText(f"<depth(first msc rev)> = {s2['mean_msc_depth']:.2f} nm")
    if s2["mean_msc_energy"] is not None:
        info2.AddText(f"<energy(first msc rev)> = {s2['mean_msc_energy']:.0f} eV")
    info2.Draw()
    keep.append(info2)

    info4 = ROOT.TPaveText(0.38, 0.08, 0.74, 0.95, "NDC")
    info4.SetFillStyle(0)
    info4.SetBorderSize(1)
    info4.SetTextFont(42)
    info4.SetTextAlign(12)
    info4.SetTextSize(0.19)
    s4 = class_stats[4]
    info4.AddText("Class 4 (Stop / no-valid-exit)")
    info4.AddText(f"Sampled events: {s4['n_samples']}")
    info4.AddText(f"With first z-reversal: {s4['n_with_rev']}")
    info4.AddText(f"With first msc reversal: {s4['n_with_msc']}")
    if s4["mean_msc_depth"] is not None:
        info4.AddText(f"<depth(first msc rev)> = {s4['mean_msc_depth']:.2f} nm")
    if s4["mean_msc_energy"] is not None:
        info4.AddText(f"<energy(first msc rev)> = {s4['mean_msc_energy']:.0f} eV")
    info4.Draw()
    keep.append(info4)

    legend = ROOT.TLegend(0.76, 0.10, 0.99, 0.95)
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.17)
    line_demo = ROOT.TLine(0, 0, 1, 1)
    line_demo.SetLineColor(ROOT.kBlue + 1)
    line_demo.SetLineWidth(2)
    legend.AddEntry(line_demo, "Sampled trajectory", "l")
    marker_msc = ROOT.TMarker(0, 0, 20)
    marker_msc.SetMarkerColor(ROOT.kRed + 1)
    legend.AddEntry(marker_msc, "First msc reversal step", "p")
    marker_rev = ROOT.TMarker(0, 0, 29)
    marker_rev.SetMarkerColor(ROOT.kMagenta + 2)
    legend.AddEntry(marker_rev, "First z-reversal step", "p")
    legend.Draw()
    keep.extend([legend, line_demo, marker_msc, marker_rev])

    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".png")
    return c, keep, out_base


def main():
    parser = argparse.ArgumentParser(
        description="Plot sampled full primary-electron trajectories (class 2 vs 4)."
    )
    parser.add_argument("--results-dir", required=True, help="Folder with ROOT files.")
    parser.add_argument("--output-dir", required=True, help="Output folder.")
    parser.add_argument("--label", default="trajdiag", help="Output label.")
    parser.add_argument(
        "--max-events-per-class",
        type=int,
        default=30,
        help="Maximum number of sampled events drawn per class (default: 30).",
    )
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    runs = _collect_runs(args.results_dir)
    if not runs:
        raise SystemExit(
            "No ROOT files with PrimaryTrajectoryDiagnostics found. "
            "Enable trajectory diagnostics in config and rerun."
        )

    os.makedirs(args.output_dir, exist_ok=True)
    out_root = os.path.join(
        args.output_dir, f"sampled_primary_trajectories_depth_energy_{args.label}.root"
    )
    fout = ROOT.TFile(out_root, "RECREATE")

    outputs = []
    all_keep = []
    for run in runs:
        canvas, keep, out_base = _draw_run_canvas(
            run,
            args.output_dir,
            args.label,
            args.max_events_per_class,
        )
        all_keep.extend(keep)
        fout.cd()
        canvas.Write(canvas.GetName())
        outputs.append(out_base + ".pdf")
        outputs.append(out_base + ".png")

    fout.Close()
    for run in runs:
        run["file"].Close()

    print("Saved:")
    for path in outputs:
        print(path)
    print(out_root)


if __name__ == "__main__":
    main()
