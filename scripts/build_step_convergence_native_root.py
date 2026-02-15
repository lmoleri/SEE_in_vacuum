#!/usr/bin/env python3
import argparse
import csv
import os

import ROOT


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
    return sorted(files)


def _normalize(hist):
    h = hist.Clone(f"{hist.GetName()}_norm")
    h.SetDirectory(0)
    integral = h.Integral()
    if integral > 0.0:
        h.Scale(1.0 / integral)
    return h


def _nonzero_xrange(hist):
    first = hist.FindFirstBinAbove(0.0)
    last = hist.FindLastBinAbove(0.0)
    xaxis = hist.GetXaxis()
    if first <= 0 or last <= 0:
        return xaxis.GetXmin(), xaxis.GetXmax()
    return xaxis.GetBinLowEdge(first), xaxis.GetBinUpEdge(last)


def _color_for_step(step_nm):
    if abs(step_nm - 0.05) < 1e-8:
        return ROOT.kRed + 1
    if abs(step_nm - 0.1) < 1e-8:
        return ROOT.kBlue + 1
    if abs(step_nm - 0.2) < 1e-8:
        return ROOT.kGreen + 2
    return ROOT.kBlack


def _step_tag(step_nm):
    return f"{step_nm:.3g}".replace(".", "p")


def _load_run_hists(results_dir):
    runs = []
    for path in _find_root_files(results_dir):
        f = ROOT.TFile.Open(path)
        if not f or not f.IsOpen():
            continue

        meta = f.Get("RunMeta")
        h_edep = f.Get("EdepPrimary")
        h_exit = f.Get("PrimaryExitEnergyOpposite")
        if not meta or not h_edep or not h_exit:
            f.Close()
            continue

        meta.GetEntry(0)
        step_nm = float(meta.maxStepNm)
        e0_eV = float(meta.primaryEnergyMeV) * 1.0e6
        events = int(meta.primaryElectrons)
        model = _meta_string(meta.emModel)

        h_edep_n = _normalize(h_edep)
        h_exit_n = _normalize(h_exit)
        h_edep_n.SetName(f"EdepPrimary_norm_step{_step_tag(step_nm)}")
        h_exit_n.SetName(f"PrimaryExitEnergyOpposite_norm_step{_step_tag(step_nm)}")

        runs.append(
            {
                "path": path,
                "step_nm": step_nm,
                "e0_eV": e0_eV,
                "events": events,
                "model": model,
                "edep": h_edep_n,
                "exit_opp": h_exit_n,
            }
        )
        f.Close()

    runs.sort(key=lambda r: r["step_nm"])
    return runs


def _draw_overlay(runs, key, canvas_name, title, x_title, out_dir, out_stem):
    c = ROOT.TCanvas(canvas_name, canvas_name, 1200, 850)
    c.SetLeftMargin(0.11)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.11)
    ROOT.gStyle.SetOptStat(0)

    xmin = None
    xmax = None
    ymax = 0.0
    for run in runs:
        h = run[key]
        x0, x1 = _nonzero_xrange(h)
        xmin = x0 if xmin is None else min(xmin, x0)
        xmax = x1 if xmax is None else max(xmax, x1)
        ymax = max(ymax, h.GetMaximum())

    if xmin is None or xmax is None:
        xmin, xmax = 0.0, 1.0
    dx = xmax - xmin
    if dx <= 0.0:
        dx = max(1.0, xmax)
    xmin = max(0.0, xmin - 0.05 * dx)
    xmax = xmax + 0.05 * dx
    if ymax <= 0.0:
        ymax = 1.0

    frame = ROOT.TH1D(f"{canvas_name}_frame", f"{title};{x_title};Arbitrary units (area=1)", 100, xmin, xmax)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.15 * ymax)
    frame.Draw()

    leg = ROOT.TLegend(0.63, 0.62, 0.96, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.032)

    keep = [frame, leg]
    for run in runs:
        h = run[key]
        color = _color_for_step(run["step_nm"])
        h.SetLineColor(color)
        h.SetLineWidth(3)
        h.SetMarkerColor(color)
        h.SetMarkerStyle(0)
        h.Draw("HIST SAME")
        leg.AddEntry(h, f"step={run['step_nm']:.3g} nm", "l")
        keep.append(h)

    leg.Draw()

    e0 = runs[0]["e0_eV"] if runs else 0.0
    model = runs[0]["model"] if runs else ""
    events = runs[0]["events"] if runs else 0
    info = ROOT.TPaveText(0.12, 0.78, 0.47, 0.90, "NDC")
    info.SetFillStyle(0)
    info.SetBorderSize(1)
    info.SetTextFont(42)
    info.SetTextSize(0.03)
    info.AddText(f"E0 = {e0:.0f} eV")
    info.AddText(f"EM model: {model}")
    info.AddText(f"MC events: {events}")
    info.Draw()
    keep.append(info)

    c.SaveAs(os.path.join(out_dir, f"{out_stem}.pdf"))
    c.SaveAs(os.path.join(out_dir, f"{out_stem}.png"))
    c._keep = keep
    return c


def _read_qoi_table(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "step_nm": float(row["step_nm"]),
                    "E0_eV": float(row["E0_eV"]),
                    "f_stop": float(row["f_stop"]),
                    "f_ent": float(row["f_ent"]),
                    "f_opp": float(row["f_opp"]),
                    "mean_edep": float(row["mean_edep"]),
                    "mode_edep": float(row["mode_edep"]),
                    "q50_eexit_opp": float(row["q50_eexit_opp"]),
                    "q90_eexit_opp": float(row["q90_eexit_opp"]),
                    "sey": float(row["sey"]),
                }
            )
    rows.sort(key=lambda r: r["step_nm"])
    return rows


def _draw_qoi_vs_step(qoi_rows, runs, out_dir, out_stem):
    c = ROOT.TCanvas("qoi_vs_step", "qoi_vs_step", 1500, 900)
    c.Divide(2, 1)
    keep = [c]

    steps = [r["step_nm"] for r in qoi_rows]
    x = ROOT.std.vector("double")()
    for s in steps:
        x.push_back(s)

    # Left panel: fractions + SEY
    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.03)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTopMargin(0.10)

    yvals_left = []
    for key in ("f_stop", "f_ent", "f_opp", "sey"):
        yvals_left.extend([r[key] for r in qoi_rows])
    ymax_left = max(yvals_left) if yvals_left else 1.0
    xmin = min(steps) - 0.02
    xmax = max(steps) + 0.02
    frame1 = ROOT.TH1D("frame_qoi_left", "Fractions and SEY vs max step;max step in Al_{2}O_{3} (nm);Value", 100, xmin, xmax)
    frame1.SetMinimum(0.0)
    frame1.SetMaximum(1.05 * ymax_left if ymax_left > 0 else 1.0)
    frame1.Draw()
    keep.append(frame1)

    defs_left = [
        ("f_stop", "f_{stop}", ROOT.kRed + 1, 20),
        ("f_ent", "f_{ent}", ROOT.kOrange + 7, 22),
        ("f_opp", "f_{opp}", ROOT.kBlue + 1, 23),
        ("sey", "SEY", ROOT.kGreen + 2, 33),
    ]
    leg1 = ROOT.TLegend(0.62, 0.18, 0.95, 0.48)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.035)
    keep.append(leg1)

    for key, label, color, marker in defs_left:
        y = ROOT.std.vector("double")()
        for r in qoi_rows:
            y.push_back(r[key])
        g = ROOT.TGraph(len(qoi_rows), x.data(), y.data())
        g.SetName(f"g_qoi_{key}")
        g.SetLineColor(color)
        g.SetMarkerColor(color)
        g.SetMarkerStyle(marker)
        g.SetMarkerSize(1.2)
        g.SetLineWidth(3)
        g.Draw("LP SAME")
        leg1.AddEntry(g, label, "lp")
        keep.append(g)
    leg1.Draw()

    # Right panel: energy metrics
    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.03)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTopMargin(0.10)
    yvals_right = []
    for key in ("mean_edep", "mode_edep", "q50_eexit_opp", "q90_eexit_opp"):
        yvals_right.extend([r[key] for r in qoi_rows])
    ymin_right = min(yvals_right) if yvals_right else 0.0
    ymax_right = max(yvals_right) if yvals_right else 1.0
    yr = ymax_right - ymin_right
    if yr <= 0:
        yr = max(1.0, ymax_right)
    frame2 = ROOT.TH1D("frame_qoi_right", "Energy metrics vs max step;max step in Al_{2}O_{3} (nm);Energy (eV)", 100, xmin, xmax)
    frame2.SetMinimum(max(0.0, ymin_right - 0.08 * yr))
    frame2.SetMaximum(ymax_right + 0.08 * yr)
    frame2.Draw()
    keep.append(frame2)

    defs_right = [
        ("mean_edep", "mean(EdepPrimary)", ROOT.kRed + 1, 20),
        ("mode_edep", "mode(EdepPrimary)", ROOT.kMagenta + 1, 22),
        ("q50_eexit_opp", "q50(E_{exit,opp})", ROOT.kBlue + 1, 23),
        ("q90_eexit_opp", "q90(E_{exit,opp})", ROOT.kGreen + 2, 33),
    ]
    leg2 = ROOT.TLegend(0.50, 0.18, 0.95, 0.50)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.032)
    keep.append(leg2)
    for key, label, color, marker in defs_right:
        y = ROOT.std.vector("double")()
        for r in qoi_rows:
            y.push_back(r[key])
        g = ROOT.TGraph(len(qoi_rows), x.data(), y.data())
        g.SetName(f"g_qoi_{key}")
        g.SetLineColor(color)
        g.SetMarkerColor(color)
        g.SetMarkerStyle(marker)
        g.SetMarkerSize(1.2)
        g.SetLineWidth(3)
        g.Draw("LP SAME")
        leg2.AddEntry(g, label, "lp")
        keep.append(g)
    leg2.Draw()

    # Global title/info
    e0 = qoi_rows[0]["E0_eV"] if qoi_rows else 0.0
    model = runs[0]["model"] if runs else ""
    events = runs[0]["events"] if runs else 0
    title = ROOT.TPaveText(0.08, 0.94, 0.92, 0.995, "NDC")
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    title.SetTextFont(42)
    title.SetTextSize(0.03)
    title.AddText(f"QoI vs max step (E_{{0}}={e0:.0f} eV, {model}, MC events={events})")
    title.Draw()
    keep.append(title)

    c.SaveAs(os.path.join(out_dir, f"{out_stem}.pdf"))
    c.SaveAs(os.path.join(out_dir, f"{out_stem}.png"))
    c._keep = keep
    return c


def _read_pairwise(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "step_coarse": float(row["step_coarse_nm"]),
                    "step_fine": float(row["step_fine_nm"]),
                    "metric": row["metric"],
                    "delta": float(row["delta"]),
                    "delta_type": row["delta_type"],
                    "threshold": float(row["threshold"]),
                }
            )
    return rows


def _draw_pairwise(pair_rows, out_dir, out_stem):
    rel_metrics = ["f_stop", "f_ent", "f_opp", "sey"]
    abs_metrics = ["mean_edep", "mode_edep", "q50_eexit_opp", "q90_eexit_opp"]

    pairs = []
    for r in pair_rows:
        p = (r["step_coarse"], r["step_fine"])
        if p not in pairs:
            pairs.append(p)

    c = ROOT.TCanvas("pairwise_deltas", "pairwise_deltas", 1500, 900)
    c.Divide(2, 1)
    keep = [c]

    colors = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1]
    markers = [20, 22, 23, 33]

    # Left: relative deltas
    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.03)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetTopMargin(0.10)

    x_left = list(range(1, len(rel_metrics) + 1))
    y_rel_max = 1e-6
    for r in pair_rows:
        if r["delta_type"] == "relative" and r["metric"] in rel_metrics:
            y_rel_max = max(y_rel_max, r["delta"])
    frame1 = ROOT.TH1D("frame_pair_rel", "#Delta_{rel};Metric;#Delta_{rel}", len(rel_metrics), 0.5, len(rel_metrics) + 0.5)
    for i, m in enumerate(rel_metrics, start=1):
        frame1.GetXaxis().SetBinLabel(i, m)
    frame1.SetMinimum(0.0)
    frame1.SetMaximum(1.15 * y_rel_max)
    frame1.Draw()
    keep.append(frame1)

    leg1 = ROOT.TLegend(0.58, 0.78, 0.95, 0.93)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.032)
    keep.append(leg1)

    for ip, pair in enumerate(pairs):
        ys = []
        for metric in rel_metrics:
            val = float("nan")
            for r in pair_rows:
                if (
                    r["delta_type"] == "relative"
                    and r["metric"] == metric
                    and abs(r["step_coarse"] - pair[0]) < 1e-12
                    and abs(r["step_fine"] - pair[1]) < 1e-12
                ):
                    val = r["delta"]
                    break
            ys.append(val if val == val else 0.0)
        g = ROOT.TGraph(len(x_left))
        for i, xv in enumerate(x_left):
            g.SetPoint(i, float(xv), float(ys[i]))
        g.SetName(f"g_rel_{ip}")
        g.SetLineColor(colors[ip % len(colors)])
        g.SetMarkerColor(colors[ip % len(colors)])
        g.SetMarkerStyle(markers[ip % len(markers)])
        g.SetMarkerSize(1.2)
        g.SetLineWidth(3)
        g.Draw("LP SAME")
        leg1.AddEntry(g, f"{pair[0]:.3g} -> {pair[1]:.3g} nm", "lp")
        keep.append(g)
    leg1.Draw()

    # Right: absolute deltas
    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.03)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetTopMargin(0.10)
    x_right = list(range(1, len(abs_metrics) + 1))
    y_abs_max = 1e-6
    thr = 10.0
    for r in pair_rows:
        if r["delta_type"] == "absolute" and r["metric"] in abs_metrics:
            y_abs_max = max(y_abs_max, r["delta"])
            thr = r["threshold"]
    frame2 = ROOT.TH1D("frame_pair_abs", "#Delta_{abs};Metric;#Delta_{abs} (eV)", len(abs_metrics), 0.5, len(abs_metrics) + 0.5)
    for i, m in enumerate(abs_metrics, start=1):
        frame2.GetXaxis().SetBinLabel(i, m)
    frame2.SetMinimum(0.0)
    frame2.SetMaximum(1.15 * y_abs_max)
    frame2.Draw()
    keep.append(frame2)

    line = ROOT.TLine(0.5, thr, len(abs_metrics) + 0.5, thr)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.SetLineColor(ROOT.kGray + 2)
    line.Draw("SAME")
    keep.append(line)

    leg2 = ROOT.TLegend(0.55, 0.74, 0.95, 0.93)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.032)
    keep.append(leg2)
    for ip, pair in enumerate(pairs):
        ys = []
        for metric in abs_metrics:
            val = float("nan")
            for r in pair_rows:
                if (
                    r["delta_type"] == "absolute"
                    and r["metric"] == metric
                    and abs(r["step_coarse"] - pair[0]) < 1e-12
                    and abs(r["step_fine"] - pair[1]) < 1e-12
                ):
                    val = r["delta"]
                    break
            ys.append(val if val == val else 0.0)
        g = ROOT.TGraph(len(x_right))
        for i, xv in enumerate(x_right):
            g.SetPoint(i, float(xv), float(ys[i]))
        g.SetName(f"g_abs_{ip}")
        g.SetLineColor(colors[ip % len(colors)])
        g.SetMarkerColor(colors[ip % len(colors)])
        g.SetMarkerStyle(markers[ip % len(markers)])
        g.SetMarkerSize(1.2)
        g.SetLineWidth(3)
        g.Draw("LP SAME")
        leg2.AddEntry(g, f"{pair[0]:.3g} -> {pair[1]:.3g} nm", "lp")
        keep.append(g)
    leg2.AddEntry(line, f"threshold ({thr:.0f} eV)", "l")
    leg2.Draw()

    c.SaveAs(os.path.join(out_dir, f"{out_stem}.pdf"))
    c.SaveAs(os.path.join(out_dir, f"{out_stem}.png"))
    c._keep = keep
    return c


def main():
    parser = argparse.ArgumentParser(
        description="Regenerate step-convergence plots as native ROOT objects and save in one ROOT file."
    )
    parser.add_argument("--results-dir", required=True, help="Directory with step-scan ROOT files.")
    parser.add_argument("--qoi-csv", required=True, help="QoI table CSV from extract_step_convergence_qoi.py.")
    parser.add_argument("--pairwise-csv", required=True, help="Pairwise deltas CSV from extract_step_convergence_qoi.py.")
    parser.add_argument("--output-dir", required=True, help="Output plot directory.")
    parser.add_argument("--label", default="penelope_800eV", help="Label suffix used in plot names.")
    parser.add_argument("--output-root", default="", help="Output ROOT file path. Default: <output-dir>/step_convergence_<label>_plots.root")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    runs = _load_run_hists(args.results_dir)
    if not runs:
        raise SystemExit("No valid ROOT runs found.")

    qoi_rows = _read_qoi_table(args.qoi_csv)
    pair_rows = _read_pairwise(args.pairwise_csv)
    if not qoi_rows:
        raise SystemExit("QoI CSV is empty.")
    if not pair_rows:
        raise SystemExit("Pairwise CSV is empty.")

    os.makedirs(args.output_dir, exist_ok=True)
    output_root = args.output_root
    if not output_root:
        output_root = os.path.join(args.output_dir, f"step_convergence_{args.label}_plots.root")

    c_edep = _draw_overlay(
        runs,
        "edep",
        "edep_primary_overlay",
        "EdepPrimary normalized shape vs step",
        "EdepPrimary (eV)",
        args.output_dir,
        f"edep_primary_overlay_{args.label}",
    )
    c_exit = _draw_overlay(
        runs,
        "exit_opp",
        "exit_energy_opp_overlay",
        "PrimaryExitEnergyOpposite normalized shape vs step",
        "Opposite-side exit energy (eV)",
        args.output_dir,
        f"exit_energy_opp_overlay_{args.label}",
    )
    c_qoi = _draw_qoi_vs_step(qoi_rows, runs, args.output_dir, f"qoi_vs_step_{args.label}")
    c_pair = _draw_pairwise(pair_rows, args.output_dir, f"pairwise_deltas_{args.label}")

    out = ROOT.TFile(output_root, "RECREATE")
    for run in runs:
        out.cd()
        run["edep"].Write()
        run["exit_opp"].Write()
    c_edep.Write("edep_primary_overlay_" + args.label)
    c_exit.Write("exit_energy_opp_overlay_" + args.label)
    c_qoi.Write("qoi_vs_step_" + args.label)
    c_pair.Write("pairwise_deltas_" + args.label)
    out.Close()

    print("Saved native ROOT plot file:")
    print(output_root)


if __name__ == "__main__":
    main()

