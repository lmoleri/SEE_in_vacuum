#!/usr/bin/env python3
import argparse
import csv
import os
from array import array
from pathlib import Path

import ROOT


CLASS_LABELS = {
    1: "Entrance exit",
    2: "Opposite exit",
    3: "Lateral exit",
    4: "Stop / no-valid-exit",
}

CLASS_COLORS = {
    1: ROOT.kRed + 1,
    2: ROOT.kBlue + 1,
    3: ROOT.kGreen + 2,
    4: ROOT.kMagenta + 1,
}

CLASS_MARKERS = {
    1: 20,
    2: 22,
    3: 33,
    4: 21,
}


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


def _safe_div(a, b):
    if b == 0:
        return float("nan")
    return float(a) / float(b)


def _has_branch(tree, branch):
    return bool(tree.GetListOfBranches().FindObject(branch))


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


def _load_rows(results_dir):
    rows = []
    for path in _find_root_files(results_dir):
        f = ROOT.TFile.Open(path)
        if not f or not f.IsOpen():
            continue

        meta = f.Get("RunMeta")
        evt = f.Get("EventDiagnostics")
        if not meta or not evt:
            f.Close()
            continue

        meta.GetEntry(0)
        step_nm = float(meta.maxStepNm)
        energy_eV = float(meta.primaryEnergyMeV) * 1.0e6
        n_events = int(meta.primaryElectrons)
        em_model = _meta_string(meta.emModel)

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
        if not all(_has_branch(evt, b) for b in needed):
            f.Close()
            continue

        stats = {}
        for cls in (1, 2, 3, 4):
            stats[cls] = {
                "n": 0,
                "rev_n": 0,
                "rev_step_sum": 0.0,
                "rev_depth_sum": 0.0,
                "rev_energy_sum": 0.0,
                "bnd_n": 0,
                "bnd_step_sum": 0.0,
                "bnd_depth_sum": 0.0,
                "bnd_energy_sum": 0.0,
                "bnd_type": {1: 0, 2: 0, 3: 0, 4: 0},
            }

        n_total = int(evt.GetEntries())
        for i in range(n_total):
            evt.GetEntry(i)
            cls = int(evt.primaryExitClass)
            if cls not in stats:
                continue
            s = stats[cls]
            s["n"] += 1

            if int(evt.firstDirectionReversalStep) > 0:
                s["rev_n"] += 1
                s["rev_step_sum"] += float(evt.firstDirectionReversalStep)
                s["rev_depth_sum"] += float(evt.firstDirectionReversalDepthNm)
                s["rev_energy_sum"] += float(evt.firstDirectionReversalEnergyEv)

            if int(evt.firstBoundaryStep) > 0:
                s["bnd_n"] += 1
                s["bnd_step_sum"] += float(evt.firstBoundaryStep)
                s["bnd_depth_sum"] += float(evt.firstBoundaryDepthNm)
                s["bnd_energy_sum"] += float(evt.firstBoundaryEnergyEv)
                btype = int(evt.firstBoundaryType)
                if btype in s["bnd_type"]:
                    s["bnd_type"][btype] += 1

        for cls in (1, 2, 3, 4):
            s = stats[cls]
            n = s["n"]
            rev_n = s["rev_n"]
            bnd_n = s["bnd_n"]
            rows.append(
                {
                    "path": path,
                    "step_nm": step_nm,
                    "energy_eV": energy_eV,
                    "em_model": em_model,
                    "n_events": n_events,
                    "exit_class": cls,
                    "class_label": CLASS_LABELS[cls],
                    "fraction": _safe_div(n, n_total),
                    "reversal_fraction": _safe_div(rev_n, n),
                    "first_reversal_step_mean": _safe_div(s["rev_step_sum"], rev_n),
                    "first_reversal_depth_nm_mean": _safe_div(s["rev_depth_sum"], rev_n),
                    "first_reversal_energy_eV_mean": _safe_div(s["rev_energy_sum"], rev_n),
                    "first_boundary_fraction": _safe_div(bnd_n, n),
                    "first_boundary_step_mean": _safe_div(s["bnd_step_sum"], bnd_n),
                    "first_boundary_depth_nm_mean": _safe_div(s["bnd_depth_sum"], bnd_n),
                    "first_boundary_energy_eV_mean": _safe_div(s["bnd_energy_sum"], bnd_n),
                    "first_boundary_type1_frac": _safe_div(s["bnd_type"][1], bnd_n),
                    "first_boundary_type2_frac": _safe_div(s["bnd_type"][2], bnd_n),
                    "first_boundary_type3_frac": _safe_div(s["bnd_type"][3], bnd_n),
                    "first_boundary_type4_frac": _safe_div(s["bnd_type"][4], bnd_n),
                }
            )

        f.Close()
    return rows


def _write_csv(rows, csv_path):
    fields = [
        "path",
        "step_nm",
        "energy_eV",
        "em_model",
        "n_events",
        "exit_class",
        "class_label",
        "fraction",
        "reversal_fraction",
        "first_reversal_step_mean",
        "first_reversal_depth_nm_mean",
        "first_reversal_energy_eV_mean",
        "first_boundary_fraction",
        "first_boundary_step_mean",
        "first_boundary_depth_nm_mean",
        "first_boundary_energy_eV_mean",
        "first_boundary_type1_frac",
        "first_boundary_type2_frac",
        "first_boundary_type3_frac",
        "first_boundary_type4_frac",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["step_nm"], x["exit_class"])):
            w.writerow(r)


def _finite(v):
    try:
        return (v == v) and (v != float("inf")) and (v != float("-inf"))
    except Exception:
        return False


def _plot_metric(pad, rows, metric_key, ytitle, ymin=None, ymax=None):
    pad.cd()
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.03)
    pad.SetBottomMargin(0.13)
    pad.SetTopMargin(0.08)

    steps = sorted({float(r["step_nm"]) for r in rows})
    xmin = min(steps) * 0.9
    xmax = max(steps) * 1.1
    vals = [float(r[metric_key]) for r in rows if _finite(float(r[metric_key]))]
    if ymin is None:
        ymin = 0.0 if not vals else min(0.0, min(vals))
    if ymax is None:
        ymax = 1.0 if not vals else max(vals) * 1.15
    if ymax <= ymin:
        ymax = ymin + 1.0

    frame = ROOT.TH2D(f"frame_{metric_key}_{id(pad)}", "", 10, xmin, xmax, 10, ymin, ymax)
    frame.SetTitle(f";max step in Al_{{2}}O_{{3}} (nm);{ytitle}")
    frame.Draw("AXIS")

    leg = ROOT.TLegend(0.52, 0.62, 0.94, 0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)

    keep = [frame, leg]
    for cls in (1, 2, 3, 4):
        data = [r for r in rows if int(r["exit_class"]) == cls]
        data.sort(key=lambda x: float(x["step_nm"]))
        x = array("d")
        y = array("d")
        for d in data:
            val = float(d[metric_key])
            if _finite(val):
                x.append(float(d["step_nm"]))
                y.append(val)
        if len(x) == 0:
            continue
        g = ROOT.TGraph(len(x), x, y)
        g.SetLineColor(CLASS_COLORS[cls])
        g.SetMarkerColor(CLASS_COLORS[cls])
        g.SetMarkerStyle(CLASS_MARKERS[cls])
        g.SetMarkerSize(1.1)
        g.SetLineWidth(2)
        g.Draw("LP SAME")
        leg.AddEntry(g, CLASS_LABELS[cls], "lp")
        keep.append(g)
    leg.Draw()
    return keep


def main():
    parser = argparse.ArgumentParser(
        description="Plot first-branch diagnostics (reversal/boundary) from EventDiagnostics."
    )
    parser.add_argument("--results-dir", required=True, help="Directory with ROOT outputs.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--label", default="transition_diag", help="Output label suffix.")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    rows = _load_rows(args.results_dir)
    if not rows:
        raise SystemExit("No ROOT files with EventDiagnostics + transition branches found.")

    energies = sorted({float(r["energy_eV"]) for r in rows})
    if len(energies) != 1:
        raise SystemExit("Mixed energies found. Use one energy campaign per call.")

    e0 = energies[0]
    model = rows[0]["em_model"]
    n_events = int(rows[0]["n_events"])

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"transition_metrics_{args.label}.csv"
    _write_csv(rows, csv_path)

    c = ROOT.TCanvas("c_transition", "Transition diagnostics", 1500, 1050)
    c.Divide(2, 2)

    keep = [c]
    keep.extend(_plot_metric(c.cd(1), rows, "fraction", "Class fraction", 0.0, 1.05))
    keep.extend(_plot_metric(c.cd(2), rows, "reversal_fraction", "P(first z-reversal | class)", 0.0, 1.05))
    keep.extend(_plot_metric(c.cd(3), rows, "first_reversal_step_mean", "Mean first-reversal step index"))
    keep.extend(_plot_metric(c.cd(4), rows, "first_boundary_step_mean", "Mean first-boundary step index"))

    title = ROOT.TPaveText(0.12, 0.955, 0.88, 0.995, "NDC")
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    title.SetTextFont(42)
    title.SetTextSize(0.03)
    title.AddText(
        f"Transition diagnostics vs step (E0={e0:.0f} eV, EM model: {model}, MC events: {n_events})"
    )
    title.Draw()
    keep.append(title)

    out_base = out_dir / f"transition_diagnostics_vs_step_{args.label}"
    c.SaveAs(str(out_base) + ".pdf")
    c.SaveAs(str(out_base) + ".png")

    print(f"Saved CSV: {csv_path}")
    print(f"Saved plots: {out_base}.pdf/.png")


if __name__ == "__main__":
    main()
