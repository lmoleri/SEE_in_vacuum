#!/usr/bin/env python3
import argparse
import csv
import math
import os

import ROOT


TARGET_CLASSES = (2, 4)
CLASS_LABEL = {
    2: "Class 2 (opposite exit)",
    4: "Class 4 (stop / no-valid-exit)",
}


def _decode_cstr(value):
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


def _color_for_step(step_nm):
    if abs(step_nm - 0.1) < 1e-9:
        return ROOT.kBlue + 1
    if abs(step_nm - 0.2) < 1e-9:
        return ROOT.kGreen + 2
    if abs(step_nm - 0.3) < 1e-9:
        return ROOT.kMagenta + 1
    return ROOT.kBlack


def _load_run(path):
    f = ROOT.TFile.Open(path)
    if not f or not f.IsOpen():
        return None
    meta = f.Get("RunMeta")
    tree = f.Get("PrimaryTrajectoryDiagnostics")
    if not meta or not tree:
        f.Close()
        return None

    needed = ["eventId", "primaryExitClass", "deltaThetaDeg", "process"]
    ok = all(bool(tree.GetListOfBranches().FindObject(b)) for b in needed)
    if not ok:
        f.Close()
        return None

    meta.GetEntry(0)
    step_nm = float(meta.maxStepNm)
    e0_eV = float(meta.primaryEnergyMeV) * 1.0e6
    model = _decode_cstr(meta.emModel)
    events = int(meta.primaryElectrons)

    data = {
        cls: {
            "all_theta": [],
            "msc_theta": [],
            "event_max_all": {},
            "event_max_msc": {},
            "event_ids": set(),
        }
        for cls in TARGET_CLASSES
    }

    nrows = tree.GetEntries()
    for i in range(nrows):
        tree.GetEntry(i)
        cls = int(tree.primaryExitClass)
        if cls not in TARGET_CLASSES:
            continue
        theta = float(tree.deltaThetaDeg)
        if not math.isfinite(theta) or theta < 0.0:
            continue
        theta = min(theta, 180.0)
        eid = int(tree.eventId)
        proc = _decode_cstr(tree.process)

        cls_data = data[cls]
        cls_data["event_ids"].add(eid)
        cls_data["all_theta"].append(theta)

        prev = cls_data["event_max_all"].get(eid, -1.0)
        if theta > prev:
            cls_data["event_max_all"][eid] = theta

        if proc == "msc":
            cls_data["msc_theta"].append(theta)
            prev_msc = cls_data["event_max_msc"].get(eid, -1.0)
            if theta > prev_msc:
                cls_data["event_max_msc"][eid] = theta

    f.Close()
    return {
        "path": path,
        "step_nm": step_nm,
        "e0_eV": e0_eV,
        "model": model,
        "events": events,
        "data": data,
    }


def _normalize(hist):
    integral = hist.Integral()
    if integral > 0:
        hist.Scale(1.0 / integral)


def _draw_delta_theta_overlay(runs, output_dir, label):
    c = ROOT.TCanvas("c_delta_theta_overlay", "Delta theta overlays", 1600, 800)
    c.Divide(2, 1)
    keep = [c]

    style_all = ROOT.TLine(0, 0, 1, 0)
    style_all.SetLineColor(ROOT.kBlack)
    style_all.SetLineWidth(3)
    style_msc = ROOT.TLine(0, 0, 1, 0)
    style_msc.SetLineColor(ROOT.kBlack)
    style_msc.SetLineStyle(2)
    style_msc.SetLineWidth(3)
    keep.extend([style_all, style_msc])

    for ipad, cls in enumerate(TARGET_CLASSES, start=1):
        c.cd(ipad)
        ROOT.gPad.SetLeftMargin(0.10)
        ROOT.gPad.SetRightMargin(0.03)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTopMargin(0.10)
        ROOT.gPad.SetLogy(True)

        h_all_list = []
        h_msc_list = []
        y_max = 0.0

        for run in runs:
            step_tag = str(run["step_nm"]).replace(".", "p")
            h_all = ROOT.TH1D(
                f"h_dtheta_all_cls{cls}_step{step_tag}",
                "",
                180,
                0.0,
                180.0,
            )
            h_msc = ROOT.TH1D(
                f"h_dtheta_msc_cls{cls}_step{step_tag}",
                "",
                180,
                0.0,
                180.0,
            )
            h_all.SetDirectory(0)
            h_msc.SetDirectory(0)

            for val in run["data"][cls]["all_theta"]:
                h_all.Fill(val)
            for val in run["data"][cls]["msc_theta"]:
                h_msc.Fill(val)

            _normalize(h_all)
            _normalize(h_msc)
            y_max = max(y_max, h_all.GetMaximum(), h_msc.GetMaximum())
            h_all_list.append((run, h_all))
            h_msc_list.append((run, h_msc))
            keep.extend([h_all, h_msc])

        if y_max <= 0:
            y_max = 1.0

        frame = ROOT.TH2D(
            f"frame_dtheta_cls{cls}",
            "",
            10,
            0.0,
            180.0,
            10,
            1.0e-5,
            y_max * 5.0,
        )
        frame.SetDirectory(0)
        frame.SetTitle(f";#Delta#theta (deg);Arbitrary units (area=1)")
        frame.Draw("AXIS")
        keep.append(frame)

        cls_txt = ROOT.TLatex()
        cls_txt.SetNDC(True)
        cls_txt.SetTextFont(42)
        cls_txt.SetTextSize(0.040)
        cls_txt.DrawLatex(0.12, 0.92, CLASS_LABEL[cls])
        keep.append(cls_txt)

        leg_step = ROOT.TLegend(0.53, 0.67, 0.96, 0.92)
        leg_step.SetBorderSize(0)
        leg_step.SetFillStyle(0)
        leg_step.SetTextFont(42)
        leg_step.SetTextSize(0.030)
        keep.append(leg_step)

        leg_style = ROOT.TLegend(0.53, 0.50, 0.96, 0.64)
        leg_style.SetBorderSize(0)
        leg_style.SetFillStyle(0)
        leg_style.SetTextFont(42)
        leg_style.SetTextSize(0.030)
        keep.append(leg_style)

        for (run, h_all), (_, h_msc) in zip(h_all_list, h_msc_list):
            color = _color_for_step(run["step_nm"])
            h_all.SetLineColor(color)
            h_all.SetLineWidth(3)
            h_all.SetLineStyle(1)
            h_msc.SetLineColor(color)
            h_msc.SetLineWidth(3)
            h_msc.SetLineStyle(2)
            h_all.Draw("HIST SAME")
            h_msc.Draw("HIST SAME")
            n_evt = len(run["data"][cls]["event_ids"])
            leg_step.AddEntry(h_all, f"step={run['step_nm']:.3g} nm (Nsample={n_evt})", "l")

        leg_style.AddEntry(style_all, "all processes", "l")
        leg_style.AddEntry(style_msc, "msc only", "l")
        leg_step.Draw()
        leg_style.Draw()

    title = ROOT.TPaveText(0.10, 0.955, 0.90, 0.995, "NDC")
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    title.SetTextFont(42)
    title.SetTextSize(0.030)
    title.AddText(
        f"#Delta#theta distributions: E0={runs[0]['e0_eV']:.0f} eV, EM={runs[0]['model']}, MC events={runs[0]['events']}"
    )
    title.Draw()
    keep.append(title)

    base = os.path.join(output_dir, f"delta_theta_all_vs_msc_class24_{label}")
    c.SaveAs(base + ".pdf")
    c.SaveAs(base + ".png")
    return c, keep, [base + ".pdf", base + ".png"]


def _draw_threshold_fraction(runs, thresholds, output_dir, label):
    c = ROOT.TCanvas("c_dtheta_threshold_fraction", "Delta theta threshold fractions", 1600, 800)
    c.Divide(2, 1)
    keep = [c]
    csv_rows = []

    for ipad, cls in enumerate(TARGET_CLASSES, start=1):
        c.cd(ipad)
        ROOT.gPad.SetLeftMargin(0.10)
        ROOT.gPad.SetRightMargin(0.03)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTopMargin(0.10)

        frame = ROOT.TH2D(
            f"frame_frac_cls{cls}",
            "",
            10,
            min(thresholds) - 5.0,
            max(thresholds) + 5.0,
            10,
            0.0,
            1.02,
        )
        frame.SetDirectory(0)
        frame.SetTitle(";#Delta#theta threshold (deg);Event fraction with #geq 1 step above threshold")
        frame.Draw("AXIS")
        keep.append(frame)

        cls_txt = ROOT.TLatex()
        cls_txt.SetNDC(True)
        cls_txt.SetTextFont(42)
        cls_txt.SetTextSize(0.040)
        cls_txt.DrawLatex(0.12, 0.92, CLASS_LABEL[cls])
        keep.append(cls_txt)

        leg = ROOT.TLegend(0.53, 0.69, 0.96, 0.92)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.030)
        keep.append(leg)

        for run in runs:
            n_evt = len(run["data"][cls]["event_ids"])
            g = ROOT.TGraph(len(thresholds))
            g.SetName(f"g_frac_cls{cls}_step{str(run['step_nm']).replace('.', 'p')}")
            g.SetLineColor(_color_for_step(run["step_nm"]))
            g.SetMarkerColor(_color_for_step(run["step_nm"]))
            g.SetLineWidth(3)
            g.SetMarkerStyle(20)
            g.SetMarkerSize(1.3)
            keep.append(g)

            for i, thr in enumerate(thresholds):
                if n_evt > 0:
                    n_pass_all = sum(1 for val in run["data"][cls]["event_max_all"].values() if val >= thr)
                    frac_all = float(n_pass_all) / float(n_evt)
                    n_pass_msc = sum(1 for val in run["data"][cls]["event_max_msc"].values() if val >= thr)
                    frac_msc = float(n_pass_msc) / float(n_evt)
                else:
                    frac_all = 0.0
                    frac_msc = 0.0
                g.SetPoint(i, float(thr), frac_all)
                csv_rows.append(
                    {
                        "step_nm": run["step_nm"],
                        "class": cls,
                        "threshold_deg": thr,
                        "fraction_all": frac_all,
                        "fraction_msc_only": frac_msc,
                        "n_sampled_events": n_evt,
                    }
                )

            g.Draw("LP SAME")
            leg.AddEntry(g, f"step={run['step_nm']:.3g} nm (Nsample={n_evt})", "lp")

        leg.Draw()

    title = ROOT.TPaveText(0.10, 0.955, 0.90, 0.995, "NDC")
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    title.SetTextFont(42)
    title.SetTextSize(0.030)
    title.AddText(
        f"Per-event #Delta#theta threshold fractions (sampled trajectories): E0={runs[0]['e0_eV']:.0f} eV, EM={runs[0]['model']}"
    )
    title.Draw()
    keep.append(title)

    base = os.path.join(output_dir, f"delta_theta_threshold_event_fraction_class24_{label}")
    c.SaveAs(base + ".pdf")
    c.SaveAs(base + ".png")

    csv_path = os.path.join(output_dir, f"delta_theta_threshold_event_fraction_class24_{label}.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=[
                "step_nm",
                "class",
                "threshold_deg",
                "fraction_all",
                "fraction_msc_only",
                "n_sampled_events",
            ],
        )
        writer.writeheader()
        for row in csv_rows:
            writer.writerow(row)

    return c, keep, [base + ".pdf", base + ".png", csv_path]


def _write_root(output_dir, label, canvases, objects):
    root_path = os.path.join(output_dir, f"delta_theta_diagnostics_class24_{label}.root")
    out = ROOT.TFile(root_path, "RECREATE")
    for obj in objects:
        if hasattr(obj, "Write"):
            obj.Write()
    for canv in canvases:
        canv.Write()
    out.Close()
    return root_path


def main():
    parser = argparse.ArgumentParser(
        description="Plot delta-theta diagnostics from PrimaryTrajectoryDiagnostics."
    )
    parser.add_argument("--results-dir", required=True, help="Directory with ROOT files.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--label", default="penelope_800eV", help="Output label suffix.")
    parser.add_argument(
        "--thresholds",
        default="10,20,30,45,60,90",
        help="Comma-separated angle thresholds in degrees.",
    )
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.TH1.AddDirectory(False)

    thresholds = [float(tok.strip()) for tok in args.thresholds.split(",") if tok.strip()]
    if not thresholds:
        raise SystemExit("No valid thresholds provided.")
    thresholds = sorted(thresholds)

    runs = []
    for path in _find_root_files(args.results_dir):
        run = _load_run(path)
        if run:
            runs.append(run)
    if not runs:
        raise SystemExit("No usable ROOT files found with PrimaryTrajectoryDiagnostics.")

    runs.sort(key=lambda r: r["step_nm"])

    os.makedirs(args.output_dir, exist_ok=True)

    c1, keep1, out1 = _draw_delta_theta_overlay(runs, args.output_dir, args.label)
    c2, keep2, out2 = _draw_threshold_fraction(runs, thresholds, args.output_dir, args.label)
    root_path = _write_root(args.output_dir, args.label, [c1, c2], keep1 + keep2)

    print("Saved:")
    for p in out1 + out2:
        print(p)
    print(root_path)


if __name__ == "__main__":
    main()
