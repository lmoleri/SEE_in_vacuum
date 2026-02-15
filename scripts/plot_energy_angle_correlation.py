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


def _pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    syy = sum(y * y for y in ys)
    sxy = sum(x * y for x, y in zip(xs, ys))
    num = n * sxy - sx * sy
    denx = n * sxx - sx * sx
    deny = n * syy - sy * sy
    den = denx * deny
    if den <= 0.0:
        return float("nan")
    return num / math.sqrt(den)


def _load_run(path):
    f = ROOT.TFile.Open(path)
    if not f or not f.IsOpen():
        return None

    meta = f.Get("RunMeta")
    tree = f.Get("PrimaryTrajectoryDiagnostics")
    if not meta or not tree:
        f.Close()
        return None

    needed = ["primaryExitClass", "preEnergyEv", "deltaThetaDeg", "process"]
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
        cls: {"all_x": [], "all_y": [], "msc_x": [], "msc_y": []}
        for cls in TARGET_CLASSES
    }

    nrows = tree.GetEntries()
    for i in range(nrows):
        tree.GetEntry(i)
        cls = int(tree.primaryExitClass)
        if cls not in TARGET_CLASSES:
            continue
        e_cur = float(tree.preEnergyEv)
        theta = float(tree.deltaThetaDeg)
        if not math.isfinite(e_cur) or e_cur < 0.0:
            continue
        if not math.isfinite(theta) or theta < 0.0:
            continue
        e_cur = min(e_cur, e0_eV)
        theta = min(theta, 180.0)
        proc = _decode_cstr(tree.process)

        data[cls]["all_x"].append(e_cur)
        data[cls]["all_y"].append(theta)
        if proc == "msc":
            data[cls]["msc_x"].append(e_cur)
            data[cls]["msc_y"].append(theta)

    f.Close()
    return {
        "path": path,
        "step_nm": step_nm,
        "e0_eV": e0_eV,
        "events": events,
        "model": model,
        "data": data,
    }


def _build_profile(xs, ys, n_bins, e_max, name):
    p = ROOT.TProfile(name, "", n_bins, 0.0, e_max)
    p.SetDirectory(0)
    for x, y in zip(xs, ys):
        p.Fill(x, y)
    return p


def main():
    parser = argparse.ArgumentParser(
        description="Current-energy vs scattering-angle diagnostic plots from PrimaryTrajectoryDiagnostics."
    )
    parser.add_argument("--results-dir", required=True, help="Directory with ROOT files.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--label", default="penelope_800eV", help="Output label suffix.")
    parser.add_argument("--n-bins", type=int, default=80, help="Number of energy bins for profiles.")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.TH1.AddDirectory(False)

    runs = []
    for path in _find_root_files(args.results_dir):
        run = _load_run(path)
        if run:
            runs.append(run)
    if not runs:
        raise SystemExit("No usable ROOT files with PrimaryTrajectoryDiagnostics found.")
    runs.sort(key=lambda r: r["step_nm"])

    os.makedirs(args.output_dir, exist_ok=True)
    e_max = max(r["e0_eV"] for r in runs)

    c = ROOT.TCanvas("c_energy_angle_corr", "Current energy vs delta theta", 1800, 900)
    c.Divide(2, 1)
    keep = [c]
    csv_rows = []

    style_all = ROOT.TLine(0, 0, 1, 0)
    style_all.SetLineWidth(3)
    style_all.SetLineStyle(1)
    style_msc = ROOT.TLine(0, 0, 1, 0)
    style_msc.SetLineWidth(3)
    style_msc.SetLineStyle(2)
    keep.extend([style_all, style_msc])

    for ipad, cls in enumerate(TARGET_CLASSES, start=1):
        c.cd(ipad)
        ROOT.gPad.SetLeftMargin(0.10)
        ROOT.gPad.SetRightMargin(0.03)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTopMargin(0.10)

        frame = ROOT.TH2D(
            f"frame_energy_angle_cls{cls}",
            "",
            10,
            0.0,
            e_max,
            10,
            0.0,
            180.0,
        )
        frame.SetDirectory(0)
        frame.SetTitle(";Current primary energy (eV);Mean #Delta#theta (deg)")
        frame.Draw("AXIS")
        keep.append(frame)

        cls_txt = ROOT.TLatex()
        cls_txt.SetNDC(True)
        cls_txt.SetTextFont(42)
        cls_txt.SetTextSize(0.040)
        cls_txt.DrawLatex(0.12, 0.92, CLASS_LABEL[cls])
        keep.append(cls_txt)

        leg_step = ROOT.TLegend(0.52, 0.66, 0.97, 0.92)
        leg_step.SetBorderSize(0)
        leg_step.SetFillStyle(0)
        leg_step.SetTextFont(42)
        leg_step.SetTextSize(0.030)
        keep.append(leg_step)

        leg_style = ROOT.TLegend(0.52, 0.52, 0.97, 0.64)
        leg_style.SetBorderSize(0)
        leg_style.SetFillStyle(0)
        leg_style.SetTextFont(42)
        leg_style.SetTextSize(0.030)
        leg_style.AddEntry(style_all, "all processes", "l")
        leg_style.AddEntry(style_msc, "msc only", "l")
        keep.append(leg_style)

        stats_box = ROOT.TPaveText(0.12, 0.64, 0.48, 0.90, "NDC")
        stats_box.SetFillColor(0)
        stats_box.SetBorderSize(1)
        stats_box.SetTextFont(42)
        stats_box.SetTextAlign(12)
        stats_box.SetTextSize(0.026)
        keep.append(stats_box)

        for run in runs:
            step_tag = str(run["step_nm"]).replace(".", "p")
            all_x = run["data"][cls]["all_x"]
            all_y = run["data"][cls]["all_y"]
            msc_x = run["data"][cls]["msc_x"]
            msc_y = run["data"][cls]["msc_y"]

            p_all = _build_profile(all_x, all_y, args.n_bins, e_max, f"p_all_cls{cls}_step{step_tag}")
            p_msc = _build_profile(msc_x, msc_y, args.n_bins, e_max, f"p_msc_cls{cls}_step{step_tag}")
            keep.extend([p_all, p_msc])

            color = _color_for_step(run["step_nm"])
            p_all.SetLineColor(color)
            p_all.SetLineWidth(3)
            p_all.SetLineStyle(1)
            p_all.SetMarkerSize(0.0)
            p_all.Draw("HIST SAME")

            p_msc.SetLineColor(color)
            p_msc.SetLineWidth(3)
            p_msc.SetLineStyle(2)
            p_msc.SetMarkerSize(0.0)
            p_msc.Draw("HIST SAME")

            leg_step.AddEntry(p_all, f"step={run['step_nm']:.3g} nm", "l")

            rho_all = _pearson(all_x, all_y)
            rho_msc = _pearson(msc_x, msc_y)
            stats_box.AddText(
                f"step={run['step_nm']:.3g}: #rho(all)={rho_all:.3f}, #rho(msc)={rho_msc:.3f}"
            )
            csv_rows.append(
                {
                    "step_nm": run["step_nm"],
                    "class": cls,
                    "pearson_all": rho_all,
                    "pearson_msc_only": rho_msc,
                    "n_steps_all": len(all_x),
                    "n_steps_msc_only": len(msc_x),
                }
            )

        leg_step.Draw()
        leg_style.Draw()
        stats_box.Draw()

    title = ROOT.TPaveText(0.10, 0.955, 0.90, 0.995, "NDC")
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    title.SetTextFont(42)
    title.SetTextSize(0.030)
    title.AddText(
        f"Current-energy vs #Delta#theta correlation: E0={runs[0]['e0_eV']:.0f} eV, EM={runs[0]['model']}, MC events={runs[0]['events']}"
    )
    title.Draw()
    keep.append(title)

    base = os.path.join(args.output_dir, f"current_energy_vs_delta_theta_correlation_class24_{args.label}")
    c.SaveAs(base + ".pdf")
    c.SaveAs(base + ".png")

    csv_path = base + ".csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=[
                "step_nm",
                "class",
                "pearson_all",
                "pearson_msc_only",
                "n_steps_all",
                "n_steps_msc_only",
            ],
        )
        writer.writeheader()
        for row in csv_rows:
            writer.writerow(row)

    root_path = base + ".root"
    fout = ROOT.TFile(root_path, "RECREATE")
    for obj in keep:
        if hasattr(obj, "Write"):
            obj.Write()
    fout.Close()

    print("Saved:")
    print(base + ".pdf")
    print(base + ".png")
    print(csv_path)
    print(root_path)


if __name__ == "__main__":
    main()
