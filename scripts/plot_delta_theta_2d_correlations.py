#!/usr/bin/env python3
import argparse
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


def _load_run(path):
    f = ROOT.TFile.Open(path)
    if not f or not f.IsOpen():
        return None

    meta = f.Get("RunMeta")
    tree = f.Get("PrimaryTrajectoryDiagnostics")
    if not meta or not tree:
        f.Close()
        return None

    needed = ["primaryExitClass", "stepNumber", "preEnergyEv", "deltaThetaDeg", "process"]
    ok = all(bool(tree.GetListOfBranches().FindObject(b)) for b in needed)
    if not ok:
        f.Close()
        return None

    meta.GetEntry(0)
    step_nm = float(meta.maxStepNm)
    e0_eV = float(meta.primaryEnergyMeV) * 1.0e6
    model = _decode_cstr(meta.emModel)
    events = int(meta.primaryElectrons)

    data = {cls: [] for cls in TARGET_CLASSES}
    max_step_idx = 1

    nrows = tree.GetEntries()
    for i in range(nrows):
        tree.GetEntry(i)
        cls = int(tree.primaryExitClass)
        if cls not in TARGET_CLASSES:
            continue
        step_idx = int(tree.stepNumber)
        e_cur = float(tree.preEnergyEv)
        theta = float(tree.deltaThetaDeg)
        if step_idx <= 0:
            continue
        if not math.isfinite(e_cur) or e_cur < 0.0:
            continue
        if not math.isfinite(theta) or theta < 0.0:
            continue
        e_cur = min(e_cur, e0_eV)
        theta = min(theta, 180.0)
        is_msc = 1 if _decode_cstr(tree.process) == "msc" else 0
        data[cls].append((step_idx, e_cur, theta, is_msc))
        max_step_idx = max(max_step_idx, step_idx)

    f.Close()
    return {
        "path": path,
        "step_nm": step_nm,
        "e0_eV": e0_eV,
        "events": events,
        "model": model,
        "max_step_idx": max_step_idx,
        "data": data,
    }


def _draw_canvas(runs, x_mode, proc_mode, output_dir, label):
    ncols = len(runs)
    nrows = len(TARGET_CLASSES)
    c = ROOT.TCanvas(
        f"c_2d_{x_mode}_{proc_mode}",
        f"2D delta theta maps ({x_mode}, {proc_mode})",
        max(1200, 560 * ncols),
        980,
    )
    c.Divide(ncols, nrows)
    keep = [c]

    x_title = "Step index" if x_mode == "step" else "Current primary energy (eV)"
    x_desc = "step index" if x_mode == "step" else "current primary energy"
    proc_desc = "all processes" if proc_mode == "all" else "msc only"

    e0 = runs[0]["e0_eV"]
    global_max_step = max(r["max_step_idx"] for r in runs)
    global_max_step = max(global_max_step, 10)

    for row, cls in enumerate(TARGET_CLASSES):
        for col, run in enumerate(runs):
            ipad = row * ncols + col + 1
            c.cd(ipad)
            ROOT.gPad.SetLeftMargin(0.12)
            ROOT.gPad.SetRightMargin(0.16)
            ROOT.gPad.SetBottomMargin(0.13)
            ROOT.gPad.SetTopMargin(0.10)
            ROOT.gPad.SetLogz(True)

            if x_mode == "step":
                x_min = 0.5
                x_max = float(global_max_step) + 0.5
                n_bins_x = int(global_max_step)
            else:
                x_min = 0.0
                x_max = max(1.0, 1.001 * e0)
                n_bins_x = 120

            h = ROOT.TH2D(
                f"h2_{x_mode}_{proc_mode}_cls{cls}_step{str(run['step_nm']).replace('.', 'p')}",
                "",
                n_bins_x,
                x_min,
                x_max,
                180,
                0.0,
                180.0,
            )
            h.SetDirectory(0)
            keep.append(h)

            n_fill = 0
            for step_idx, e_cur, theta, is_msc in run["data"][cls]:
                if proc_mode == "msc" and not is_msc:
                    continue
                x_val = float(step_idx) if x_mode == "step" else e_cur
                h.Fill(x_val, theta)
                n_fill += 1

            h.SetTitle(f";{x_title};#Delta#theta (deg)")
            if n_fill > 0:
                h.SetMinimum(1.0)
                h.Draw("COLZ")
            else:
                ROOT.gPad.SetLogz(False)
                h.Draw("COLZ")
                txt_empty = ROOT.TLatex()
                txt_empty.SetNDC(True)
                txt_empty.SetTextFont(42)
                txt_empty.SetTextSize(0.05)
                txt_empty.DrawLatex(0.35, 0.50, "No entries")
                keep.append(txt_empty)

            txt = ROOT.TLatex()
            txt.SetNDC(True)
            txt.SetTextFont(42)
            txt.SetTextSize(0.035)
            txt.DrawLatex(0.14, 0.93, CLASS_LABEL[cls])
            txt.DrawLatex(0.14, 0.87, f"step = {run['step_nm']:.3g} nm")
            txt.DrawLatex(0.14, 0.81, f"N steps = {n_fill}")
            keep.append(txt)

    top = ROOT.TPaveText(0.06, 0.955, 0.94, 0.995, "NDC")
    top.SetFillStyle(0)
    top.SetBorderSize(0)
    top.SetTextFont(42)
    top.SetTextSize(0.028)
    top.AddText(
        f"2D #Delta#theta map vs {x_desc} ({proc_desc}): E0={e0:.0f} eV, EM={runs[0]['model']}, MC events={runs[0]['events']}"
    )
    top.Draw()
    keep.append(top)

    base = os.path.join(output_dir, f"delta_theta_2d_vs_{x_mode}_{proc_mode}_class24_{label}")
    c.SaveAs(base + ".pdf")
    c.SaveAs(base + ".png")
    return c, keep, [base + ".pdf", base + ".png"]


def main():
    parser = argparse.ArgumentParser(
        description="2D correlation plots: delta theta vs step index or current energy."
    )
    parser.add_argument("--results-dir", required=True, help="Directory with ROOT files.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--label", default="penelope_800eV", help="Output label suffix.")
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

    canvases = []
    all_keep = []
    outputs = []
    for x_mode in ("step", "energy"):
        for proc_mode in ("all", "msc"):
            c, keep, outs = _draw_canvas(runs, x_mode, proc_mode, args.output_dir, args.label)
            canvases.append(c)
            all_keep.extend(keep)
            outputs.extend(outs)

    root_path = os.path.join(args.output_dir, f"delta_theta_2d_correlations_class24_{args.label}.root")
    fout = ROOT.TFile(root_path, "RECREATE")
    for obj in all_keep:
        if hasattr(obj, "Write"):
            obj.Write()
    for c in canvases:
        c.Write()
    fout.Close()

    print("Saved:")
    for p in outputs:
        print(p)
    print(root_path)


if __name__ == "__main__":
    main()
