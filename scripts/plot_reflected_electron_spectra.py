#!/usr/bin/env python3
import argparse
import math
import os

import ROOT


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
    if abs(step_nm - 0.4) < 1e-9:
        return ROOT.kOrange + 7
    if abs(step_nm - 0.5) < 1e-9:
        return ROOT.kRed + 1
    if abs(step_nm - 0.05) < 1e-9:
        return ROOT.kRed + 1
    return ROOT.kBlack


def _normalize(hist):
    integral = hist.Integral()
    if integral > 0:
        hist.Scale(1.0 / integral)


def _get_attr(obj, name, default):
    try:
        return getattr(obj, name)
    except Exception:
        return default


def _build_eloss_hist(h_exit, e0_eV, name):
    n_bins = h_exit.GetNbinsX()
    h_loss = ROOT.TH1D(name, "", n_bins, 0.0, e0_eV)
    h_loss.SetDirectory(0)
    for i in range(1, n_bins + 1):
        count = h_exit.GetBinContent(i)
        if count <= 0:
            continue
        e_exit = h_exit.GetBinCenter(i)
        e_loss = e0_eV - e_exit
        if e_loss < 0.0:
            e_loss = 0.0
        if e_loss > e0_eV:
            e_loss = e0_eV
        h_loss.Fill(e_loss, count)
    return h_loss


def _load_runs(results_dir):
    runs = []
    for path in _find_root_files(results_dir):
        f = ROOT.TFile.Open(path)
        if not f or not f.IsOpen():
            continue
        meta = f.Get("RunMeta")
        h_exit = f.Get("PrimaryExitEnergyEntrance")
        h_exit_spec = f.Get("PrimaryExitEnergyEntranceSpecular")
        h_class = f.Get("PrimaryExitClass")
        if not meta or not h_exit or not h_class:
            f.Close()
            continue
        meta.GetEntry(0)
        step_nm = float(meta.maxStepNm)
        e0_eV = float(meta.primaryEnergyMeV) * 1.0e6
        n_primary = int(meta.primaryElectrons)
        model = _decode_cstr(meta.emModel)
        incidence_surface_deg = float(
            _get_attr(
                meta,
                "incidenceAngleSurfaceDeg",
                math.degrees(math.asin(min(1.0, abs(float(_get_attr(meta, "primaryDirZ", 1.0)))))),
            )
        )
        spec_enabled = bool(int(_get_attr(meta, "specularAcceptanceEnabled", 0)))
        spec_half_angle_deg = float(_get_attr(meta, "specularAcceptanceHalfAngleDeg", 0.0))

        use_specular = bool(h_exit_spec and h_exit_spec.Integral() > 0.0)
        h_reflected_source = h_exit_spec if use_specular else h_exit

        h_exit_c = h_reflected_source.Clone(f"h_exit_step{str(step_nm).replace('.', 'p')}")
        h_exit_c.SetDirectory(0)
        h_loss_c = _build_eloss_hist(h_exit_c, e0_eV, f"h_loss_step{str(step_nm).replace('.', 'p')}")

        if use_specular:
            n_ref = float(h_reflected_source.Integral())
        else:
            # class=1 is ROOT bin 2 for PrimaryExitClass histogram with [0,4], 4 bins
            n_ref = float(h_class.GetBinContent(2))
        f_ref = n_ref / float(n_primary) if n_primary > 0 else 0.0

        runs.append(
            {
                "path": path,
                "step_nm": step_nm,
                "e0_eV": e0_eV,
                "n_primary": n_primary,
                "model": model,
                "n_ref": n_ref,
                "f_ref": f_ref,
                "h_exit": h_exit_c,
                "h_loss": h_loss_c,
                "uses_specular": use_specular,
                "incidence_surface_deg": incidence_surface_deg,
                "spec_enabled": spec_enabled,
                "spec_half_angle_deg": spec_half_angle_deg,
            }
        )
        f.Close()

    runs.sort(key=lambda r: r["step_nm"])
    return runs


def _draw_overlay(runs, key, x_title, title, out_base):
    c = ROOT.TCanvas(f"c_{key}", title, 1000, 760)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.08)

    keep = [c]
    hists = []
    y_max = 0.0
    x_max = max(r["e0_eV"] for r in runs) * 1.001
    for run in runs:
        h = run[key].Clone(f"{run[key].GetName()}_norm")
        h.SetDirectory(0)
        _normalize(h)
        y_max = max(y_max, h.GetMaximum())
        hists.append((run, h))
        keep.append(h)

    if y_max <= 0:
        y_max = 1.0

    frame = ROOT.TH2D(
        f"frame_{key}",
        "",
        10,
        0.0,
        x_max,
        10,
        0.0,
        1.2 * y_max,
    )
    frame.SetDirectory(0)
    frame.SetTitle(f"{title};{x_title};Arbitrary units (area=1)")
    frame.Draw("AXIS")
    keep.append(frame)

    leg = ROOT.TLegend(0.56, 0.64, 0.94, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.032)
    keep.append(leg)

    for run, h in hists:
        color = _color_for_step(run["step_nm"])
        h.SetLineColor(color)
        h.SetLineWidth(3)
        h.Draw("HIST SAME")
        leg.AddEntry(h, f"step={run['step_nm']:.3g} nm  (f_ref={run['f_ref']:.3f})", "l")
    leg.Draw()

    info = ROOT.TPaveText(0.14, 0.78, 0.46, 0.90, "NDC")
    info.SetFillColor(0)
    info.SetBorderSize(1)
    info.SetTextFont(42)
    info.SetTextSize(0.028)
    info.AddText(f"E0 = {runs[0]['e0_eV']:.0f} eV")
    info.AddText(f"EM model: {runs[0]['model']}")
    info.AddText(f"MC events: {runs[0]['n_primary']}")
    info.AddText(f"Incidence: {runs[0]['incidence_surface_deg']:.1f}#circ to surface")
    if any(r.get("uses_specular", False) for r in runs):
        info.AddText(
            f"Selection: entrance-side, spec cone #pm{runs[0]['spec_half_angle_deg']:.1f}#circ"
        )
    info.Draw()
    keep.append(info)

    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".png")
    return c, keep


def _draw_reflected_fraction(runs, out_base):
    c = ROOT.TCanvas("c_reflected_fraction", "Reflected fraction vs step", 1000, 760)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.08)
    keep = [c]

    x_min = min(r["step_nm"] for r in runs) * 0.9
    x_max = max(r["step_nm"] for r in runs) * 1.1
    y_max = max(r["f_ref"] for r in runs) * 1.25
    if y_max <= 0:
        y_max = 1.0

    frame = ROOT.TH2D("frame_reflected_fraction", "", 10, x_min, x_max, 10, 0.0, y_max)
    frame.SetDirectory(0)
    frame.SetTitle("Reflected-electron fraction vs max step;max step in Al_{2}O_{3} (nm);Reflected fraction (class 1)")
    frame.Draw("AXIS")
    keep.append(frame)

    g = ROOT.TGraph(len(runs))
    g.SetLineColor(ROOT.kBlue + 1)
    g.SetMarkerColor(ROOT.kBlue + 1)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.4)
    g.SetLineWidth(3)
    for i, run in enumerate(runs):
        g.SetPoint(i, run["step_nm"], run["f_ref"])
    g.Draw("LP SAME")
    keep.append(g)

    info = ROOT.TPaveText(0.14, 0.78, 0.46, 0.90, "NDC")
    info.SetFillColor(0)
    info.SetBorderSize(1)
    info.SetTextFont(42)
    info.SetTextSize(0.028)
    info.AddText(f"E0 = {runs[0]['e0_eV']:.0f} eV")
    info.AddText(f"EM model: {runs[0]['model']}")
    info.AddText(f"MC events: {runs[0]['n_primary']}")
    info.AddText(f"Incidence: {runs[0]['incidence_surface_deg']:.1f}#circ to surface")
    if any(r.get("uses_specular", False) for r in runs):
        info.AddText(f"Selection: spec cone #pm{runs[0]['spec_half_angle_deg']:.1f}#circ")
    info.Draw()
    keep.append(info)

    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".png")
    return c, keep


def main():
    parser = argparse.ArgumentParser(
        description="Produce reflected-electron-only spectra (class 1 / entrance-exit) for REELS comparison."
    )
    parser.add_argument("--results-dir", required=True, help="Directory containing ROOT outputs.")
    parser.add_argument("--output-dir", required=True, help="Output directory for plots.")
    parser.add_argument("--label", default="reflected", help="Label suffix for output filenames.")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.AddDirectory(False)

    runs = _load_runs(args.results_dir)
    if not runs:
        raise SystemExit(
            "No usable ROOT files found with PrimaryExitEnergyEntrance and PrimaryExitClass."
        )

    os.makedirs(args.output_dir, exist_ok=True)
    outputs = []
    all_keep = []
    canvases = []

    out_base_exit = os.path.join(args.output_dir, f"reflected_exit_energy_shape_{args.label}")
    c_exit, keep_exit = _draw_overlay(
        runs,
        "h_exit",
        "Reflected primary energy at exit (eV)",
        "Reflected-electron exit-energy spectra (class 1 only)",
        out_base_exit,
    )
    outputs.extend([out_base_exit + ".pdf", out_base_exit + ".png"])
    canvases.append(c_exit)
    all_keep.extend(keep_exit)

    out_base_loss = os.path.join(args.output_dir, f"reels_reflected_energy_loss_shape_{args.label}")
    c_loss, keep_loss = _draw_overlay(
        runs,
        "h_loss",
        "Energy loss of reflected primary (eV)",
        "REELS-like reflected energy-loss spectra (class 1 only)",
        out_base_loss,
    )
    outputs.extend([out_base_loss + ".pdf", out_base_loss + ".png"])
    canvases.append(c_loss)
    all_keep.extend(keep_loss)

    out_base_frac = os.path.join(args.output_dir, f"reflected_fraction_vs_step_{args.label}")
    c_frac, keep_frac = _draw_reflected_fraction(runs, out_base_frac)
    outputs.extend([out_base_frac + ".pdf", out_base_frac + ".png"])
    canvases.append(c_frac)
    all_keep.extend(keep_frac)

    root_path = os.path.join(args.output_dir, f"reflected_only_plots_{args.label}.root")
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
