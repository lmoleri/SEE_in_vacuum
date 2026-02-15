#!/usr/bin/env python3
import argparse
import os

import ROOT


CLASS_DEFS = [
    (1, "Entrance exit", ROOT.kRed + 1),
    (2, "Opposite exit", ROOT.kBlue + 1),
    (3, "Lateral exit", ROOT.kGreen + 2),
    (4, "Stop / no-valid-exit", ROOT.kMagenta + 1),
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


def _normalize(h):
    integ = h.Integral()
    if integ > 0:
        h.Scale(1.0 / integ)


def _draw_metric(tree, metric, title, x_title, bins, x_min, x_max, normalize=False):
    is_logy = bool(ROOT.gPad.GetLogy())
    y0 = 1e-4 if is_logy else 0.0
    frame = ROOT.TH2D(
        f"frame_{metric}_{ROOT.gRandom.Integer(10_000_000)}",
        "",
        10,
        x_min,
        x_max,
        10,
        y0,
        1.0,
    )
    frame.SetTitle(f"{title};{x_title};{'Arbitrary units (area=1)' if normalize else 'Entries'}")
    frame.Draw("AXIS")

    keep = [frame]
    legend = ROOT.TLegend(0.50, 0.62, 0.94, 0.92)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)
    keep.append(legend)

    y_max = 0.0
    drawn = False
    for cls, label, color in CLASS_DEFS:
        h = ROOT.TH1D(
            f"h_{metric}_c{cls}_{ROOT.gRandom.Integer(10_000_000)}",
            "",
            bins,
            x_min,
            x_max,
        )
        sel = f"primaryExitClass=={cls} && {metric}>-0.5"
        if metric in ("firstDirectionReversalStep", "firstBoundaryStep"):
            sel = f"primaryExitClass=={cls} && {metric}>0"
        elif metric in ("firstDirectionReversalDepthNm", "firstBoundaryDepthNm"):
            sel = f"primaryExitClass=={cls} && {metric}>=0"
        elif metric in ("firstDirectionReversalEnergyEv", "firstBoundaryEnergyEv"):
            sel = f"primaryExitClass=={cls} && {metric}>0"
        elif metric == "firstBoundaryType":
            sel = f"primaryExitClass=={cls} && {metric}>0"

        draw_expr = f"{metric}>>{h.GetName()}"
        n = tree.Draw(draw_expr, sel, "goff")
        h.SetDirectory(0)
        if n <= 0:
            keep.append(h)
            continue
        if normalize:
            _normalize(h)
        y_max = max(y_max, h.GetMaximum())
        h.SetLineColor(color)
        h.SetMarkerColor(color)
        h.SetLineWidth(3)
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.8)
        h.Draw("HIST SAME")
        legend.AddEntry(h, f"{label} (N={n})", "l")
        keep.append(h)
        drawn = True

    if drawn:
        y_low = 0.0
        if is_logy:
            y_low = max(1e-4, y_max * 1e-4)
        frame.GetYaxis().SetRangeUser(y_low, max(1.25 * y_max, y_low * 10.0))
        frame.Draw("AXIS SAME")
        for obj in keep:
            if isinstance(obj, ROOT.TH1):
                if obj.GetEntries() > 0:
                    obj.Draw("HIST SAME")
        legend.Draw()
    else:
        if is_logy:
            ROOT.gPad.SetLogy(False)
        msg = ROOT.TLatex()
        msg.SetNDC(True)
        msg.SetTextFont(42)
        msg.SetTextSize(0.04)
        msg.DrawLatex(0.22, 0.50, "No entries")
        keep.append(msg)

    return keep


def _plot_file(path, output_dir, label, normalize):
    f = ROOT.TFile.Open(path)
    if not f or not f.IsOpen():
        return []
    meta = f.Get("RunMeta")
    tree = f.Get("EventDiagnostics")
    if not meta or not tree:
        f.Close()
        return []

    meta.GetEntry(0)
    step_nm = float(meta.maxStepNm)
    e0_eV = float(meta.primaryEnergyMeV) * 1.0e6
    thickness_nm = float(meta.sampleThicknessNm)
    model = _meta_string(meta.emModel)
    n_events = int(meta.primaryElectrons)

    c = ROOT.TCanvas(
        f"c_diag_dist_{ROOT.gRandom.Integer(10_000_000)}",
        "Transition diagnostic distributions",
        1900,
        1300,
    )
    c.Divide(3, 3)
    keep = [c]

    metrics = [
        ("firstDirectionReversalStep", "First z-reversal step", "Step index", 120, 0.0, 240.0),
        (
            "firstDirectionReversalDepthNm",
            "First z-reversal depth",
            "Depth from entrance (nm)",
            80,
            0.0,
            max(5.0, thickness_nm),
        ),
        (
            "firstDirectionReversalEnergyEv",
            "First z-reversal energy",
            "Kinetic energy (eV)",
            120,
            0.0,
            max(1.0, 1.02 * e0_eV),
        ),
        ("firstBoundaryStep", "First outward boundary step", "Step index", 120, 0.0, 240.0),
        (
            "firstBoundaryDepthNm",
            "First outward boundary depth",
            "Depth from entrance (nm)",
            80,
            0.0,
            max(5.0, thickness_nm),
        ),
        (
            "firstBoundaryEnergyEv",
            "First outward boundary energy",
            "Kinetic energy (eV)",
            120,
            0.0,
            max(1.0, 1.02 * e0_eV),
        ),
        ("firstBoundaryType", "First outward boundary type", "Type code", 4, 0.5, 4.5),
    ]

    for i, (metric, title, xtitle, bins, xmin, xmax) in enumerate(metrics, start=1):
        c.cd(i)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.03)
        ROOT.gPad.SetBottomMargin(0.14)
        ROOT.gPad.SetTopMargin(0.08)
        # Sparse class-conditioned distributions are easier to inspect in log scale.
        if metric != "firstBoundaryType":
            ROOT.gPad.SetLogy(True)
        else:
            ROOT.gPad.SetLogy(False)
        keep.extend(_draw_metric(tree, metric, title, xtitle, bins, xmin, xmax, normalize=normalize))

    c.cd(8)
    ROOT.gPad.SetLeftMargin(0.10)
    ROOT.gPad.SetRightMargin(0.02)
    ROOT.gPad.SetBottomMargin(0.10)
    ROOT.gPad.SetTopMargin(0.10)
    info = ROOT.TPaveText(0.05, 0.05, 0.95, 0.95, "NDC")
    info.SetFillColor(0)
    info.SetBorderSize(1)
    info.SetTextFont(42)
    info.SetTextSize(0.05)
    info.AddText(f"E0: {e0_eV:.0f} eV")
    info.AddText(f"max step: {step_nm:.3g} nm")
    info.AddText(f"thickness: {thickness_nm:.1f} nm")
    info.AddText(f"EM model: {model}")
    info.AddText(f"MC events: {n_events}")
    info.AddText("Boundary types: 1=Al2O3->World, 3=Al2O3->other")
    info.Draw()
    keep.append(info)

    c.cd(9)
    ROOT.gPad.SetLeftMargin(0.10)
    ROOT.gPad.SetRightMargin(0.02)
    ROOT.gPad.SetBottomMargin(0.10)
    ROOT.gPad.SetTopMargin(0.10)
    note = ROOT.TPaveText(0.05, 0.05, 0.95, 0.95, "NDC")
    note.SetFillColor(0)
    note.SetBorderSize(1)
    note.SetTextFont(42)
    note.SetTextSize(0.05)
    note.AddText("Reversal vars use first z-sign flip")
    note.AddText("Boundary vars use first outward crossing")
    note.AddText("Entries with no such transition are excluded")
    note.Draw()
    keep.append(note)

    os.makedirs(output_dir, exist_ok=True)
    out_base = os.path.join(output_dir, f"transition_diag_distributions_step{str(step_nm).replace('.', 'p')}nm_{label}")
    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".png")

    f.Close()
    return [out_base + ".pdf", out_base + ".png"]


def main():
    parser = argparse.ArgumentParser(description="Plot distributions of EventDiagnostics transition variables.")
    parser.add_argument("--results-dir", required=True, help="Directory with ROOT files.")
    parser.add_argument("--output-dir", required=True, help="Output directory for plots.")
    parser.add_argument("--label", default="diag", help="Suffix label for output file names.")
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize class distributions to unit area (default is raw entries).",
    )
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    files = _find_root_files(args.results_dir)
    if not files:
        raise SystemExit("No ROOT files found.")

    outputs = []
    for path in files:
        outputs.extend(_plot_file(path, args.output_dir, args.label, normalize=args.normalize))

    print("Saved:")
    for out in outputs:
        print(out)


if __name__ == "__main__":
    main()
