#!/usr/bin/env python3
"""
Toy model: events with many crossings of the shell.

Each event has N crossings (configurable). For each crossing, energy deposition
ΔE is sampled from EdepPrimary (Geant4 histogram with good statistics), then
the number of secondary electrons is drawn from Poisson(μ) with μ = (ΔE/ε)*P_esc.
Total SE per event = sum over crossings.

Configuration is read from a JSON file (all parameters that were CLI in the
single-histogram toy MC, plus n_events, crossings_per_event, primary input).
"""

import json
import os
import sys
from math import exp

import numpy as np

# Optional: use shared Poisson sampler from calculate_muon_sey
try:
    from calculate_muon_sey import sample_poisson
except ImportError:
    # Standalone fallback: use ROOT's Poisson
    def sample_poisson(mu, random_gen=None):
        import ROOT
        if random_gen is None:
            random_gen = ROOT.gRandom
        if mu <= 0:
            return 0
        return int(random_gen.Poisson(mu))


def load_config(config_path):
    """Load and validate toy model config (JSON)."""
    with open(config_path, "r") as f:
        cfg = json.load(f)
    # Required
    for key in ("edep_root_file", "histogram", "n_events", "crossings_per_event", "seed",
                "epsilon", "B", "alpha", "depth", "output_dir"):
        if key not in cfg:
            raise ValueError("Config missing key: {}".format(key))
    # crossings_per_event can be int or list
    cpe = cfg["crossings_per_event"]
    if isinstance(cpe, (list, tuple)):
        cfg["_crossings_scan"] = list(cpe)
    else:
        cfg["_crossings_scan"] = [int(cpe)]
    return cfg


def run_toy_events(config_path, cli_overrides=None):
    """
    Run toy events: for each crossings_per_event value, generate n_events events,
    each with that many crossings; sample ΔE from EdepPrimary, SEE from Poisson(μ(ΔE)).
    """
    cfg = load_config(config_path)
    if cli_overrides:
        for k, v in cli_overrides.items():
            if k in cfg:
                cfg[k] = v

    import ROOT
    from ROOT import TFile, TH1D, TCanvas, TLatex, TPaveText, gStyle, TGraph, TGraphErrors

    ROOT.gRandom.SetSeed(cfg["seed"])
    epsilon = float(cfg["epsilon"])
    B = float(cfg["B"])
    alpha = float(cfg["alpha"])
    depth = float(cfg["depth"])
    P_esc = B * exp(-alpha * depth)

    edep_path = cfg["edep_root_file"]
    if not os.path.isabs(edep_path):
        base = os.path.dirname(os.path.abspath(config_path))
        # Resolve relative paths from project root (config may be in config/toy_model/ or config/geant4/)
        while base and os.path.basename(base) in ("toy_model", "geant4", "config"):
            base = os.path.dirname(base)
        edep_path = os.path.normpath(os.path.join(base, edep_path))
    if not os.path.isfile(edep_path):
        raise FileNotFoundError("Edep ROOT file not found: {}".format(edep_path))

    rf = TFile.Open(edep_path, "READ")
    if not rf or rf.IsZombie():
        raise RuntimeError("Cannot open ROOT file: {}".format(edep_path))
    edep_hist = rf.Get(cfg["histogram"])
    if not edep_hist:
        rf.Close()
        raise RuntimeError("Histogram '{}' not found in {}".format(cfg["histogram"], edep_path))
    # Build CDF from histogram for sampling (avoid GetRandom() which can segfault with PyROOT)
    n_bins_edep = edep_hist.GetNbinsX()
    bin_edges_lo = []
    bin_edges_hi = []
    bin_contents = []
    total_content = 0.0
    for i in range(1, n_bins_edep + 1):
        c = edep_hist.GetBinContent(i)
        if c > 0:
            bin_edges_lo.append(edep_hist.GetXaxis().GetBinLowEdge(i))
            bin_edges_hi.append(edep_hist.GetXaxis().GetBinUpEdge(i))
            bin_contents.append(c)
            total_content += c
    rf.Close()
    if total_content <= 0:
        raise RuntimeError("EdepPrimary histogram has no positive bins")
    bin_probs = [x / total_content for x in bin_contents]

    n_events = int(cfg["n_events"])
    crossings_scan = cfg["_crossings_scan"]
    out_dir = cfg["output_dir"]
    os.makedirs(out_dir, exist_ok=True)

    def sample_delta_e():
        u = ROOT.gRandom.Uniform(0.0, 1.0)
        cum = 0.0
        for j, p in enumerate(bin_probs):
            cum += p
            if u <= cum:
                return ROOT.gRandom.Uniform(bin_edges_lo[j], bin_edges_hi[j])
        return bin_edges_hi[-1] if bin_edges_hi else 0.0

    results = {}
    for n_crossings in crossings_scan:
        total_se_per_event = []
        for _ in range(n_events):
            event_se = 0
            for _ in range(n_crossings):
                delta_e = sample_delta_e()
                if delta_e <= 0:
                    continue
                mu = (delta_e / epsilon) * P_esc
                n_se = sample_poisson(mu, ROOT.gRandom)
                event_se += n_se
            total_se_per_event.append(event_se)

        arr = np.array(total_se_per_event)
        n_with_se = int(np.sum(arr >= 1))
        efficiency = n_with_se / n_events if n_events > 0 else 0.0
        max_se = int(np.max(arr))
        n_bins = max_se + 1
        h = TH1D("TotalSE_per_event_n{}".format(n_crossings),
                 "Total SE per event ({} crossings)".format(n_crossings),
                 n_bins, -0.5, max_se + 0.5)
        h.SetXTitle("Total secondary electrons per event")
        h.SetYTitle("Number of events")
        for v in total_se_per_event:
            h.Fill(v)
        results[n_crossings] = {
            "histogram": h,
            "mean": float(np.mean(arr)),
            "std": float(np.std(arr)),
            "min": int(np.min(arr)),
            "max": int(np.max(arr)),
            "n_with_se": n_with_se,
            "efficiency": efficiency,
        }

    summary = "n_events={} crossings_scan={} epsilon={} B={} alpha={} depth={} P_esc={}\n".format(
        n_events, crossings_scan, epsilon, B, alpha, depth, P_esc)
    for n_crossings, data in results.items():
        summary += "crossings={} mean_SE={:.4f} std_SE={:.4f} min={} max={} efficiency={:.4f} (n_with_SE={})\n".format(
            n_crossings, data["mean"], data["std"], data["min"], data["max"],
            data["efficiency"], data["n_with_se"])

    # Save ROOT file with all histograms
    out_root = os.path.join(out_dir, "toy_events_SE_per_event.root")
    of = TFile.Open(out_root, "RECREATE")
    for n_crossings, data in results.items():
        data["histogram"].Write()
    of.Close()

    # Save config copy and text summary
    out_cfg = os.path.join(out_dir, "toy_model_config_used.json")
    with open(out_cfg, "w") as f:
        json.dump({k: v for k, v in cfg.items() if not k.startswith("_")}, f, indent=2)
    out_txt = os.path.join(out_dir, "toy_events_summary.txt")
    with open(out_txt, "w") as f:
        f.write(summary)

    # Simple plot: one canvas per crossings value (before histograms go out of scope)
    if out_dir.startswith("results/"):
        plots_dir = os.path.join("plots", os.path.basename(out_dir))
    else:
        plots_dir = os.path.join("plots", os.path.basename(out_dir.rstrip(os.sep)))
    os.makedirs(plots_dir, exist_ok=True)
    gStyle.SetOptStat(110)
    colors = [1, 2, 4, 6, 8, 9]  # kBlack, kRed, kBlue, kMagenta, kGreen+2, kBlue+1
    for idx, (n_crossings, data) in enumerate(results.items()):
        c = TCanvas("c_{}".format(n_crossings), "Total SE per event ({} crossings)".format(n_crossings), 800, 600)
        c.SetGrid()
        c.SetLogy(1)
        h = data["histogram"]
        h.SetMinimum(0.1)
        h.Draw("HIST")
        lat = TLatex()
        lat.SetNDC(True)
        lat.SetTextSize(0.03)
        lat.DrawLatex(0.15, 0.85, "Toy events: {} crossings per event".format(n_crossings))
        lat.DrawLatex(0.15, 0.80, "Mean total SE = {:.2f}".format(data["mean"]))
        lat.DrawLatex(0.15, 0.75, "Std = {:.2f}".format(data["std"]))
        lat.DrawLatex(0.15, 0.70, "Efficiency (#geq 1 SE) = {:.2f}% ({}/{})".format(
            data["efficiency"] * 100.0, data["n_with_se"], n_events))
        c.Update()
        pdf_path = os.path.join(plots_dir, "TotalSE_per_event_{}_crossings.pdf".format(n_crossings))
        root_path = os.path.join(plots_dir, "TotalSE_per_event_{}_crossings.root".format(n_crossings))
        c.SaveAs(pdf_path)
        c.SaveAs(root_path)
        print("Plot saved: {} and {}".format(pdf_path, root_path))

    # Summary plot: overlay all crossing counts on one canvas
    if len(results) > 1:
        # Use outline (no fill) for overlay so all visible. Fix frame first so ROOT does not replace it when drawing the second histogram.
        order = sorted(results.keys(), key=lambda n: results[n]["max"])
        all_histos = [results[n]["histogram"] for n in order]
        ymax_global = max(hi.GetMaximum() for hi in all_histos)
        xmin_global = min(hi.GetXaxis().GetBinLowEdge(1) for hi in all_histos)
        xmax_global = max(hi.GetXaxis().GetBinUpEdge(hi.GetNbinsX()) for hi in all_histos)
        ymin_frame = 0.1
        ymax_frame = ymax_global * 1.15
        c_summary = TCanvas("c_summary", "Total SE per event (all crossings)", 800, 600)
        c_summary.SetGrid()
        c_summary.SetLogy(1)
        # Set frame first so both histograms are drawn into the same fixed range (avoids only second histogram visible)
        frame_hist = c_summary.DrawFrame(xmin_global, ymin_frame, xmax_global, ymax_frame)
        frame_hist.SetTitle("Total SE per event (%d events)" % n_events)
        frame_hist.SetXTitle("Total secondary electrons per event")
        frame_hist.SetYTitle("Number of events")
        # Legend box in upper-right where histograms tail off (empty space)
        # Ensure box and text stay within plot area (right edge of plot ~0.90 NDC)
        leg_box = TPaveText(0.55, 0.68, 0.88, 0.82, "NDC")
        leg_box.SetBorderSize(1)
        leg_box.SetFillStyle(1001)  # solid fill so legend is readable over grid
        leg_box.SetFillColor(0)    # white background
        leg_box.SetTextAlign(12)   # left-aligned
        leg_box.SetTextFont(42)    # Helvetica, precision 2
        leg_box.SetTextSize(0.028)
        # Keep references to clones so PyROOT does not garbage-collect them before SaveAs
        summary_histos = []
        for idx, n_crossings in enumerate(order):
            data = results[n_crossings]
            h = data["histogram"].Clone("TotalSE_per_event_n{}_clone".format(n_crossings))
            h.SetDirectory(0)
            summary_histos.append(h)
            col = colors[idx % len(colors)]
            h.SetLineColor(col)
            h.SetLineWidth(2)
            h.SetFillStyle(0)  # no fill so overlaid histograms don't cover each other
            h.Draw("HIST SAME")
            line = leg_box.AddText("N=%d: #mu=%.1f, eff=%.1f%%" % (
                n_crossings, data["mean"], data["efficiency"] * 100.0))
            line.SetTextColor(col)
        leg_box.Draw()
        c_summary.Update()
        summary_pdf = os.path.join(plots_dir, "TotalSE_per_event_summary.pdf")
        summary_root = os.path.join(plots_dir, "TotalSE_per_event_summary.root")
        c_summary.SaveAs(summary_pdf)
        c_summary.SaveAs(summary_root)
        print("Summary plot saved: {} and {}".format(summary_pdf, summary_root))

        # Prepare arrays for efficiency and mean SE vs crossings plots
        import array
        crossings_sorted = sorted(results.keys())
        n_points = len(crossings_sorted)
        arr_crossings = array.array('d', [float(n) for n in crossings_sorted])
        arr_efficiency = array.array('d', [results[n]["efficiency"] * 100.0 for n in crossings_sorted])
        arr_mean = array.array('d', [results[n]["mean"] for n in crossings_sorted])
        arr_std = array.array('d', [results[n]["std"] for n in crossings_sorted])
        arr_zeros = array.array('d', [0.0] * n_points)  # no x error

        # Plot 1: Efficiency vs crossings
        c_eff = TCanvas("c_efficiency", "Efficiency vs crossings", 800, 600)
        c_eff.SetGrid()
        gr_eff = TGraph(n_points, arr_crossings, arr_efficiency)
        gr_eff.SetTitle("Efficiency vs shell crossings (%d events);Number of shell crossings;Efficiency (%%)" % n_events)
        gr_eff.SetMarkerStyle(20)
        gr_eff.SetMarkerSize(1.2)
        gr_eff.SetMarkerColor(4)
        gr_eff.SetLineColor(4)
        gr_eff.SetLineWidth(2)
        gr_eff.GetYaxis().SetRangeUser(0, 105)
        gr_eff.Draw("APL")
        c_eff.Update()
        eff_pdf = os.path.join(plots_dir, "Efficiency_vs_crossings.pdf")
        eff_root = os.path.join(plots_dir, "Efficiency_vs_crossings.root")
        c_eff.SaveAs(eff_pdf)
        c_eff.SaveAs(eff_root)
        print("Efficiency plot saved: {} and {}".format(eff_pdf, eff_root))

        # Plot 2: Mean SE vs crossings with std dev as error bars
        c_mean = TCanvas("c_mean_se", "Mean SE vs crossings", 800, 600)
        c_mean.SetGrid()
        gr_mean = TGraphErrors(n_points, arr_crossings, arr_mean, arr_zeros, arr_std)
        gr_mean.SetTitle("Mean SE vs shell crossings (%d events);Number of shell crossings;Mean secondary electrons per event" % n_events)
        gr_mean.SetMarkerStyle(20)
        gr_mean.SetMarkerSize(1.2)
        gr_mean.SetMarkerColor(2)
        gr_mean.SetLineColor(2)
        gr_mean.SetLineWidth(2)
        gr_mean.Draw("APL")
        c_mean.Update()
        mean_pdf = os.path.join(plots_dir, "MeanSE_vs_crossings.pdf")
        mean_root = os.path.join(plots_dir, "MeanSE_vs_crossings.root")
        c_mean.SaveAs(mean_pdf)
        c_mean.SaveAs(mean_root)
        print("Mean SE plot saved: {} and {}".format(mean_pdf, mean_root))

    print("Results written to {}".format(out_dir))
    print("ROOT file: {}".format(out_root))
    print("Summary:\n{}".format(summary))
    return results


def main():
    import argparse
    p = argparse.ArgumentParser(description="Run toy model: events with many crossings (config-driven)")
    p.add_argument("config_file", nargs="?", default="config/toy_model/toy_model_config.json",
                   help="Path to toy model config JSON (default: config/toy_model/toy_model_config.json)")
    p.add_argument("--n-events", type=int, default=None, help="Override n_events from config")
    p.add_argument("--seed", type=int, default=None, help="Override seed from config")
    args = p.parse_args()
    overrides = {}
    if args.n_events is not None:
        overrides["n_events"] = args.n_events
    if args.seed is not None:
        overrides["seed"] = args.seed
    run_toy_events(args.config_file, cli_overrides=overrides if overrides else None)


if __name__ == "__main__":
    main()
