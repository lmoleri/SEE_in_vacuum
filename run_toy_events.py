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
    # Standalone fallback
    def sample_poisson(mu, random_gen=None):
        if random_gen is None:
            import ROOT
            random_gen = ROOT.gRandom
        if mu <= 0:
            return 0
        if mu > 100:
            sample = random_gen.Gaus(mu, np.sqrt(mu))
            return max(0, int(round(sample)))
        u = random_gen.Uniform(0.0, 1.0)
        k = 0
        cumulative = exp(-mu)
        term = exp(-mu)
        while cumulative < u:
            k += 1
            term = term * mu / k
            cumulative += term
            if k > 10000 or term < 1e-100:
                break
        return k


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
    from ROOT import TFile, TH1D, TCanvas, TLatex, gStyle

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
        }

    summary = "n_events={} crossings_scan={} epsilon={} B={} alpha={} depth={} P_esc={}\n".format(
        n_events, crossings_scan, epsilon, B, alpha, depth, P_esc)
    for n_crossings, data in results.items():
        summary += "crossings={} mean_SE={:.4f} std_SE={:.4f} min={} max={}\n".format(
            n_crossings, data["mean"], data["std"], data["min"], data["max"])

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
    for n_crossings, data in results.items():
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
        c.Update()
        pdf_path = os.path.join(plots_dir, "TotalSE_per_event_{}_crossings.pdf".format(n_crossings))
        c.SaveAs(pdf_path)
        print("Plot saved: {}".format(pdf_path))

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
