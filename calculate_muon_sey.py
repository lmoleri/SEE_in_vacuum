#!/usr/bin/env python3
"""
Monte Carlo calculation for Secondary Electron Emission from Muon Interactions

This script reads energy deposition data from Geant4 ROOT output files and
calculates the probability of extracting secondary electrons using a custom
Monte Carlo method.

Usage:
    python calculate_muon_sey.py <input_root_file> [options]

Example:
    python calculate_muon_sey.py results/scan_mu4GeV/SEE_in_vacuum_thick5nm_particlemu-_energy4000MeV_events10000.root
"""

import sys
import argparse
import json
import os
import numpy as np
from math import exp, factorial
import ROOT
from ROOT import TFile, TH1D, TCanvas, TLegend, TLatex, gStyle, gRandom

# Physical constants
EPSILON = 27.0  # eV - average energy per internal free electron production
B = 0.46  # surface escape probability
ALPHA = 0.0075  # Å^-1 - attenuation coefficient
Z_DEPTH = 25.0  # Å - production depth (2.5 nm = 25 Å)

# Calculate escape probability at depth z
P_ESC = B * exp(-ALPHA * Z_DEPTH)


def poisson_cdf(k, mu):
    """
    Calculate cumulative distribution function for Poisson distribution.
    
    Parameters:
    -----------
    k : int
        Number of events
    mu : float
        Poisson parameter (mean)
    
    Returns:
    --------
    float
        Cumulative probability P(X <= k)
    """
    if mu < 0:
        return 0.0
    if k < 0:
        return 0.0
    
    cdf = 0.0
    for n in range(k + 1):
        if n == 0:
            term = exp(-mu)
        else:
            term = (mu ** n) / factorial(n) * exp(-mu)
        cdf += term
    
    return cdf


def sample_poisson(mu, random_gen=None):
    """
    Sample from Poisson distribution using inverse CDF method.
    
    Parameters:
    -----------
    mu : float
        Poisson parameter (mean)
    random_gen : TRandom3 or None
        ROOT random number generator (if None, uses gRandom)
    
    Returns:
    --------
    int
        Sampled value k from Poisson(mu)
    """
    if random_gen is None:
        random_gen = gRandom
    
    if mu <= 0:
        return 0
    
    # For very large mu, use normal approximation
    if mu > 100:
        # Normal approximation: N(mu, sqrt(mu))
        sample = random_gen.Gaus(mu, np.sqrt(mu))
        return max(0, int(round(sample)))
    
    # Generate uniform random number
    u = random_gen.Uniform(0.0, 1.0)
    
    # Find smallest k such that CDF(k) >= u
    k = 0
    cumulative = exp(-mu)  # P(k=0)
    term = exp(-mu)  # Current term P(k)
    
    while cumulative < u:
        k += 1
        # Calculate P(k) = (mu^k / k!) * exp(-mu)
        # Use recurrence: P(k) = P(k-1) * mu / k
        term = term * mu / k
        cumulative += term
        
        # Safety check for numerical issues
        if k > 10000:
            print(f"Warning: Poisson sampling reached k={k} for mu={mu:.2f}, u={u:.6f}")
            break
        # Check for numerical underflow
        if term < 1e-100:
            # Term is negligible, we've converged
            break
    
    return k


def calculate_sey_monte_carlo(edep_hist, random_gen=None, use_histogram_sampling=True):
    """
    Calculate secondary electron yield using Monte Carlo method.
    
    Parameters:
    -----------
    edep_hist : TH1D
        Histogram containing energy deposition per event
    random_gen : TRandom3 or None
        ROOT random number generator
    use_histogram_sampling : bool
        If True, sample from histogram distribution (for aggregated data).
        If False, process each bin as representing actual events.
    
    Returns:
    --------
    tuple : (total_SE, per_event_SE, statistics_dict)
        - total_SE: total number of secondary electrons
        - per_event_SE: list of SE counts per event
        - statistics_dict: dictionary with mean, std, etc.
    """
    if random_gen is None:
        random_gen = gRandom
    
    n_events = int(edep_hist.GetEntries())
    per_event_SE = []
    total_SE = 0
    
    print(f"Processing {n_events} events...")
    print(f"Physical parameters:")
    print(f"  ε (energy per free electron) = {EPSILON:.1f} eV")
    print(f"  B (surface escape prob) = {B:.2f}")
    print(f"  α (attenuation coeff) = {ALPHA:.4f} Å^-1")
    print(f"  z (production depth) = {Z_DEPTH:.1f} Å ({Z_DEPTH/10:.1f} nm)")
    print(f"  P_esc(z) = {P_ESC:.4f}")
    print()
    
    if use_histogram_sampling:
        # Method 1: Sample energy deposition from histogram distribution
        # This treats the histogram as a probability distribution
        print("Using histogram sampling method...")
        for i in range(n_events):
            # Sample energy deposition from the histogram distribution
            delta_E = edep_hist.GetRandom()  # Sample energy deposition in eV
            
            if delta_E <= 0:
                per_event_SE.append(0)
                continue
            
            # Calculate average number of free electrons
            N_int = delta_E / EPSILON
            
            # Calculate Poisson parameter
            mu = N_int * P_ESC
            
            # Sample from Poisson distribution
            N_SE_i = sample_poisson(mu, random_gen)
            
            per_event_SE.append(N_SE_i)
            total_SE += N_SE_i
            
            # Progress indicator
            if (i + 1) % 1000 == 0:
                print(f"  Processed {i+1}/{n_events} events...", end='\r')
    else:
        # Method 2: Process each bin as representing actual events
        # This iterates through histogram bins and processes events based on bin content
        print("Using bin-by-bin processing method...")
        n_bins = edep_hist.GetNbinsX()
        events_processed = 0
        
        for bin_i in range(1, n_bins + 1):  # ROOT bins are 1-indexed
            bin_content = edep_hist.GetBinContent(bin_i)
            if bin_content <= 0:
                continue
            
            # Get bin center (energy deposition value)
            delta_E = edep_hist.GetBinCenter(bin_i)  # in eV
            
            if delta_E <= 0:
                continue
            
            # Process each event in this bin
            n_events_in_bin = int(bin_content)
            for _ in range(n_events_in_bin):
                # Calculate average number of free electrons
                N_int = delta_E / EPSILON
                
                # Calculate Poisson parameter
                mu = N_int * P_ESC
                
                # Sample from Poisson distribution
                N_SE_i = sample_poisson(mu, random_gen)
                
                per_event_SE.append(N_SE_i)
                total_SE += N_SE_i
                events_processed += 1
                
                # Progress indicator
                if events_processed % 1000 == 0:
                    print(f"  Processed {events_processed}/{n_events} events...", end='\r')
        
        print(f"\nCompleted processing {events_processed} events from histogram bins")
    
    # Calculate statistics
    per_event_SE_array = np.array(per_event_SE)
    mean_SE = np.mean(per_event_SE_array)
    std_SE = np.std(per_event_SE_array)
    median_SE = np.median(per_event_SE_array)
    
    # Calculate fraction of events with at least one secondary electron
    n_events_with_SE = np.sum(per_event_SE_array > 0)
    fraction_with_SE = n_events_with_SE / n_events if n_events > 0 else 0.0
    
    # Calculate expected mean (for validation)
    # Expected mean = average of mu over all events
    # We approximate this by: mean(ΔE) / ε * P_esc
    mean_edep = edep_hist.GetMean()
    expected_mean = (mean_edep / EPSILON) * P_ESC
    
    statistics = {
        'total_SE': total_SE,
        'mean_SE': mean_SE,
        'std_SE': std_SE,
        'median_SE': median_SE,
        'min_SE': int(np.min(per_event_SE_array)),
        'max_SE': int(np.max(per_event_SE_array)),
        'expected_mean': expected_mean,
        'mean_edep': mean_edep,
        'n_events': n_events,
        'n_events_with_SE': n_events_with_SE,
        'fraction_with_SE': fraction_with_SE
    }
    
    return total_SE, per_event_SE, statistics


def process_histogram_from_file(root_file_path, hist_name="EdepPrimary", seed=42, use_histogram_sampling=True):
    """
    Process energy deposition histogram from ROOT file.
    
    Parameters:
    -----------
    root_file_path : str
        Path to input ROOT file
    hist_name : str
        Name of the histogram to read
    seed : int
        Random number generator seed
    
    Returns:
    --------
    tuple : (total_SE, per_event_SE, statistics, output_file_path)
    """
    # Open ROOT file
    root_file = TFile.Open(root_file_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Error: Cannot open ROOT file {root_file_path}")
        return None, None, None, None
    
    print(f"Opened ROOT file: {root_file_path}")
    
    # Get histogram
    edep_hist = root_file.Get(hist_name)
    if not edep_hist:
        print(f"Error: Histogram '{hist_name}' not found in file")
        print("Available objects:")
        root_file.ls()
        root_file.Close()
        return None, None, None, None
    
    print(f"Found histogram: {hist_name}")
    print(f"  Entries: {edep_hist.GetEntries()}")
    print(f"  Mean: {edep_hist.GetMean():.2f} eV")
    print(f"  RMS: {edep_hist.GetRMS():.2f} eV")
    print()
    
    # Initialize random number generator
    random_gen = ROOT.TRandom3(seed)
    gRandom.SetSeed(seed)
    
    # Calculate SEY
    total_SE, per_event_SE, statistics = calculate_sey_monte_carlo(
        edep_hist, random_gen, use_histogram_sampling)
    
    # Create output histograms: exactly one bin per integer (0, 1, 2, ...)
    max_SE = max(per_event_SE) if per_event_SE else 0
    n_bins = int(max_SE) + 1  # bins for 0, 1, 2, ..., max_SE
    xmin = -0.5
    xmax = max_SE + 0.5
    
    sey_hist = TH1D("SEY_MonteCarlo", "Secondary Electron Yield (Monte Carlo)", 
                    n_bins, xmin, xmax)
    sey_hist.SetXTitle("Number of Secondary Electrons per Event")
    sey_hist.SetYTitle("Number of Events")
    sey_hist.SetFillColor(9)  # kBlue = 9
    sey_hist.SetFillStyle(3004)
    sey_hist.SetLineColor(9)  # kBlue = 9
    sey_hist.SetLineWidth(2)
    sey_hist.SetMarkerStyle(20)
    sey_hist.SetMarkerColor(9)  # kBlue = 9
    
    # Fill histogram
    for n_se in per_event_SE:
        sey_hist.Fill(n_se)
    
    # Ensure histogram has proper range
    if sey_hist.GetEntries() > 0:
        sey_hist.SetMinimum(0)  # Start Y-axis at 0
    
    # Create plot FIRST (before writing to file, to keep histogram valid)
    try:
        create_plot(sey_hist, statistics, root_file_path)
    except Exception as e:
        print("Warning: Could not create plot: {}".format(e))
        import traceback
        traceback.print_exc()
    
    # Create output file
    output_file_path = root_file_path.replace(".root", "_SEY_MonteCarlo.root")
    output_file = TFile.Open(output_file_path, "RECREATE")
    
    # Write histograms
    sey_hist.Write()
    edep_hist.Write("EdepPrimary_original")  # Save original for reference
    
    # Create summary text
    summary_text = f"""
Monte Carlo SEY Calculation Summary
===================================
Input file: {root_file_path}
Histogram: {hist_name}

Physical Parameters:
  ε (energy per free electron) = {EPSILON:.1f} eV
  B (surface escape prob) = {B:.2f}
  α (attenuation coeff) = {ALPHA:.4f} Å^-1
  z (production depth) = {Z_DEPTH:.1f} Å ({Z_DEPTH/10:.1f} nm)
  P_esc(z) = {P_ESC:.4f}

Results:
  Total events processed: {statistics['n_events']}
  Total secondary electrons: {statistics['total_SE']}
  Mean SEY per event: {statistics['mean_SE']:.4f}
  Expected mean (theoretical): {statistics['expected_mean']:.4f}
  Standard deviation: {statistics['std_SE']:.4f}
  Median SEY: {statistics['median_SE']:.1f}
  Min SEY: {statistics['min_SE']}
  Max SEY: {statistics['max_SE']}
  
  Mean energy deposition: {statistics['mean_edep']:.2f} eV
  Mean free electrons per event: {statistics['mean_edep']/EPSILON:.2f}
  
  Events with at least one SE: {statistics['n_events_with_SE']}/{statistics['n_events']} ({statistics['fraction_with_SE']*100:.2f}%)
"""
    
    print(summary_text)
    
    # Write summary to file
    summary_hist = TH1D("Summary", "Summary", 1, 0, 1)
    summary_hist.SetTitle(summary_text)
    summary_hist.Write()
    
    # Close files
    output_file.Close()
    root_file.Close()
    
    return total_SE, per_event_SE, statistics, output_file_path


def create_plot(sey_hist, statistics, input_file_path):
    """
    Create and save plot of SEY distribution.
    
    Parameters:
    -----------
    sey_hist : TH1D
        Histogram of SEY values
    statistics : dict
        Statistics dictionary
    input_file_path : str
        Input file path (for output naming)
    """
    gStyle.SetOptStat(1110)  # Show mean, RMS, entries (top-right)
    gStyle.SetOptFit(0)
    
    canvas = TCanvas("SEYCanvas", "Secondary Electron Yield (Monte Carlo)", 1000, 700)
    canvas.SetGrid()
    canvas.cd()
    
    # Set Y-axis to logarithmic scale
    canvas.SetLogy(1)
    
    # Ensure histogram is properly styled (in case it wasn't set earlier)
    # Wrap in try-except to avoid crashes if color setting fails
    try:
        sey_hist.SetFillColor(9)  # kBlue = 9
        sey_hist.SetFillStyle(3004)
        sey_hist.SetLineColor(9)  # kBlue = 9
        sey_hist.SetMarkerColor(9)  # kBlue = 9
    except:
        pass  # Skip color setting if it fails
    
    # These should always work
    sey_hist.SetLineWidth(2)
    sey_hist.SetMarkerStyle(20)
    
    # For log scale, set minimum to a small positive value (not 0)
    min_y = sey_hist.GetMinimum()
    if min_y <= 0:
        sey_hist.SetMinimum(0.1)  # Small positive value for log scale
    
    # Set reasonable Y-axis range for log scale
    max_y = sey_hist.GetMaximum()
    if max_y > 0:
        sey_hist.SetMaximum(max_y * 10)  # Larger range for log scale
    else:
        sey_hist.SetMaximum(1000)  # Default if empty
    
    # Draw histogram - use "HIST" option for better visibility
    sey_hist.Draw("HIST")
    canvas.Update()
    
    # Verify histogram is visible
    if sey_hist.GetEntries() == 0:
        print("Warning: Histogram has no entries!")
    else:
        print("Histogram has {} entries, max bin: {}".format(
            sey_hist.GetEntries(), sey_hist.GetMaximum()))
    
    # Determine output directory in plots folder
    # Extract run identifier from input file path
    import os
    input_dir = os.path.dirname(input_file_path)
    input_basename = os.path.basename(input_file_path)
    
    # Create plots subdirectory based on input directory structure
    # e.g., results/scan_.../file.root -> plots/scan_.../
    if "results" in input_dir:
        plots_subdir = input_dir.replace("results", "plots")
    else:
        # Fallback: use input directory name
        plots_subdir = os.path.join("plots", os.path.basename(input_dir))
    
    # Create plots directory if it doesn't exist
    os.makedirs(plots_subdir, exist_ok=True)
    
    # Generate plot filename
    plot_basename = input_basename.replace(".root", "_SEY_MonteCarlo")
    plot_file_path = os.path.join(plots_subdir, plot_basename + ".pdf")
    
    # Save plot first (before adding text that might cause issues)
    canvas.SaveAs(plot_file_path)
    print("\nPlot saved to: {}".format(plot_file_path))
    
    # Try to add text, but don't fail if it doesn't work
    try:
        text = TLatex()
        text.SetNDC(True)
        text.SetTextSize(0.025)
        text.SetTextAlign(12)
        
        # Position text on the LEFT side to avoid overlap with:
        # - Title (top-center, y ~ 0.95)
        # - Stat box (top-right, x ~ 0.6-0.95, y ~ 0.7-0.95)
        x_pos = 0.15  # Left side
        y_pos = 0.85  # Below title area, above stat box
        
        # Draw text one by one with explicit string conversion
        text.DrawLatex(float(x_pos), float(y_pos), "Monte Carlo SEY Calculation")
        y_pos -= 0.035
        mean_str = "Mean SEY: %.4f" % statistics['mean_SE']
        text.DrawLatex(float(x_pos), float(y_pos), mean_str)
        y_pos -= 0.035
        # Expected = theoretical mean SEY = (mean_ΔE / ε) * P_esc
        exp_str = "Expected (theoretical): %.4f" % statistics['expected_mean']
        text.DrawLatex(float(x_pos), float(y_pos), exp_str)
        y_pos -= 0.035
        events_str = "Events: %d" % statistics['n_events']
        text.DrawLatex(float(x_pos), float(y_pos), events_str)
        y_pos -= 0.035
        total_str = "Total SE: %d" % statistics['total_SE']
        text.DrawLatex(float(x_pos), float(y_pos), total_str)
        y_pos -= 0.035
        # Fraction of events with at least one SE
        frac_str = "Events with SE: %.2f%% (%d/%d)" % (
            statistics['fraction_with_SE'] * 100.0,
            statistics['n_events_with_SE'],
            statistics['n_events']
        )
        text.DrawLatex(float(x_pos), float(y_pos), frac_str)
        y_pos -= 0.035
        # P_esc = escape probability at production depth z
        # P_esc = B * exp(-alpha * z)
        pesc_str = "P_{esc}(z) = %.4f" % P_ESC
        text.DrawLatex(float(x_pos), float(y_pos), pesc_str)
        y_pos -= 0.035
        depth_str = "Production depth: z = %.1f nm" % (Z_DEPTH/10)
        text.DrawLatex(float(x_pos), float(y_pos), depth_str)
        
        canvas.Update()
        # Save again with text
        canvas.SaveAs(plot_file_path)
    except Exception as e:
        print("Note: Text annotations skipped (histogram saved successfully)")
        # PDF is already saved, so we're good
    
    # Also save as ROOT file in plots folder
    root_plot_path = os.path.join(plots_subdir, plot_basename + "_plot.root")
    canvas.SaveAs(root_plot_path)
    print("Plot (ROOT) saved to: {}".format(root_plot_path))


def load_toy_config(config_path):
    """Load toy model config (JSON). Returns dict with keys used by calculate_muon_sey."""
    with open(config_path, "r") as f:
        cfg = json.load(f)
    # Input file: edep_root_file (run_toy_events) or input_file
    raw_input = cfg.get("edep_root_file") or cfg.get("input_file")
    return {
        "input_file": raw_input,
        "histogram": cfg.get("histogram", "EdepPrimary"),
        "seed": int(cfg.get("seed", 42)),
        "epsilon": float(cfg.get("epsilon", EPSILON)),
        "B": float(cfg.get("B", B)),
        "alpha": float(cfg.get("alpha", ALPHA)),
        "depth": float(cfg.get("depth", Z_DEPTH)),
        "bin_by_bin": bool(cfg.get("bin_by_bin", False)),
    }


def main():
    """Main function."""
    # Declare globals at the start to avoid syntax errors
    global EPSILON, B, ALPHA, Z_DEPTH, P_ESC
    
    # Store original values for help strings
    epsilon_default = EPSILON
    B_default = B
    alpha_default = ALPHA
    depth_default = Z_DEPTH
    
    parser = argparse.ArgumentParser(
        description="Calculate Secondary Electron Yield from Muon Interactions using Monte Carlo method",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python calculate_muon_sey.py results/scan_mu4GeV/SEE_in_vacuum_thick5nm_particlemu-_energy4000MeV_events10000.root
  python calculate_muon_sey.py --config config/toy_model/toy_model_config.json
  python calculate_muon_sey.py input.root --histogram EdepPrimary --seed 12345
        """
    )
    
    parser.add_argument("input_file", nargs="?", default=None,
                       help="Input ROOT file from Geant4 simulation (optional if --config is used)")
    parser.add_argument("--config", "-c", default=None,
                       help="Path to toy model config JSON (provides input_file, histogram, seed, epsilon, B, alpha, depth, bin_by_bin)")
    parser.add_argument("--histogram", "-H",
                       default=None,
                       help="Name of energy deposition histogram (default: EdepPrimary or from config)")
    parser.add_argument("--seed", "-s",
                       type=int,
                       default=None,
                       help="Random number generator seed (default: 42 or from config)")
    parser.add_argument("--epsilon", "-e",
                       type=float,
                       default=None,
                       help=f"Energy per free electron in eV (default: {epsilon_default} or from config)")
    parser.add_argument("--B", 
                       type=float,
                       default=None,
                       help=f"Surface escape probability (default: {B_default} or from config)")
    parser.add_argument("--alpha", "-a",
                       type=float,
                       default=None,
                       help=f"Attenuation coefficient in Å^-1 (default: {alpha_default} or from config)")
    parser.add_argument("--depth", "-d",
                       type=float,
                       default=None,
                       help=f"Production depth in Å (default: {depth_default}, i.e., {depth_default/10} nm, or from config)")
    parser.add_argument("--bin-by-bin",
                       action="store_true",
                       help="Process histogram bin-by-bin instead of sampling (overrides config if set)")
    
    args = parser.parse_args()
    
    # Load config if provided
    cfg = None
    config_path = args.config
    if config_path:
        if not os.path.isfile(config_path) and not os.path.dirname(config_path):
            for alt in (os.path.join("config", "toy_model", config_path), os.path.join("config", config_path)):
                if os.path.isfile(alt):
                    config_path = alt
                    break
        if not os.path.isfile(config_path):
            print("Error: Config file not found: {}".format(args.config))
            sys.exit(1)
        args.config = config_path
        cfg = load_toy_config(config_path)
    
    # Resolve input file: CLI > config
    input_file = args.input_file
    if input_file is None and cfg is not None:
        input_file = cfg.get("input_file")
    if not input_file:
        print("Error: No input ROOT file. Provide input_file as positional argument or in --config (edep_root_file / input_file).")
        sys.exit(1)
    # If path is relative and we have a config file, resolve relative to project root (config may be in config/toy_model/ or config/geant4/)
    if not os.path.isabs(input_file) and args.config and os.path.isfile(args.config):
        base = os.path.dirname(os.path.abspath(args.config))
        while base and os.path.basename(base) in ("toy_model", "geant4", "config"):
            base = os.path.dirname(base)
        input_file = os.path.normpath(os.path.join(base, input_file))
    
    # Resolve other options: CLI overrides config, config overrides defaults
    histogram = args.histogram if args.histogram is not None else (cfg.get("histogram") if cfg else "EdepPrimary")
    seed = args.seed if args.seed is not None else (cfg.get("seed") if cfg else 42)
    epsilon_val = args.epsilon if args.epsilon is not None else (cfg.get("epsilon") if cfg else EPSILON)
    B_val = args.B if args.B is not None else (cfg.get("B") if cfg else B)
    alpha_val = args.alpha if args.alpha is not None else (cfg.get("alpha") if cfg else ALPHA)
    depth_val = args.depth if args.depth is not None else (cfg.get("depth") if cfg else Z_DEPTH)
    bin_by_bin = args.bin_by_bin or (cfg.get("bin_by_bin") if cfg else False)
    
    p_esc_val = B_val * exp(-alpha_val * depth_val)
    
    # Update module-level constants for use in functions
    EPSILON = epsilon_val
    B = B_val
    ALPHA = alpha_val
    Z_DEPTH = depth_val
    P_ESC = p_esc_val
    
    print("=" * 70)
    print("Monte Carlo SEY Calculation for Muon Interactions")
    print("=" * 70)
    if args.config:
        print("Config file: {}".format(args.config))
    print()
    
    # Process file
    use_sampling = not bin_by_bin
    result = process_histogram_from_file(
        input_file, histogram, seed, use_sampling)
    
    if result[0] is not None:
        total_SE, per_event_SE, statistics, output_file = result
        print()
        print("=" * 70)
        print(f"Results saved to: {output_file}")
        print("=" * 70)
    else:
        print("Error: Failed to process file")
        sys.exit(1)


if __name__ == "__main__":
    main()
