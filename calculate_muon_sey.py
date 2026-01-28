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
        'n_events': n_events
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
    
    # Create output histograms
    max_SE = max(per_event_SE) if per_event_SE else 0
    n_bins = max(50, max_SE + 10)  # At least 50 bins, or enough to cover max
    
    sey_hist = TH1D("SEY_MonteCarlo", "Secondary Electron Yield (Monte Carlo)", 
                    n_bins, -0.5, max_SE + 0.5)
    sey_hist.SetXTitle("Number of Secondary Electrons per Event")
    sey_hist.SetYTitle("Number of Events")
    sey_hist.SetFillColor(ROOT.kBlue)
    sey_hist.SetFillStyle(3004)
    sey_hist.SetLineColor(ROOT.kBlue)
    sey_hist.SetLineWidth(2)
    
    # Fill histogram
    for n_se in per_event_SE:
        sey_hist.Fill(n_se)
    
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
"""
    
    print(summary_text)
    
    # Write summary to file
    summary_hist = TH1D("Summary", "Summary", 1, 0, 1)
    summary_hist.SetTitle(summary_text)
    summary_hist.Write()
    
    output_file.Close()
    root_file.Close()
    
    # Create plot
    create_plot(sey_hist, statistics, root_file_path)
    
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
    gStyle.SetOptStat(1110)  # Show mean, RMS, entries
    gStyle.SetOptFit(0)
    
    canvas = TCanvas("SEYCanvas", "Secondary Electron Yield (Monte Carlo)", 1000, 700)
    canvas.SetGrid()
    
    sey_hist.Draw("hist")
    
    # Add text with statistics
    text = TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.SetTextAlign(12)
    
    y_pos = 0.95
    text.DrawLatex(0.15, y_pos, "Monte Carlo SEY Calculation")
    y_pos -= 0.04
    text.DrawLatex(0.15, y_pos, f"Mean SEY: {statistics['mean_SE']:.4f} #pm {statistics['std_SE']:.4f}")
    y_pos -= 0.04
    text.DrawLatex(0.15, y_pos, f"Expected (theoretical): {statistics['expected_mean']:.4f}")
    y_pos -= 0.04
    text.DrawLatex(0.15, y_pos, f"Total events: {statistics['n_events']}")
    y_pos -= 0.04
    text.DrawLatex(0.15, y_pos, f"Total SE: {statistics['total_SE']}")
    y_pos -= 0.04
    text.DrawLatex(0.15, y_pos, f"P_{{esc}} = {P_ESC:.4f} (z = {Z_DEPTH/10:.1f} nm)")
    
    canvas.Update()
    
    # Save plot
    plot_file_path = input_file_path.replace(".root", "_SEY_MonteCarlo.pdf")
    canvas.SaveAs(plot_file_path)
    print(f"\nPlot saved to: {plot_file_path}")
    
    # Also save as ROOT file
    root_plot_path = input_file_path.replace(".root", "_SEY_MonteCarlo_plot.root")
    canvas.SaveAs(root_plot_path)
    print(f"Plot (ROOT) saved to: {root_plot_path}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Calculate Secondary Electron Yield from Muon Interactions using Monte Carlo method",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python calculate_muon_sey.py results/scan_mu4GeV/SEE_in_vacuum_thick5nm_particlemu-_energy4000MeV_events10000.root
  
  python calculate_muon_sey.py input.root --histogram EdepPrimary --seed 12345
        """
    )
    
    parser.add_argument("input_file", 
                       help="Input ROOT file from Geant4 simulation")
    parser.add_argument("--histogram", "-H",
                       default="EdepPrimary",
                       help="Name of energy deposition histogram (default: EdepPrimary)")
    parser.add_argument("--seed", "-s",
                       type=int,
                       default=42,
                       help="Random number generator seed (default: 42)")
    parser.add_argument("--epsilon", "-e",
                       type=float,
                       default=EPSILON,
                       help=f"Energy per free electron in eV (default: {EPSILON})")
    parser.add_argument("--B", "-b",
                       type=float,
                       default=B,
                       help=f"Surface escape probability (default: {B})")
    parser.add_argument("--alpha", "-a",
                       type=float,
                       default=ALPHA,
                       help=f"Attenuation coefficient in Å^-1 (default: {ALPHA})")
    parser.add_argument("--depth", "-d",
                       type=float,
                       default=Z_DEPTH,
                       help=f"Production depth in Å (default: {Z_DEPTH}, i.e., {Z_DEPTH/10} nm)")
    parser.add_argument("--bin-by-bin", "-b",
                       action="store_true",
                       help="Process histogram bin-by-bin instead of sampling (default: sampling)")
    
    args = parser.parse_args()
    
    # Update global constants if provided
    global EPSILON, B, ALPHA, Z_DEPTH, P_ESC
    EPSILON = args.epsilon
    B = args.B
    ALPHA = args.alpha
    Z_DEPTH = args.depth
    P_ESC = B * exp(-ALPHA * Z_DEPTH)
    
    print("=" * 70)
    print("Monte Carlo SEY Calculation for Muon Interactions")
    print("=" * 70)
    print()
    
    # Process file
    use_sampling = not args.bin_by_bin
    result = process_histogram_from_file(
        args.input_file, args.histogram, args.seed, use_sampling)
    
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
