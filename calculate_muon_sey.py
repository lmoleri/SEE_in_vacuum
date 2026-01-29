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
from ROOT import TFile, TH1D, TCanvas, TLegend, TLatex, TPaveText, gStyle, gRandom

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
    Sample from Poisson(μ) using ROOT's TRandom::Poisson.
    
    Parameters:
    -----------
    mu : float
        Poisson parameter (mean)
    random_gen : TRandom or None
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
    return int(random_gen.Poisson(mu))


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
    
    # Upper edge of bin containing 0: "Edep > 0" means above this (match histogram definition)
    bin_zero = edep_hist.GetXaxis().FindBin(0.0)
    if 1 <= bin_zero <= edep_hist.GetNbinsX():
        edep_positive_threshold = edep_hist.GetXaxis().GetBinUpEdge(bin_zero)
    else:
        edep_positive_threshold = 0.0
    
    print(f"Processing {n_events} events...")
    print(f"Physical parameters:")
    print(f"  ε (energy per free electron) = {EPSILON:.1f} eV")
    print(f"  B (surface escape prob) = {B:.2f}")
    print(f"  α (attenuation coeff) = {ALPHA:.4f} Å^-1")
    print(f"  z (production depth) = {Z_DEPTH:.1f} Å ({Z_DEPTH/10:.1f} nm)")
    print(f"  P_esc(z) = {P_ESC:.4f}")
    print()
    
    # Track sampled Edep when Edep > threshold (for comparison with histogram)
    sampled_edep_when_positive = []
    n_MC_edep_positive = 0   # MC count with Edep > threshold (match histogram "Edep > 0")
    n_MC_edep_positive_with_SE = 0
    sampled_edep_debug = []  # All sampled Edep values (for debug plot, histogram sampling only)
    
    if use_histogram_sampling:
        # Method 1: Sample energy deposition from histogram distribution
        # This treats the histogram as a probability distribution
        print("Using histogram sampling method...")
        for i in range(n_events):
            # Sample energy deposition from the histogram distribution
            delta_E = edep_hist.GetRandom()  # Sample energy deposition in eV
            sampled_edep_debug.append(delta_E)
            
            if delta_E <= 0:
                per_event_SE.append(0)
                continue
            
            # Count "Edep > 0" only when above zero bin (match histogram definition)
            if delta_E > edep_positive_threshold:
                sampled_edep_when_positive.append(delta_E)
                n_MC_edep_positive += 1
            # Poisson parameter: μ = (ΔE/ε) * P_esc (must apply P_esc)
            mu = (delta_E / EPSILON) * P_ESC
            # Sample from Poisson distribution (ROOT's Poisson; custom sample_poisson had a bug)
            N_SE_i = int(random_gen.Poisson(mu)) if mu > 0 else 0
            
            if delta_E > edep_positive_threshold and N_SE_i >= 1:
                n_MC_edep_positive_with_SE += 1
            
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
                # Sample from Poisson distribution (use ROOT's Poisson)
                N_SE_i = int(random_gen.Poisson(mu)) if mu > 0 else 0
                
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
    
    # Fraction of events with Edep > 0: exclude the bin that contains 0 (not bins with center > 0)
    integral_total = edep_hist.Integral()
    bin_containing_zero = edep_hist.GetXaxis().FindBin(0.0)
    if 1 <= bin_containing_zero <= edep_hist.GetNbinsX():
        content_zero_bin = edep_hist.GetBinContent(bin_containing_zero)
        integral_positive = integral_total - content_zero_bin
    else:
        integral_positive = integral_total
    fraction_with_edep = (integral_positive / integral_total) if integral_total > 0 else 0.0
    
    # From histogram: mean Edep given Edep > 0, and theoretical P(≥1 SE | Edep > 0)
    # "Edep > 0" = exclude the bin that contains 0 (same as integral_positive above).
    # For each bin i (center E_i, content n_i): P(≥1 SE | E_i) = 1 - exp(-(E_i/ε)*P_esc).
    # Theoretical P(≥1 SE | Edep>0) = Σ (n_i * P(≥1 SE | E_i)) / integral_positive  (weighted by bin content).
    sum_edep_positive = 0.0
    sum_p_at_least_one_SE = 0.0
    for bin_i in range(1, edep_hist.GetNbinsX() + 1):
        if bin_i == bin_containing_zero:
            continue
        center = edep_hist.GetBinCenter(bin_i)
        content = edep_hist.GetBinContent(bin_i)
        if content <= 0:
            continue
        sum_edep_positive += center * content
        mu_bin = (center / EPSILON) * P_ESC
        p_at_least_one = 1.0 - exp(-mu_bin)
        sum_p_at_least_one_SE += content * p_at_least_one
    mean_edep_given_positive = (sum_edep_positive / integral_positive) if integral_positive > 0 else 0.0
    theoretical_P_SE_given_edep = (sum_p_at_least_one_SE / integral_positive) if integral_positive > 0 else 0.0
    # P(≥1 SE | Edep > 0): use theoretical (histogram-based) as the physics expectation
    fraction_SE_given_edep = theoretical_P_SE_given_edep if integral_positive > 0 else 0.0
    # Expected fraction with SEE: sum over ALL bins so it is comparable to actual MC fraction.
    # E[P(≥1 SE)] = Σ_i (n_i / N_total) * P(≥1 SE | E_i). Includes the zero bin (small E can still yield SE).
    expected_fraction_with_SEE = 0.0
    for bin_i in range(1, edep_hist.GetNbinsX() + 1):
        center = edep_hist.GetBinCenter(bin_i)
        content = edep_hist.GetBinContent(bin_i)
        if content <= 0:
            continue
        mu_bin = (center / EPSILON) * P_ESC
        p_at_least_one = (1.0 - exp(-mu_bin)) if center > 0 else 0.0
        expected_fraction_with_SEE += (content / integral_total) * p_at_least_one
    
    # Mean Edep (given Edep > 0) from MC samples (histogram sampling only)
    mean_edep_MC_given_positive = float(np.mean(sampled_edep_when_positive)) if sampled_edep_when_positive else 0.0
    
    # For summary: use n_MC_edep_positive when histogram sampling so reported fraction is consistent
    n_events_with_edep_MC = n_MC_edep_positive if (use_histogram_sampling and n_MC_edep_positive > 0) else None
    
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
        'fraction_with_SE': fraction_with_SE,
        'fraction_with_edep': fraction_with_edep,
        'fraction_SE_given_edep': fraction_SE_given_edep,
        'mean_edep_given_positive': mean_edep_given_positive,
        'theoretical_P_SE_given_edep': theoretical_P_SE_given_edep,
        'expected_fraction_with_SEE': expected_fraction_with_SEE,
        'mean_edep_MC_given_positive': mean_edep_MC_given_positive,
        'n_MC_edep_positive': n_MC_edep_positive,
        'n_MC_edep_positive_with_SE': n_MC_edep_positive_with_SE,
        'sampled_edep_debug': sampled_edep_debug if use_histogram_sampling else None,
    }
    # Consistency: ⟨N_SE⟩ = ⟨μ⟩ = ⟨ΔE⟩/ε * P_esc. From sampled Edep, expected mean must match mean SEY.
    if use_histogram_sampling and sampled_edep_debug:
        mean_sampled_edep = float(np.mean(sampled_edep_debug))
        expected_mean_from_sampled = (mean_sampled_edep / EPSILON) * P_ESC
        statistics['mean_sampled_edep'] = mean_sampled_edep
        statistics['expected_mean_from_sampled'] = expected_mean_from_sampled
        # Sanity: Mean SEY (MC) must equal Expected (from sampled Edep) up to statistics
        rel_diff = abs(mean_SE - expected_mean_from_sampled) / (expected_mean_from_sampled + 1e-12)
        if rel_diff > 0.1:
            print("WARNING: Mean SEY (MC) = {:.4f} vs Expected (from sampled Edep) = {:.4f} (rel diff {:.1%}). Check P_esc applied in sampling.".format(
                mean_SE, expected_mean_from_sampled, rel_diff))
    
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
    
    # Build histogram of number of ionized electrons per event: N_int = Edep/ε
    # Fill per event from MC-sampled Edep so Entries = n_events (one entry per event).
    edep_xmax = edep_hist.GetXaxis().GetXmax()
    n_int_max = max(1.0, edep_xmax / EPSILON * 1.05)
    n_int_bins = min(400, max(100, int(n_int_max)))
    n_int_hist = TH1D("N_int", "Number of ionized electrons per event (N_{int} = #DeltaE/#varepsilon)",
                      n_int_bins, 0.0, n_int_max)
    n_int_hist.SetXTitle("N_{int} = #DeltaE / #varepsilon")
    n_int_hist.SetYTitle("Number of Events")
    n_int_hist.SetLineColor(ROOT.kOrange + 2)
    n_int_hist.SetFillColor(ROOT.kOrange + 2)
    n_int_hist.SetFillStyle(3004)
    if statistics.get("sampled_edep_debug"):
        for edep in statistics["sampled_edep_debug"]:
            if edep > 0:
                n_int_hist.Fill(edep / EPSILON)
    else:
        # Bin-by-bin mode: fill once per event using bin center and content
        for bin_i in range(1, edep_hist.GetNbinsX() + 1):
            center = edep_hist.GetBinCenter(bin_i)
            content = edep_hist.GetBinContent(bin_i)
            if content <= 0 or center <= 0:
                continue
            n_int_val = center / EPSILON
            for _ in range(int(content)):
                n_int_hist.Fill(n_int_val)
    
    # Create plot FIRST (before writing to file, to keep histogram valid)
    try:
        create_plot(sey_hist, statistics, root_file_path, n_int_hist=n_int_hist)
    except Exception as e:
        print("Warning: Could not create plot: {}".format(e))
        import traceback
        traceback.print_exc()
    
    # Create output file
    output_file_path = root_file_path.replace(".root", "_SEY_MonteCarlo.root")
    output_file = TFile.Open(output_file_path, "RECREATE")
    
    # Write histograms
    sey_hist.Write()
    n_int_hist.Write()
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
  Expected mean (theoretical, from histogram GetMean): {statistics['expected_mean']:.4f}
  Expected mean from sampled Edep (must match Mean SEY): {statistics.get('expected_mean_from_sampled', statistics['expected_mean']):.4f}
  Mean sampled Edep (MC): {statistics.get('mean_sampled_edep', statistics['mean_edep']):.2f} eV
  Standard deviation: {statistics['std_SE']:.4f}
  Median SEY: {statistics['median_SE']:.1f}
  Min SEY: {statistics['min_SE']}
  Max SEY: {statistics['max_SE']}
  
  Mean energy deposition: {statistics['mean_edep']:.2f} eV
  Mean free electrons per event: {statistics['mean_edep']/EPSILON:.2f}
  
  Events with at least one SE: {statistics['n_events_with_SE']}/{statistics['n_events']} ({statistics['fraction_with_SE']*100:.2f}%)
  
  Decomposition (histogram-based vs actual MC):
  - Fraction of events with Edep > 0 (histogram): {statistics['fraction_with_edep']*100:.2f}%
  - Mean Edep (given Edep > 0) from histogram: {statistics['mean_edep_given_positive']:.1f} eV
  - P(≥1 SE | Edep > 0) [theoretical from histogram]: {statistics['fraction_SE_given_edep']*100:.1f}%
  - Expected fraction with SEE (all bins): Σ (n_i/N)*P(≥1 SE|E_i) = {statistics['expected_fraction_with_SEE']*100:.2f}%
  - Actual fraction with SEE (MC): {statistics['fraction_with_SE']*100:.2f}%
  - Expected and actual are directly comparable (same population, same definition). Small differences are MC sampling noise. See doc/MUON_SEY_MONTE_CARLO.md.
  - For one ionization (Edep~ε=27 eV): P(≥1 SE)≈1-exp(-P_esc)≈32%.
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


def create_edep_debug_plot(sampled_edep_list, output_path):
    """
    Create a debug plot of MC-sampled energy deposition values with fine binning near 0.
    
    Parameters:
    -----------
    sampled_edep_list : list of float
        All Edep values sampled by the MC (GetRandom() from EdepPrimary).
    output_path : str
        Path for the output PDF (e.g. plots/.../file_Edep_debug.pdf).
    """
    if not sampled_edep_list:
        return
    # Fine binning: 0.5 eV from 0 to 250 eV (500 bins) to see behavior near 0.
    # Underflow captures values <= 0; overflow captures > 250 eV.
    nbins = 500
    xmin = 0.0
    xmax = 250.0  # eV
    edep_debug_hist = TH1D("Edep_MC_debug", "MC-sampled energy deposition (from EdepPrimary);Edep [eV];Events",
                            nbins, xmin, xmax)
    edep_debug_hist.SetDirectory(0)
    edep_debug_hist.SetFillColor(ROOT.kAzure + 2)
    edep_debug_hist.SetLineColor(ROOT.kAzure + 2)
    edep_debug_hist.SetLineWidth(1)
    for edep in sampled_edep_list:
        edep_debug_hist.Fill(edep)
    # Create canvas and draw
    c_debug = TCanvas("EdepDebug", "Edep debug", 900, 600)
    c_debug.SetGrid(1, 1)
    gStyle.SetOptStat(1110)  # entries, mean, RMS
    edep_debug_hist.Draw("HIST")
    # Add text: underflow/overflow if any
    n_under = edep_debug_hist.GetBinContent(0)
    n_over = edep_debug_hist.GetBinContent(nbins + 1)
    n_tot = len(sampled_edep_list)
    text = TLatex()
    text.SetNDC(True)
    text.SetTextSize(0.03)
    text.DrawLatex(0.15, 0.85, "Fine binning: 0.5 eV per bin (0--250 eV)")
    if n_under > 0 or n_over > 0:
        text.DrawLatex(0.15, 0.81, "Underflow (Edep #leq 0): %d   Overflow (>250 eV): %d" % (n_under, n_over))
    text.DrawLatex(0.15, 0.77, "Total sampled: %d" % n_tot)
    c_debug.SaveAs(output_path)
    root_path = output_path.replace(".pdf", ".root")
    c_debug.SaveAs(root_path)
    print("Edep debug plot saved to: {} and {}".format(output_path, root_path))


def create_plot(sey_hist, statistics, input_file_path, n_int_hist=None):
    """
    Create and save plot of SEY distribution (and optionally N_int distribution).
    
    Parameters:
    -----------
    sey_hist : TH1D
        Histogram of SEY values
    statistics : dict
        Statistics dictionary
    input_file_path : str
        Input file path (for output naming)
    n_int_hist : TH1D or None
        Optional histogram of number of ionized electrons per event (N_int = Edep/ε)
    """
    gStyle.SetOptStat(1110)  # Show mean, RMS, entries (top-right)
    gStyle.SetOptFit(0)
    
    # Remove any fit (e.g. Poisson) attached to the histogram so it is not drawn
    flist = sey_hist.GetListOfFunctions()
    if flist:
        flist.Clear()
    
    canvas = TCanvas("SEYCanvas", "Secondary Electron Yield (Monte Carlo)", 1000, 700)
    canvas.SetGrid()
    canvas.cd()
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
    
    # Debug plot: MC-sampled Edep values (fine binning near 0)
    if statistics.get("sampled_edep_debug"):
        debug_basename = input_basename.replace(".root", "_Edep_debug")
        debug_plot_path = os.path.join(plots_subdir, debug_basename + ".pdf")
        create_edep_debug_plot(statistics["sampled_edep_debug"], debug_plot_path)
    
    # Save plot (with text added below); multi-page if n_int_hist provided
    print("\nPlot saved to: {}".format(plot_file_path))
    
    # Try to add text, but don't fail if it doesn't work
    try:
        # Results: box with titled sections, right of histogram and large enough for text
        rbox = TPaveText(0.21, 0.40, 0.56, 0.88, "NDC")
        rbox.SetFillColor(0)
        rbox.SetBorderSize(1)
        rbox.SetTextAlign(12)
        rbox.SetTextSize(0.022)
        rbox.AddText("Monte Carlo SEY Calculation")
        rbox.AddText(" ")
        rbox.AddText("Results")
        rbox.AddText("  Mean SEY: %.4f" % statistics['mean_SE'])
        rbox.AddText("  Total SE: %d" % statistics['total_SE'])
        rbox.AddText("  Events with SE (actual MC): %.2f%% (%d/%d)" % (
            statistics['fraction_with_SE'] * 100.0,
            statistics['n_events_with_SE'],
            statistics['n_events']))
        rbox.AddText(" ")
        rbox.AddText("Check: Mean SEY = Expected")
        rbox.AddText("  Expected (histogram): %.4f" % statistics['expected_mean'])
        if statistics.get('expected_mean_from_sampled') is not None:
            rbox.AddText("  Expected (from sampled Edep): %.4f" % statistics['expected_mean_from_sampled'])
        rbox.AddText(" ")
        if 'fraction_with_edep' in statistics and 'fraction_SE_given_edep' in statistics and 'expected_fraction_with_SEE' in statistics:
            rbox.AddText("Check: fraction with SEE")
            rbox.AddText("  Edep>0: %.2f%%,  P(#geq 1 SE|Edep>0): %.1f%%" % (
                statistics['fraction_with_edep'] * 100.0,
                statistics['fraction_SE_given_edep'] * 100.0))
            rbox.AddText("  Expected (all bins #Sigma): %.2f%%" % (statistics['expected_fraction_with_SEE'] * 100.0))
            rbox.AddText("  Actual (MC): %.2f%%" % (statistics['fraction_with_SE'] * 100.0))
        rbox.Draw()
        
        # Parameters: right box, just below the stat box (top-right)
        pbox = TPaveText(0.58, 0.52, 0.88, 0.72, "NDC")
        pbox.SetFillColor(0)
        pbox.SetBorderSize(1)
        pbox.SetTextAlign(12)
        pbox.SetTextSize(0.022)
        pbox.AddText("Parameters:")
        pbox.AddText("P_{esc}(z) = %.4f" % P_ESC)
        pbox.AddText("Production depth: z = %.1f nm" % (Z_DEPTH/10))
        pbox.AddText("Events: %d" % statistics['n_events'])
        pbox.Draw()
        
        canvas.Update()
        # PDF: single page (SaveAs) or multi-page (Print): page 1 = SEY, page 2 = N_int
        if n_int_hist and n_int_hist.GetEntries() > 0:
            canvas.Print(plot_file_path + "(")
            c2 = TCanvas("cNint", "Number of ionized electrons per event", 1000, 700)
            c2.SetGrid()
            c2.SetLogy(1)
            c2.SetLeftMargin(0.12)
            n_int_hist.GetYaxis().SetTitleOffset(1.2)
            n_int_hist.Draw("HIST")
            rbox2 = TPaveText(0.15, 0.75, 0.48, 0.88, "NDC")
            rbox2.SetFillColor(0)
            rbox2.SetBorderSize(1)
            rbox2.SetTextAlign(12)
            rbox2.SetTextSize(0.022)
            rbox2.AddText("N_{int} = #DeltaE / #varepsilon  (#varepsilon = %.0f eV)" % EPSILON)
            rbox2.Draw()
            pbox2 = TPaveText(0.58, 0.52, 0.88, 0.72, "NDC")
            pbox2.SetFillColor(0)
            pbox2.SetBorderSize(1)
            pbox2.SetTextAlign(12)
            pbox2.SetTextSize(0.022)
            pbox2.AddText("Parameters:")
            pbox2.AddText("P_{esc}(z) = %.4f" % P_ESC)
            pbox2.AddText("Production depth: z = %.1f nm" % (Z_DEPTH/10))
            pbox2.AddText("Events: %d" % statistics['n_events'])
            pbox2.Draw()
            c2.Update()
            c2.Print(plot_file_path + ")")
        else:
            canvas.SaveAs(plot_file_path)
    except Exception as e:
        print("Note: Text annotations skipped (histogram saved successfully)")
        # PDF may be partial; try to close if multi-page was started
        try:
            canvas.Print(plot_file_path + ")")
        except Exception:
            pass
    
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
