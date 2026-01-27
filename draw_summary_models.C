// ROOT macro to compare EM models per energy value.
// Usage:
//   root -l -b -q 'draw_summary_models.C("dirPAI","dirLiv","dirPen","summary_models.root")'

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLatex.h"

#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cmath>

namespace {

std::vector<std::string> ListRootFiles(const char* dirPath) {
    std::vector<std::string> files;
    TSystemDirectory dir("scan_dir", dirPath);
    TList* list = dir.GetListOfFiles();
    if (!list) {
        return files;
    }
    TIter next(list);
    while (auto* file = dynamic_cast<TSystemFile*>(next())) {
        TString name = file->GetName();
        if (file->IsDirectory()) {
            continue;
        }
        if (name.EndsWith(".root") && !name.BeginsWith("summary")) {
            files.push_back(std::string(dirPath) + "/" + name.Data());
        }
    }
    files.shrink_to_fit();
    return files;
}

std::string ModelLabel(const std::string& name) {
    if (name.empty()) return "unknown";
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (lower == "pai") return "PAI";
    if (lower == "livermore" || lower == "livermorephysics" ||
        lower == "g4emlivermorephysics") return "Livermore";
    if (lower == "penelope" || lower == "penelopephysics" ||
        lower == "g4empenelopephysics") return "Penelope";
    return name;
}

std::string ParticleLabel(const std::string& name) {
    if (name == "e-" || name == "e") return "e^{-}";
    if (name == "e+") return "e^{+}";
    if (name == "mu-" || name == "mu") return "#mu^{-}";
    if (name == "mu+") return "#mu^{+}";
    if (!name.empty()) return name;
    return "e^{-}";
}

TString EnergyLabel(double energyMeV) {
    if (energyMeV <= 0.) {
        return "E=n/a";
    }
    if (energyMeV < 1e-3) {
        return Form("E=%.0f eV", energyMeV * 1.0e6);
    }
    if (energyMeV < 1.0) {
        return Form("E=%.3g keV", energyMeV * 1.0e3);
    }
    return Form("E=%.3g MeV", energyMeV);
}

int ModelOrder(const std::string& model) {
    std::string label = ModelLabel(model);
    if (label == "PAI") return 0;
    if (label == "Livermore") return 1;
    if (label == "Penelope") return 2;
    return 3;
}

void RebinToMaxBins(TH1* hist, int maxBins) {
    if (!hist || maxBins <= 0) {
        return;
    }
    const int nbins = hist->GetNbinsX();
    if (nbins <= maxBins) {
        return;
    }
    const int factor = static_cast<int>(std::ceil(nbins / static_cast<double>(maxBins)));
    if (factor > 1) {
        hist->Rebin(factor);
    }
}

} // namespace

// Dynamic rebinning and range setting for EdepPrimary
// Creates a new histogram with appropriate binning for the actual data range
TH1* OptimizeEdepPrimaryAxis(TH1* hist) {
    if (!hist || hist->GetEntries() == 0) return hist;
    
    // Find actual data range
    // Use the histogram's mean and RMS to determine the data range
    const double mean = hist->GetMean();
    const double rms = hist->GetRMS();
    const double entries = hist->GetEntries();
    
    if (entries == 0 || mean <= 0) return hist;
    
    // If RMS is very small (all data in one bin or very narrow), use mean-based range
    double dataMin, dataMax, range;
    if (rms < mean * 0.01) {
        // All data clustered around mean - use mean ± reasonable range
        // For very narrow distributions, ensure we capture the actual data range
        // Find the actual min and max from non-empty bins
        double actualMin = hist->GetXaxis()->GetXmax();
        double actualMax = hist->GetXaxis()->GetXmin();
        bool foundData = false;
        for (int i = 1; i <= hist->GetNbinsX(); i++) {
            if (hist->GetBinContent(i) > 0) {
                foundData = true;
                double binLow = hist->GetXaxis()->GetBinLowEdge(i);
                double binUp = hist->GetXaxis()->GetBinUpEdge(i);
                if (binLow < actualMin) actualMin = binLow;
                if (binUp > actualMax) actualMax = binUp;
            }
        }
        if (foundData) {
            // Use actual data range with padding
            range = std::max((actualMax - actualMin) * 1.2, std::max(mean * 0.2, 1.0));
            dataMin = std::max(0.0, actualMin - range * 0.1);
            dataMax = actualMax + range * 0.1;
        } else {
            // Fallback to mean-based range
            range = std::max(mean * 0.2, 1.0);
            dataMin = std::max(0.0, mean - range);
            dataMax = mean + range;
        }
    } else {
        // Data has spread - use mean ± 3*RMS for range
        dataMin = std::max(0.0, mean - 3.0 * rms);
        dataMax = mean + 3.0 * rms;
        range = dataMax - dataMin;
    }
    const double currentBinWidth = hist->GetBinWidth(1);
    const double currentRange = hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin();
    
    // If the histogram range is much larger than the data range (e.g., muons with 0-4e9 eV range but ~3 eV data)
    // or if binning is too coarse, create a new histogram with better binning
    const bool needsRebinning = (currentRange > range * 10.0) || (currentBinWidth > range / 50.0);
    
    if (needsRebinning) {
        // Determine appropriate binning: aim for ~200 bins in the data range
        const int targetBins = 200;
        const int minBins = 50; // Minimum bins to ensure good visualization
        double idealBinWidth = range / targetBins;
        
        // Round to a nice bin width (powers of 10 or common fractions)
        double binWidth = idealBinWidth;
        if (binWidth < 0.01) {
            binWidth = 0.01; // 0.01 eV bins for very small ranges
        } else if (binWidth < 0.1) {
            binWidth = 0.1; // 0.1 eV bins
        } else if (binWidth < 1.0) {
            binWidth = 1.0; // 1 eV bins
        } else if (binWidth < 10.0) {
            binWidth = 10.0; // 10 eV bins
        } else {
            // Round to nearest power of 10
            double log10 = std::log10(binWidth);
            binWidth = std::pow(10.0, std::floor(log10 + 0.5));
        }
        
        // Add padding (20% on each side, but at least 2*binWidth)
        const double padding = std::max(range * 0.2, binWidth * 2.0);
        const double newMin = std::max(0.0, dataMin - padding);
        const double newMax = dataMax + padding;
        const double newRange = newMax - newMin;
        
        // Calculate number of bins, ensuring minimum
        int newNbins = static_cast<int>(std::ceil(newRange / binWidth));
        if (newNbins < minBins) {
            // Adjust bin width to get at least minBins
            binWidth = newRange / minBins;
            newNbins = minBins;
        }
        
        // Ensure we have at least 1 bin
        if (newNbins < 1) newNbins = 1;
        
        // Create new histogram with appropriate range and binning
        TH1* newHist = new TH1D(Form("%s_optimized", hist->GetName()),
                               hist->GetTitle(),
                               newNbins, newMin, newMax);
        newHist->SetDirectory(nullptr);
        
        // Copy style
        newHist->SetLineColor(hist->GetLineColor());
        newHist->SetLineWidth(hist->GetLineWidth());
        newHist->SetFillColor(hist->GetFillColor());
        newHist->SetFillStyle(hist->GetFillStyle());
        
        // Copy axis titles
        newHist->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
        newHist->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
        
        // Fill new histogram from original histogram's actual data
        // Properly map bins from original to new histogram to avoid gaps
        // For very narrow distributions (all data in one bin), use mean-based filling to preserve data location
        const double mean = hist->GetMean();
        const double rms = hist->GetRMS();
        const bool isVeryNarrow = (rms < mean * 0.01) && (hist->GetNbinsX() > 1000); // Very narrow and original histogram has many bins
        
        if (isVeryNarrow) {
            // For very narrow distributions, fill at the mean location to preserve the actual data position
            // This avoids artificial spreading when a single wide bin maps to multiple narrow bins
            if (mean >= newMin && mean <= newMax) {
                int bin = newHist->GetXaxis()->FindBin(mean);
                if (bin >= 1 && bin <= newHist->GetNbinsX()) {
                    newHist->SetBinContent(bin, hist->GetEntries());
                    newHist->SetEntries(hist->GetEntries());
                }
            }
        } else {
            // For wider distributions, use proper bin mapping
            for (int i = 1; i <= hist->GetNbinsX(); i++) {
                const double content = hist->GetBinContent(i);
                if (content > 0) {
                    const double binLow = hist->GetXaxis()->GetBinLowEdge(i);
                    const double binUp = hist->GetXaxis()->GetBinUpEdge(i);
                    const double binCenter = hist->GetXaxis()->GetBinCenter(i);
                    
                    // Only process if bin overlaps with new histogram range
                    if (binUp >= newMin && binLow <= newMax) {
                        // Find which bins in the new histogram this bin overlaps with
                        int newBinLow = newHist->GetXaxis()->FindBin(binLow);
                        int newBinUp = newHist->GetXaxis()->FindBin(binUp);
                        
                        // Ensure bins are within valid range
                        if (newBinLow < 1) newBinLow = 1;
                        if (newBinUp > newHist->GetNbinsX()) newBinUp = newHist->GetNbinsX();
                        
                        if (newBinLow == newBinUp) {
                            // Single bin - fill directly
                            newHist->SetBinContent(newBinLow, newHist->GetBinContent(newBinLow) + content);
                        } else {
                            // Multiple bins - distribute content proportionally
                            // Calculate overlap fraction for each bin
                            for (int j = newBinLow; j <= newBinUp; j++) {
                                double newBinLowEdge = newHist->GetXaxis()->GetBinLowEdge(j);
                                double newBinUpEdge = newHist->GetXaxis()->GetBinUpEdge(j);
                                
                                // Calculate overlap
                                double overlapLow = std::max(binLow, newBinLowEdge);
                                double overlapUp = std::min(binUp, newBinUpEdge);
                                double overlap = std::max(0.0, overlapUp - overlapLow);
                                double originalBinWidth = binUp - binLow;
                                
                                if (overlap > 0 && originalBinWidth > 0) {
                                    double fraction = overlap / originalBinWidth;
                                    newHist->SetBinContent(j, newHist->GetBinContent(j) + content * fraction);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Update entries to match sum of bin contents
        double totalContent = 0;
        for (int i = 1; i <= newHist->GetNbinsX(); i++) {
            totalContent += newHist->GetBinContent(i);
        }
        newHist->SetEntries(totalContent);
        
        // If no data was filled (shouldn't happen, but safety check),
        // fall back to mean-based filling
        if (totalContent == 0) {
            const double histMean = hist->GetMean();
            const double totalEntries = hist->GetEntries();
            if (totalEntries > 0 && histMean > 0 && histMean >= newMin && histMean <= newMax) {
                int bin = newHist->GetXaxis()->FindBin(histMean);
                if (bin >= 1 && bin <= newHist->GetNbinsX()) {
                    newHist->SetBinContent(bin, totalEntries);
                    newHist->SetEntries(totalEntries);
                }
            }
        }
        
        return newHist;
    } else {
        // Just set the range appropriately
        const double padding = range * 0.1;
        hist->GetXaxis()->SetRangeUser(std::max(0.0, dataMin - padding), dataMax + padding);
        return hist;
    }
}

void draw_summary_models(const char* dirPai,
                         const char* dirLivermore,
                         const char* dirPenelope,
                         const char* outputFile = "summary_models.root") {
    if (!dirPai || !dirLivermore || !dirPenelope) {
        printf("Error: must provide three model directories.\n");
        return;
    }

    std::vector<std::string> dirs = {dirPai, dirLivermore, dirPenelope};
    std::vector<std::string> allFiles;
    for (const auto& dir : dirs) {
        auto files = ListRootFiles(dir.c_str());
        allFiles.insert(allFiles.end(), files.begin(), files.end());
    }
    if (allFiles.empty()) {
        printf("No .root files found in provided directories.\n");
        return;
    }

    struct Entry {
        std::string path;
        std::string model;
        std::string particle;
        double energyMeV;
        double thicknessNm;
    };

    std::vector<Entry> entries;
    entries.reserve(allFiles.size());

    for (const auto& path : allFiles) {
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) {
            continue;
        }
        double energyMeV = -1.0;
        double thicknessNm = -1.0;
        char particle[64] = "";
        char emModel[64] = "";
        TTree* meta = (TTree*)f->Get("RunMeta");
        if (meta && meta->GetEntries() > 0) {
            if (meta->GetBranch("primaryEnergyMeV")) {
                meta->SetBranchAddress("primaryEnergyMeV", &energyMeV);
            }
            if (meta->GetBranch("sampleThicknessNm")) {
                meta->SetBranchAddress("sampleThicknessNm", &thicknessNm);
            }
            if (meta->GetBranch("primaryParticle")) {
                meta->SetBranchAddress("primaryParticle", particle);
            }
            if (meta->GetBranch("emModel")) {
                meta->SetBranchAddress("emModel", emModel);
            }
            meta->GetEntry(0);
        }
        f->Close();
        if (energyMeV <= 0.) {
            continue;
        }
        entries.push_back({path, ModelLabel(emModel), particle, energyMeV, thicknessNm});
    }

    if (entries.empty()) {
        printf("No entries with valid RunMeta found.\n");
        return;
    }

    std::map<long long, std::vector<Entry>> byEnergy;
    for (const auto& entry : entries) {
        const long long energyEv = std::llround(entry.energyMeV * 1.0e6);
        byEnergy[energyEv].push_back(entry);
    }

    gStyle->SetOptStat(0);

    std::map<std::string, int> modelColor = {
        {"PAI", kMagenta + 1},
        {"Livermore", kAzure + 2},
        {"Penelope", kOrange + 7}
    };

    TFile* out = TFile::Open(outputFile, "RECREATE");
    if (!out || out->IsZombie()) {
        printf("Error: cannot create output file %s\n", outputFile);
        return;
    }

    std::string plotsDir;
    if (outputFile && outputFile[0] != '\0') {
        std::string outPath = outputFile;
        std::string baseDir;
        const std::string marker = "/results/";
        const auto pos = outPath.rfind(marker);
        if (pos != std::string::npos) {
            baseDir = outPath.substr(0, pos);
        } else {
            baseDir = gSystem->DirName(outPath.c_str());
            if (baseDir == "results") {
                baseDir = ".";
            }
        }
        std::string fileBase = gSystem->BaseName(outPath.c_str());
        if (fileBase.size() >= 5 && fileBase.substr(fileBase.size() - 5) == ".root") {
            fileBase = fileBase.substr(0, fileBase.size() - 5);
        }
        if (baseDir.empty() || baseDir == ".") {
            plotsDir = "plots/" + fileBase;
        } else {
            plotsDir = baseDir + "/plots/" + fileBase;
        }
        gSystem->mkdir(plotsDir.c_str(), true);
    }
    auto saveCanvas = [&](TCanvas* canvas, const TString& name) {
        if (!canvas || plotsDir.empty() || name.IsNull()) {
            return;
        }
        std::string rootPath = plotsDir + "/" + name.Data() + ".root";
        std::string pdfPath = plotsDir + "/" + name.Data() + ".pdf";
        canvas->SaveAs(rootPath.c_str());
        canvas->SaveAs(pdfPath.c_str());
    };

    for (const auto& kv : byEnergy) {
        const auto& group = kv.second;
        if (group.empty()) {
            continue;
        }

        std::vector<Entry> sorted = group;
        std::sort(sorted.begin(), sorted.end(),
                  [](const Entry& a, const Entry& b) {
                      int ao = ModelOrder(a.model);
                      int bo = ModelOrder(b.model);
                      if (ao != bo) return ao < bo;
                      return a.model < b.model;
                  });

        const std::string particleText = ParticleLabel(sorted.front().particle);
        const double thicknessNm = sorted.front().thicknessNm;
        const double energyMeV = sorted.front().energyMeV;

        TString titleSuffix = Form("%s, t=%.2f nm, %s",
                                   particleText.c_str(),
                                   thicknessNm > 0. ? thicknessNm : 0.0,
                                   EnergyLabel(energyMeV).Data());

        TString canvasName = Form("EdepPrimaryModels_E%lld", kv.first);
        TCanvas* c1 = new TCanvas(canvasName, "EdepPrimary models", 900, 700);
        c1->SetGrid();
        c1->SetLogy();
        c1->SetLeftMargin(0.15);

        TLegend* leg1 = new TLegend(0.45, 0.65, 0.78, 0.88);
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        leg1->SetTextSize(0.03);

        // Find max range from optimized histograms
        // First pass: optimize all histograms and find the actual last non-empty bin
        double maxPrimaryEdgeEv = 0.0;
        std::vector<TH1*> tempOptimized;
        for (const auto& entry : sorted) {
            TFile* f = TFile::Open(entry.path.c_str(), "READ");
            if (!f || f->IsZombie()) {
                continue;
            }
            TH1* h = (TH1*)f->Get("EdepPrimary");
            if (h) {
                // Clone first to detach from file
                TH1* hClone = (TH1*)h->Clone("temp");
                hClone->SetDirectory(nullptr);
                f->Close();
                
                // Now optimize the cloned histogram
                TH1* hOpt = OptimizeEdepPrimaryAxis(hClone);
                bool hOptCreated = (hOpt != nullptr && hOpt != hClone);
                
                if (hOpt) {
                    // Find the actual last non-empty bin in the optimized histogram
                    const int lastBin = hOpt->FindLastBinAbove(0.0);
                    if (lastBin > 0) {
                        const double edge = hOpt->GetXaxis()->GetBinUpEdge(lastBin);
                        if (edge > maxPrimaryEdgeEv) {
                            maxPrimaryEdgeEv = edge;
                        }
                    }
                    // Store for cleanup
                    if (hOptCreated) {
                        tempOptimized.push_back(hOpt);
                    }
                }
                
                // Clean up clone
                delete hClone;
            } else {
                f->Close();
            }
        }
        
        // Clean up temporary optimized histograms
        for (auto* h : tempOptimized) {
            delete h;
        }

        bool first = true;
        for (const auto& entry : sorted) {
            TFile* f = TFile::Open(entry.path.c_str(), "READ");
            if (!f || f->IsZombie()) {
                continue;
            }
            TH1* h = (TH1*)f->Get("EdepPrimary");
            if (!h) {
                f->Close();
                continue;
            }
            
            // Clone first to detach from file
            TH1* hClone = (TH1*)h->Clone("temp");
            hClone->SetDirectory(nullptr);
            f->Close();
            
            // Optimize the histogram
            TH1* hOpt = OptimizeEdepPrimaryAxis(hClone);
            bool hOptCreated = (hOpt != nullptr && hOpt != hClone);
            
            TString cloneName = Form("EdepPrimary_%s_E%lld", entry.model.c_str(), kv.first);
            TH1* hc = (TH1*)hOpt->Clone(cloneName);
            if (!hc) {
                if (hOptCreated && hOpt) delete hOpt;
                delete hClone;
                continue;
            }
            hc->SetDirectory(nullptr);
            
            // Clean up temporary histograms
            if (hOptCreated && hOpt) delete hOpt;
            delete hClone;
            
            RebinToMaxBins(hc, 2000);
            int color = modelColor.count(entry.model) ? modelColor[entry.model] : kBlack;
            hc->SetLineColor(color);
            hc->SetLineWidth(4);
            
            if (first) {
                hc->SetTitle(Form("EdepPrimary model comparison (%s)", titleSuffix.Data()));
                hc->GetYaxis()->SetTitleOffset(1.35);
                hc->GetYaxis()->SetLabelSize(0.035);
                hc->Draw("HIST");
                first = false;
            } else {
                hc->Draw("HIST SAME");
            }
            leg1->AddEntry(hc, entry.model.c_str(), "l");
        }
        
        // Now find the actual last non-empty bin across all drawn histograms and trim
        double actualMaxEdge = 0.0;
        TList* prims = c1->GetListOfPrimitives();
        if (prims) {
            TObject* obj = prims->First();
            while (obj) {
                if (obj->InheritsFrom("TH1")) {
                    TH1* h = (TH1*)obj;
                    if (TString(h->GetName()).Contains("EdepPrimary")) {
                        const int lastBin = h->FindLastBinAbove(0.0);
                        if (lastBin > 0) {
                            const double edge = h->GetXaxis()->GetBinUpEdge(lastBin);
                            if (edge > actualMaxEdge) {
                                actualMaxEdge = edge;
                            }
                        }
                    }
                }
                obj = prims->After(obj);
            }
            
            // Set range to trim after the last non-empty bin (with 5% padding)
            if (actualMaxEdge > 0.0) {
                const double paddedMax = actualMaxEdge * 1.05;
                obj = prims->First();
                while (obj) {
                    if (obj->InheritsFrom("TH1")) {
                        TH1* h = (TH1*)obj;
                        if (TString(h->GetName()).Contains("EdepPrimary")) {
                            h->GetXaxis()->SetRangeUser(0.0, paddedMax);
                        }
                    }
                    obj = prims->After(obj);
                }
                // Also set the range on the pad to ensure it's applied
                TPad* pad = (TPad*)c1->GetPad(1);
                if (pad) {
                    pad->Modified();
                    pad->Update();
                }
            }
            
            // Find maximum bin content across all histograms to ensure y-axis accommodates saturated bins
            double maxBinContent = 0.0;
            obj = prims->First();
            while (obj) {
                if (obj->InheritsFrom("TH1")) {
                    TH1* h = (TH1*)obj;
                    if (TString(h->GetName()).Contains("EdepPrimary")) {
                        // Check all bins in the visible range (including overflow if visible)
                        const int firstBin = h->GetXaxis()->GetFirst();
                        const int lastBin = h->GetXaxis()->GetLast();
                        for (int i = firstBin; i <= lastBin; i++) {
                            const double content = h->GetBinContent(i);
                            if (content > maxBinContent) {
                                maxBinContent = content;
                            }
                        }
                        // Also check overflow bin if it's within the visible range
                        const int overflowBin = h->GetNbinsX() + 1;
                        if (overflowBin <= lastBin) {
                            const double overflowContent = h->GetBinContent(overflowBin);
                            if (overflowContent > maxBinContent) {
                                maxBinContent = overflowContent;
                            }
                        }
                    }
                }
                obj = prims->After(obj);
            }
            
            // Set y-axis range to accommodate maximum bin content with padding
            if (maxBinContent > 0.0) {
                // Add 20% padding above the maximum to ensure saturated bins are visible
                const double yMax = maxBinContent * 1.2;
                // For log scale, ensure minimum is appropriate
                const double yMin = c1->GetLogy() ? std::max(0.1, maxBinContent * 1e-6) : 0.0;
                
                obj = prims->First();
                while (obj) {
                    if (obj->InheritsFrom("TH1")) {
                        TH1* h = (TH1*)obj;
                        if (TString(h->GetName()).Contains("EdepPrimary")) {
                            h->SetMinimum(yMin);
                            h->SetMaximum(yMax);
                        }
                    }
                    obj = prims->After(obj);
                }
            }
        }
        
        leg1->Draw();
        c1->Update();

        // Create a zoomed version of EdepPrimaryModels (0-100eV range) for 1MeV electrons and 4GeV muons
        TCanvas* c1zoom = nullptr;
        TString canvasNameEdepZoom;
        const long long energyEv = kv.first;
        // Check if this is 1MeV (1000000 eV) for electrons or 4GeV (4000000000 eV) for muons
        if (energyEv == 1000000 || energyEv == 4000000000) {
            canvasNameEdepZoom = Form("EdepPrimaryModelsZoom_E%lld", kv.first);
            c1zoom = new TCanvas(canvasNameEdepZoom, "EdepPrimary models (zoom 0-100eV)", 900, 700);
            c1zoom->SetGrid();
            c1zoom->SetLogy();
            c1zoom->SetLeftMargin(0.15);

            TLegend* leg1zoom = new TLegend(0.45, 0.65, 0.78, 0.88);
            leg1zoom->SetBorderSize(0);
            leg1zoom->SetFillStyle(0);
            leg1zoom->SetTextSize(0.03);

            // Reuse the same histograms but with zoomed range
            bool firstEdepZoom = true;
            for (const auto& entry : sorted) {
                TFile* f = TFile::Open(entry.path.c_str(), "READ");
                if (!f || f->IsZombie()) {
                    continue;
                }
                TH1* h = (TH1*)f->Get("EdepPrimary");
                if (!h) {
                    f->Close();
                    continue;
                }
                
                // Clone first to detach from file
                TH1* hClone = (TH1*)h->Clone("temp");
                hClone->SetDirectory(nullptr);
                f->Close();
                
                // Optimize the histogram
                TH1* hOpt = OptimizeEdepPrimaryAxis(hClone);
                bool hOptCreated = (hOpt != nullptr && hOpt != hClone);
                
                TString cloneName = Form("EdepPrimaryZoom_%s_E%lld", entry.model.c_str(), kv.first);
                TH1* hc = (TH1*)hOpt->Clone(cloneName);
                if (!hc) {
                    if (hOptCreated && hOpt) delete hOpt;
                    delete hClone;
                    continue;
                }
                hc->SetDirectory(nullptr);
                
                // Clean up temporary histograms
                if (hOptCreated && hOpt) delete hOpt;
                delete hClone;
                
                RebinToMaxBins(hc, 2000);
                int color = modelColor.count(entry.model) ? modelColor[entry.model] : kBlack;
                hc->SetLineColor(color);
                hc->SetLineWidth(4);
                
                // Set zoomed range (0-100eV)
                hc->GetXaxis()->SetRangeUser(0.0, 100.0);
                
                if (firstEdepZoom) {
                    hc->SetTitle(Form("Primary energy deposition in Al_{2}O_{3} - Zoom 0-100eV (%s)", titleSuffix.Data()));
                    hc->GetYaxis()->SetTitleOffset(1.35);
                    hc->GetYaxis()->SetLabelSize(0.035);
                    hc->Draw("HIST");
                    firstEdepZoom = false;
                } else {
                    hc->Draw("HIST SAME");
                }
                leg1zoom->AddEntry(hc, entry.model.c_str(), "l");
            }
            
            leg1zoom->Draw();
            c1zoom->Update();
        }

        TString canvasNameSteps = Form("EdepInteractionsModels_E%lld", kv.first);
        TCanvas* c2 = new TCanvas(canvasNameSteps, "EdepInteractions models", 900, 700);
        c2->SetGrid();
        c2->SetLogy();
        c2->SetLeftMargin(0.15);

        TLegend* leg2 = new TLegend(0.65, 0.65, 0.88, 0.88);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetTextSize(0.03);

        double maxStepsEdge = 0.0;
        for (const auto& entry : sorted) {
            TFile* f = TFile::Open(entry.path.c_str(), "READ");
            if (!f || f->IsZombie()) {
                continue;
            }
            TH1* h = (TH1*)f->Get("EdepInteractions");
            if (h) {
                const int lastBin = h->FindLastBinAbove(0.0);
                if (lastBin > 0) {
                    const double edge = h->GetXaxis()->GetBinUpEdge(lastBin);
                    if (edge > maxStepsEdge) {
                        maxStepsEdge = edge;
                    }
                }
            }
            f->Close();
        }

        bool firstSteps = true;
        for (const auto& entry : sorted) {
            TFile* f = TFile::Open(entry.path.c_str(), "READ");
            if (!f || f->IsZombie()) {
                continue;
            }
            TH1* h = (TH1*)f->Get("EdepInteractions");
            if (!h) {
                f->Close();
                continue;
            }
            TString cloneName = Form("EdepInteractions_%s_E%lld", entry.model.c_str(), kv.first);
            TH1* hc = (TH1*)h->Clone(cloneName);
            if (!hc) {
                f->Close();
                continue;
            }
            hc->SetDirectory(nullptr);
            f->Close();
            RebinToMaxBins(hc, 2000);
            int color = modelColor.count(entry.model) ? modelColor[entry.model] : kBlack;
            hc->SetLineColor(color);
            hc->SetLineWidth(4);
            if (firstSteps) {
                hc->SetTitle(Form("EdepInteractions model comparison (%s)", titleSuffix.Data()));
                hc->GetYaxis()->SetTitleOffset(1.35);
                hc->GetYaxis()->SetLabelSize(0.035);
                if (maxStepsEdge > 0.0) {
                    hc->GetXaxis()->SetRangeUser(0.0, maxStepsEdge);
                }
                hc->Draw("HIST");
                firstSteps = false;
            } else {
                hc->Draw("HIST SAME");
            }
            leg2->AddEntry(hc, entry.model.c_str(), "l");
        }
        leg2->Draw();
        c2->Update();

        TString canvasNameStepLen = Form("StepLengthModels_E%lld", kv.first);
        TCanvas* c3 = new TCanvas(canvasNameStepLen, "Step length models", 900, 700);
        c3->SetGrid();
        c3->SetLogy();
        c3->SetLeftMargin(0.15);

        TLegend* leg3 = new TLegend(0.45, 0.65, 0.78, 0.88);
        leg3->SetBorderSize(0);
        leg3->SetFillStyle(0);
        leg3->SetTextSize(0.03);

        // First pass: find the actual last non-empty bin across all histograms
        // First pass: find the maximum last non-empty bin edge across all histograms
        double maxStepLenEdge = 0.0;
        for (const auto& entry : sorted) {
            TFile* f = TFile::Open(entry.path.c_str(), "READ");
            if (!f || f->IsZombie()) {
                continue;
            }
            TH1* h = (TH1*)f->Get("StepLengthAl2O3");
            if (h) {
                const int lastBin = h->FindLastBinAbove(0.0);
                if (lastBin > 0) {
                    const double edge = h->GetXaxis()->GetBinUpEdge(lastBin);
                    if (edge > maxStepLenEdge) {
                        maxStepLenEdge = edge;
                    }
                }
            }
            f->Close();
        }

        // Calculate the range to use (trim after last non-empty bin with 5% padding)
        double stepLenRangeMax = 0.0;
        if (maxStepLenEdge > 0.0) {
            stepLenRangeMax = maxStepLenEdge * 1.05;
        }

        bool firstStepLen = true;
        for (const auto& entry : sorted) {
            TFile* f = TFile::Open(entry.path.c_str(), "READ");
            if (!f || f->IsZombie()) {
                continue;
            }
            TH1* h = (TH1*)f->Get("StepLengthAl2O3");
            if (!h) {
                f->Close();
                continue;
            }
            TString cloneName = Form("StepLength_%s_E%lld", entry.model.c_str(), kv.first);
            TH1* hc = (TH1*)h->Clone(cloneName);
            if (!hc) {
                f->Close();
                continue;
            }
            hc->SetDirectory(nullptr);
            f->Close();
            RebinToMaxBins(hc, 2000);
            int color = modelColor.count(entry.model) ? modelColor[entry.model] : kBlack;
            hc->SetLineColor(color);
            hc->SetLineWidth(4);
            if (firstStepLen) {
                hc->SetTitle(Form("Step length in Al_{2}O_{3} (%s)", titleSuffix.Data()));
                hc->GetYaxis()->SetTitleOffset(1.35);
                hc->GetYaxis()->SetLabelSize(0.035);
                hc->Draw("HIST");
                firstStepLen = false;
            } else {
                hc->Draw("HIST SAME");
            }
            leg3->AddEntry(hc, entry.model.c_str(), "l");
        }
        
        // After drawing, find the actual last non-empty bin across all drawn histograms and trim
        double actualMaxStepLenEdge = 0.0;
        TList* prims3 = c3->GetListOfPrimitives();
        if (prims3) {
            TObject* obj = prims3->First();
            while (obj) {
                if (obj->InheritsFrom("TH1")) {
                    TH1* h = (TH1*)obj;
                    if (TString(h->GetName()).Contains("StepLength")) {
                        const int lastBin = h->FindLastBinAbove(0.0);
                        if (lastBin > 0) {
                            const double edge = h->GetXaxis()->GetBinUpEdge(lastBin);
                            if (edge > actualMaxStepLenEdge) {
                                actualMaxStepLenEdge = edge;
                            }
                        }
                    }
                }
                obj = prims3->After(obj);
            }
            
            // Set range to trim after the last non-empty bin (with 5% padding)
            if (actualMaxStepLenEdge > 0.0) {
                const double paddedMax = actualMaxStepLenEdge * 1.05;
                // Find maximum histogram range
                obj = prims3->First();
                double maxHistRange = 0.0;
                while (obj) {
                    if (obj->InheritsFrom("TH1")) {
                        TH1* h = (TH1*)obj;
                        if (TString(h->GetName()).Contains("StepLength")) {
                            const double histMax = h->GetXaxis()->GetXmax();
                            if (histMax > maxHistRange) {
                                maxHistRange = histMax;
                            }
                        }
                    }
                    obj = prims3->After(obj);
                }
                const double finalMax = std::min(paddedMax, maxHistRange);
                
                // Set range on all histograms
                obj = prims3->First();
                while (obj) {
                    if (obj->InheritsFrom("TH1")) {
                        TH1* h = (TH1*)obj;
                        if (TString(h->GetName()).Contains("StepLength")) {
                            h->GetXaxis()->SetRangeUser(0.0, finalMax);
                        }
                    }
                    obj = prims3->After(obj);
                }
                
                // Force canvas update to apply range changes
                c3->Modified();
                c3->Update();
            }
        }
        
        leg3->Draw();
        c3->Update();

        // Create a zoomed version of StepLengthModels (0-10nm range)
        TString canvasNameStepLenZoom = Form("StepLengthModelsZoom_E%lld", kv.first);
        TCanvas* c3zoom = new TCanvas(canvasNameStepLenZoom, "Step length models (zoom 0-10nm)", 900, 700);
        c3zoom->SetGrid();
        c3zoom->SetLogy();
        c3zoom->SetLeftMargin(0.15);

        TLegend* leg3zoom = new TLegend(0.45, 0.65, 0.78, 0.88);
        leg3zoom->SetBorderSize(0);
        leg3zoom->SetFillStyle(0);
        leg3zoom->SetTextSize(0.03);

        // Reuse the same histograms but with zoomed range
        bool firstStepLenZoom = true;
        for (const auto& entry : sorted) {
            TFile* f = TFile::Open(entry.path.c_str(), "READ");
            if (!f || f->IsZombie()) {
                continue;
            }
            TH1* h = (TH1*)f->Get("StepLengthAl2O3");
            if (!h) {
                f->Close();
                continue;
            }
            TString cloneName = Form("StepLengthZoom_%s_E%lld", entry.model.c_str(), kv.first);
            TH1* hc = (TH1*)h->Clone(cloneName);
            if (!hc) {
                f->Close();
                continue;
            }
            hc->SetDirectory(nullptr);
            f->Close();
            RebinToMaxBins(hc, 2000);
            int color = modelColor.count(entry.model) ? modelColor[entry.model] : kBlack;
            hc->SetLineColor(color);
            hc->SetLineWidth(4);
            
            // Set zoomed range (0-10nm)
            hc->GetXaxis()->SetRangeUser(0.0, 10.0);
            
            if (firstStepLenZoom) {
                hc->SetTitle(Form("Step length in Al_{2}O_{3} - Zoom 0-10nm (%s)", titleSuffix.Data()));
                hc->GetYaxis()->SetTitleOffset(1.35);
                hc->GetYaxis()->SetLabelSize(0.035);
                hc->Draw("HIST");
                firstStepLenZoom = false;
            } else {
                hc->Draw("HIST SAME");
            }
            leg3zoom->AddEntry(hc, entry.model.c_str(), "l");
        }
        
        leg3zoom->Draw();
        c3zoom->Update();

        out->cd();
        c1->Write(canvasName, TObject::kOverwrite);
        c2->Write(canvasNameSteps, TObject::kOverwrite);
        c3->Write(canvasNameStepLen, TObject::kOverwrite);
        c3zoom->Write(canvasNameStepLenZoom, TObject::kOverwrite);
        if (c1zoom) {
            c1zoom->Write(canvasNameEdepZoom, TObject::kOverwrite);
            saveCanvas(c1zoom, canvasNameEdepZoom);
        }
        saveCanvas(c1, canvasName);
        saveCanvas(c2, canvasNameSteps);
        saveCanvas(c3, canvasNameStepLen);
        saveCanvas(c3zoom, canvasNameStepLenZoom);
    }

    out->Write("", TObject::kOverwrite);
    out->Close();
}
