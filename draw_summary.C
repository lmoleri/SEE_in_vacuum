// ROOT macro to build summary plots across a scan directory.
// Usage:
//   root -l -b -q 'draw_summary.C("results/scan_...","summary.root")'

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TString.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMultiGraph.h"

#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cmath>

namespace {

std::vector<TString> ListRootFiles(const char* dirPath) {
    std::vector<TString> files;
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
            files.push_back(TString(dirPath) + "/" + name);
        }
    }
    files.shrink_to_fit();
    return files;
}

std::vector<int> BuildColors() {
    return {
        kMagenta + 1, kAzure + 2, kOrange + 7, kGreen + 2, kRed + 1,
        kViolet + 6, kCyan + 2, kPink + 2, kTeal + 3, kSpring + 4
    };
}

double ParseParam(const TString& value) {
    if (value.Length() == 0) {
        return -1.0;
    }
    TString tmp = value;
    tmp.ReplaceAll("p", ".");
    return tmp.Atof();
}

TString BaseNameNoExt(const TString& path) {
    TString name = path;
    Ssiz_t slash = name.Last('/');
    if (slash != kNPOS) {
        name = name(slash + 1, name.Length() - slash - 1);
    }
    if (name.EndsWith(".root")) {
        name = name(0, name.Length() - 5);
    }
    return name;
}

TString ExtractBetween(const TString& src, const TString& startKey, const TString& endKey) {
    Ssiz_t start = src.Index(startKey);
    if (start == kNPOS) {
        return "";
    }
    start += startKey.Length();
    if (endKey.Length() == 0) {
        return src(start, src.Length() - start);
    }
    Ssiz_t end = src.Index(endKey, start);
    if (end == kNPOS) {
        return "";
    }
    return src(start, end - start);
}

TString BuildLabelFromName(const TString& baseName) {
    TString thickness = ExtractBetween(baseName, "thick", "nm");
    TString energy = ExtractBetween(baseName, "energy", "MeV");
    TString label;
    if (thickness.Length() > 0) {
        label += "t=" + thickness + " nm";
    }
    if (energy.Length() > 0) {
        if (label.Length() > 0) {
            label += ", ";
        }
        label += "E=" + energy + " MeV";
    }
    if (label.Length() == 0) {
        label = baseName;
    }
    return label;
}

TString EnergyLabel(double energy) {
    if (energy <= 0.) {
        return "E=n/a";
    }
    if (energy < 1e-3) {
        return Form("E=%.0f eV", energy * 1.0e6);
    }
    if (energy < 1.0) {
        return Form("E=%.3g keV", energy * 1.0e3);
    }
    return Form("E=%.3g MeV", energy);
}

TString ParticleLabel(const TString& name) {
    if (name == "e-" || name == "e") return "e^{-}";
    if (name == "e+") return "e^{+}";
    if (name == "mu-" || name == "mu") return "#mu^{-}";
    if (name == "mu+") return "#mu^{+}";
    if (!name.IsNull()) return name;
    return "e^{-}";
}

} // namespace

void draw_summary(const char* scanDir,
                  const char* outputFile = "summary.root") {
    if (!scanDir || TString(scanDir).Length() == 0) {
        printf("Error: scan directory not provided.\n");
        return;
    }

    auto files = ListRootFiles(scanDir);
    if (files.empty()) {
        printf("No .root files found in %s\n", scanDir);
        return;
    }

    gStyle->SetOptStat(0);

    auto colors = BuildColors();
    const int nColors = static_cast<int>(colors.size());

    struct FileEntry {
        TString path;
        double thickness;
        double energy;
        double events;
        TString particle;
    };
    std::vector<FileEntry> entries;
    entries.reserve(files.size());
    for (const auto& filePath : files) {
        TString baseName = BaseNameNoExt(filePath);
        double thickness = ParseParam(ExtractBetween(baseName, "thick", "nm"));
        double energy = ParseParam(ExtractBetween(baseName, "energy", "MeV"));
        double events = ParseParam(ExtractBetween(baseName, "events", ""));
        TString particle = ExtractBetween(baseName, "particle", "_energy");
        entries.push_back({filePath, thickness, energy, events, particle});
    }

    for (auto& entry : entries) {
        TFile* f = TFile::Open(entry.path, "READ");
        if (!f || f->IsZombie()) {
            continue;
        }
        TTree* meta = dynamic_cast<TTree*>(f->Get("RunMeta"));
        if (meta) {
            char particle[64] = "";
            if (meta->GetBranch("primaryParticle")) {
                meta->SetBranchAddress("primaryParticle", particle);
            }
            if (meta->GetEntries() > 0) {
                meta->GetEntry(0);
                if (particle[0] != '\0') {
                    entry.particle = particle;
                }
            }
        }
        f->Close();
    }
    std::sort(entries.begin(), entries.end(), [](const FileEntry& a, const FileEntry& b) {
        if (a.thickness != b.thickness) {
            return a.thickness < b.thickness;
        }
        if (a.energy != b.energy) {
            return a.energy < b.energy;
        }
        return a.events < b.events;
    });

    std::map<TString, std::vector<FileEntry>> entriesByParticle;
    for (const auto& entry : entries) {
        TString key = entry.particle;
        if (key.IsNull()) {
            key = "unknown";
        }
        entriesByParticle[key].push_back(entry);
    }

    for (const auto& particleGroup : entriesByParticle) {
        const TString particleKey = particleGroup.first;
        const auto& particleEntries = particleGroup.second;
        if (particleEntries.empty()) {
            continue;
        }

        TString particleText = ParticleLabel(particleKey);
        TString eventsTitle;
        if (particleEntries.front().events > 0.) {
            eventsTitle = Form("N=%.0f", particleEntries.front().events);
        }

        TString particleTag = particleKey;
        particleTag.ReplaceAll("+", "plus");
        particleTag.ReplaceAll("-", "minus");
        particleTag.ReplaceAll("/", "_");

        TString outPath = outputFile;
        if (outPath.EndsWith(".root")) {
            outPath = outPath(0, outPath.Length() - 5);
        }
        outPath += "_" + particleTag + ".root";

        TFile* out = TFile::Open(outPath, "RECREATE");
        if (!out || out->IsZombie()) {
            printf("Error: cannot create output file %s\n", outPath.Data());
            continue;
        }

        TCanvas* c1 = new TCanvas("EdepPrimarySummary_" + particleTag, "EdepPrimary Summary", 900, 700);
        c1->SetGrid();
        c1->SetLogy();
        c1->SetLeftMargin(0.15);

        TCanvas* c2 = new TCanvas("EdepInteractionsSummary_" + particleTag, "EdepInteractions Summary", 900, 700);
        c2->SetGrid();
        c2->SetLeftMargin(0.15);

    TLegend* leg1 = new TLegend(0.45, 0.65, 0.78, 0.88);
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
    leg1->SetTextSize(0.03);
        if (!eventsTitle.IsNull()) {
            leg1->SetHeader(eventsTitle, "C");
        }

    TLegend* leg2 = new TLegend(0.45, 0.65, 0.78, 0.88);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
    leg2->SetTextSize(0.03);
        if (!eventsTitle.IsNull()) {
            leg2->SetHeader(eventsTitle, "C");
        }

        bool firstPrimary = true;
        bool firstSteps = true;
        int colorIndex = 0;
        std::map<double, std::vector<std::pair<double, double>>> fracByEnergy;
        std::map<double, std::vector<std::pair<double, double>>> mpvByEnergy;

        double maxPrimaryEnergyMeV = 0.0;
        for (const auto& entry : particleEntries) {
            if (entry.energy > maxPrimaryEnergyMeV) {
                maxPrimaryEnergyMeV = entry.energy;
            }
        }

        for (const auto& entry : particleEntries) {
        const auto& filePath = entry.path;
        TFile* f = TFile::Open(filePath, "READ");
        if (!f || f->IsZombie()) {
            printf("Warning: cannot open %s\n", filePath.Data());
            continue;
        }

        TH1* hPrimary = dynamic_cast<TH1*>(f->Get("EdepPrimary"));
        TH1* hSteps = dynamic_cast<TH1*>(f->Get("EdepInteractions"));
        if (!hPrimary || !hSteps) {
            printf("Warning: missing histograms in %s\n", filePath.Data());
            f->Close();
            continue;
        }

        TString baseName = BaseNameNoExt(filePath);
        TString label = BuildLabelFromName(baseName);
        int color = colors[colorIndex % nColors];
        colorIndex++;

        TH1* hPrimaryClone = dynamic_cast<TH1*>(hPrimary->Clone("primary_" + label));
        TH1* hStepsClone = dynamic_cast<TH1*>(hSteps->Clone("steps_" + label));
        if (!hPrimaryClone || !hStepsClone) {
            f->Close();
            continue;
        }
        hPrimaryClone->SetDirectory(nullptr);
        hStepsClone->SetDirectory(nullptr);
        hPrimaryClone->SetLineColor(color);
        hStepsClone->SetLineColor(color);
        hPrimaryClone->SetLineWidth(5);
        hStepsClone->SetLineWidth(4);

        // Fraction of events with non-zero interactions (steps > 0)
        const int stepsBins = hSteps->GetNbinsX();
        const double totalStepsEvents = hSteps->Integral(1, stepsBins);
        const double zeroStepsEvents = hSteps->GetBinContent(1);
        const double fracNonZero =
            (totalStepsEvents > 0.) ? (totalStepsEvents - zeroStepsEvents) / totalStepsEvents : 0.;
        fracByEnergy[entry.energy].push_back({entry.thickness, fracNonZero});

        // MPV of Landau fit to energy deposition histogram
        double mpv = std::nan("");
        if (hPrimary->GetEntries() > 0) {
            TF1 landauFit("landauFit", "landau", hPrimary->GetXaxis()->GetXmin(),
                          hPrimary->GetXaxis()->GetXmax());
            int fitStatus = hPrimary->Fit(&landauFit, "Q0");
            if (fitStatus == 0) {
                mpv = landauFit.GetParameter(1);
            }
        }
        if (!std::isnan(mpv)) {
            mpvByEnergy[entry.energy].push_back({entry.thickness, mpv});
        }

        c1->cd();
        if (firstPrimary) {
            hPrimaryClone->SetTitle(Form("EdepPrimary Summary (%s)", particleText.Data()));
            hPrimaryClone->GetYaxis()->SetTitleOffset(1.35);
            hPrimaryClone->GetYaxis()->SetLabelSize(0.035);
            if (maxPrimaryEnergyMeV > 0.) {
                hPrimaryClone->GetXaxis()->SetRangeUser(0.0, maxPrimaryEnergyMeV * 1.0e6);
            }
            hPrimaryClone->Draw("HIST");
            firstPrimary = false;
        } else {
            hPrimaryClone->Draw("HIST SAME");
        }
        leg1->AddEntry(hPrimaryClone, label, "l");

        c2->cd();
        if (firstSteps) {
            hStepsClone->SetTitle(Form("EdepInteractions Summary (%s)", particleText.Data()));
            hStepsClone->GetYaxis()->SetTitleOffset(1.35);
            hStepsClone->GetYaxis()->SetLabelSize(0.035);
            hStepsClone->Draw("HIST");
            firstSteps = false;
        } else {
            hStepsClone->Draw("HIST SAME");
        }
        leg2->AddEntry(hStepsClone, label, "l");

        f->Close();
    }

        c1->cd();
        leg1->Draw();
        c1->Update();

        c2->cd();
        leg2->Draw();
        c2->Update();

        // Summary graph: fraction of events with non-zero steps
        TCanvas* c3 = new TCanvas("NonZeroStepsSummary_" + particleTag, "Non-zero steps fraction", 900, 700);
        c3->SetGrid();
        c3->SetLeftMargin(0.15);
        TMultiGraph* mgFrac = new TMultiGraph();
        TLegend* leg3 = new TLegend(0.45, 0.65, 0.78, 0.88);
        leg3->SetBorderSize(0);
        leg3->SetFillStyle(0);
        leg3->SetTextSize(0.03);
        if (!eventsTitle.IsNull()) {
            leg3->SetHeader(eventsTitle, "C");
        }

        int fracColorIndex = 0;
        for (const auto& kv : fracByEnergy) {
            const auto& points = kv.second;
            if (points.empty()) {
                continue;
            }
            TGraph* g = new TGraph(static_cast<int>(points.size()));
            for (int i = 0; i < static_cast<int>(points.size()); ++i) {
                g->SetPoint(i, points[i].first, points[i].second);
            }
            int color = colors[fracColorIndex % nColors];
            fracColorIndex++;
            g->SetLineColor(color);
            g->SetMarkerColor(color);
            g->SetMarkerStyle(20);
            g->SetLineWidth(3);
            mgFrac->Add(g, "LP");
            leg3->AddEntry(g, EnergyLabel(kv.first), "lp");
        }
        mgFrac->SetTitle(Form("Non-zero steps fraction (%s)", particleText.Data()));
        mgFrac->Draw("A");
        mgFrac->GetXaxis()->SetTitle("Sample thickness (nm)");
        mgFrac->GetYaxis()->SetTitle("Fraction of events with non-zero steps");
        mgFrac->GetYaxis()->SetTitleOffset(1.35);
        mgFrac->GetYaxis()->SetLabelSize(0.035);
        mgFrac->GetYaxis()->SetRangeUser(0.0, 1.0);
        leg3->Draw();
        c3->Update();

        // Summary graph: MPV of Landau fit to EdepPrimary
        TCanvas* c4 = new TCanvas("EdepPrimaryMPVSummary_" + particleTag, "EdepPrimary MPV (Landau)", 900, 700);
        c4->SetGrid();
        c4->SetLeftMargin(0.15);
        TMultiGraph* mgMpv = new TMultiGraph();
        TLegend* leg4 = new TLegend(0.45, 0.65, 0.78, 0.88);
        leg4->SetBorderSize(0);
        leg4->SetFillStyle(0);
        leg4->SetTextSize(0.03);
        if (!eventsTitle.IsNull()) {
            leg4->SetHeader(eventsTitle, "C");
        }

        int mpvColorIndex = 0;
        for (const auto& kv : mpvByEnergy) {
            const auto& points = kv.second;
            if (points.empty()) {
                continue;
            }
            TGraph* g = new TGraph(static_cast<int>(points.size()));
            for (int i = 0; i < static_cast<int>(points.size()); ++i) {
                g->SetPoint(i, points[i].first, points[i].second);
            }
            int color = colors[mpvColorIndex % nColors];
            mpvColorIndex++;
            g->SetLineColor(color);
            g->SetMarkerColor(color);
            g->SetMarkerStyle(21);
            g->SetLineWidth(3);
            mgMpv->Add(g, "LP");
            leg4->AddEntry(g, EnergyLabel(kv.first), "lp");
        }
        mgMpv->SetTitle(Form("EdepPrimary MPV (Landau) (%s)", particleText.Data()));
        mgMpv->Draw("A");
        mgMpv->GetXaxis()->SetTitle("Sample thickness (nm)");
        mgMpv->GetYaxis()->SetTitle("MPV of Landau fit (eV)");
        mgMpv->GetYaxis()->SetTitleOffset(1.35);
        mgMpv->GetYaxis()->SetLabelSize(0.035);
        leg4->Draw();
        c4->Update();

        out->cd();
        c1->Write("EdepPrimarySummary_" + particleTag, TObject::kOverwrite);
        c2->Write("EdepInteractionsSummary_" + particleTag, TObject::kOverwrite);
        c3->Write("NonZeroStepsSummary_" + particleTag, TObject::kOverwrite);
        c4->Write("EdepPrimaryMPVSummary_" + particleTag, TObject::kOverwrite);
        out->Write("", TObject::kOverwrite);
        out->Close();
    }
}
