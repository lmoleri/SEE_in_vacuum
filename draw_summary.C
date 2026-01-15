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

#include <vector>
#include <string>

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
        if (name.EndsWith(".root")) {
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

    TFile* out = TFile::Open(outputFile, "RECREATE");
    if (!out || out->IsZombie()) {
        printf("Error: cannot create output file %s\n", outputFile);
        return;
    }

    gStyle->SetOptStat(0);

    auto colors = BuildColors();
    const int nColors = static_cast<int>(colors.size());

    TCanvas* c1 = new TCanvas("EdepPrimarySummary", "EdepPrimary Summary", 900, 700);
    c1->SetGrid();
    c1->SetLogy();
    c1->SetLeftMargin(0.15);

    TCanvas* c2 = new TCanvas("EdepInteractionsSummary", "EdepInteractions Summary", 900, 700);
    c2->SetGrid();
    c2->SetLeftMargin(0.15);

    TLegend* leg1 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);

    TLegend* leg2 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    bool firstPrimary = true;
    bool firstSteps = true;

    int colorIndex = 0;
    for (const auto& filePath : files) {
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

        TString label = BaseNameNoExt(filePath);
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
        hPrimaryClone->SetLineWidth(2);
        hStepsClone->SetLineWidth(2);

        c1->cd();
        if (firstPrimary) {
            hPrimaryClone->GetYaxis()->SetTitleOffset(1.35);
            hPrimaryClone->GetYaxis()->SetLabelSize(0.035);
            hPrimaryClone->Draw("HIST");
            firstPrimary = false;
        } else {
            hPrimaryClone->Draw("HIST SAME");
        }
        leg1->AddEntry(hPrimaryClone, label, "l");

        c2->cd();
        if (firstSteps) {
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

    out->cd();
    c1->Write("EdepPrimarySummary", TObject::kOverwrite);
    c2->Write("EdepInteractionsSummary", TObject::kOverwrite);
    out->Write("", TObject::kOverwrite);
    out->Close();
}
