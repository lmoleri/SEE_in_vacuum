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

} // namespace

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
            TString cloneName = Form("EdepPrimary_%s_E%lld", entry.model.c_str(), kv.first);
            TH1* hc = (TH1*)h->Clone(cloneName);
            if (!hc) {
                f->Close();
                continue;
            }
            hc->SetDirectory(nullptr);
            f->Close();
            int color = modelColor.count(entry.model) ? modelColor[entry.model] : kBlack;
            hc->SetLineColor(color);
            hc->SetLineWidth(4);
            if (first) {
                hc->SetTitle(Form("EdepPrimary model comparison (%s)", titleSuffix.Data()));
                hc->GetYaxis()->SetTitleOffset(1.35);
                hc->GetYaxis()->SetLabelSize(0.035);
                if (energyMeV > 0.) {
                    hc->GetXaxis()->SetRangeUser(0.0, energyMeV * 1.0e6);
                }
                hc->Draw("HIST");
                first = false;
            } else {
                hc->Draw("HIST SAME");
            }
            leg1->AddEntry(hc, entry.model.c_str(), "l");
        }
        leg1->Draw();
        c1->Update();

        TString canvasNameSteps = Form("EdepInteractionsModels_E%lld", kv.first);
        TCanvas* c2 = new TCanvas(canvasNameSteps, "EdepInteractions models", 900, 700);
        c2->SetGrid();
        c2->SetLeftMargin(0.15);

        TLegend* leg2 = new TLegend(0.45, 0.65, 0.78, 0.88);
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

        out->cd();
        c1->Write(canvasName, TObject::kOverwrite);
        c2->Write(canvasNameSteps, TObject::kOverwrite);
    }

    out->Write("", TObject::kOverwrite);
    out->Close();
}
