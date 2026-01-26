// ROOT macro to draw the energy deposition histogram
// Usage: root -l draw_histo.C
// Or from build directory: root -l ../draw_histo.C
// Optional: root -l 'draw_histo.C("file.root")'

#include "TSystem.h"

#include <algorithm>
#include <cctype>
#include <string>

void draw_histo(const char* fileName = "SEE_in_vacuum.root") {
    // Open the ROOT file
    TFile* f = TFile::Open(fileName, "UPDATE");
    if (!f || f->IsZombie()) {
        cout << "Error: Cannot open SEE_in_vacuum.root" << endl;
        cout << "Make sure you're running from the build directory" << endl;
        return;
    }
    
    // Get the histogram
    TH1* h = (TH1*)f->Get("EdepPrimary");
    if (!h) {
        cout << "Error: Histogram EdepPrimary not found" << endl;
        f->Close();
        return;
    }

    TH1* hSteps = (TH1*)f->Get("EdepInteractions");
    if (!hSteps) {
        cout << "Error: Histogram EdepInteractions not found" << endl;
        f->Close();
        return;
    }

    TH1* hStep = (TH1*)f->Get("EdepStep");
    if (!hStep) {
        cout << "Warning: Histogram EdepStep not found" << endl;
    }

    TH1* hPai = (TH1*)f->Get("PAITransfer");
    if (!hPai) {
        cout << "Warning: Histogram PAITransfer not found" << endl;
    }

    TH1* hResidual = (TH1*)f->Get("PrimaryResidualEnergy");
    if (!hResidual) {
        cout << "Warning: Histogram PrimaryResidualEnergy not found" << endl;
    }

    TH1* hEndVol = (TH1*)f->Get("PrimaryEndVolume");
    if (!hEndVol) {
        cout << "Warning: Histogram PrimaryEndVolume not found" << endl;
    }

    TH1* hStepLen = (TH1*)f->Get("StepLengthAl2O3");
    if (!hStepLen) {
        cout << "Warning: Histogram StepLengthAl2O3 not found" << endl;
    }

    TH2* hEdepVsSteps = (TH2*)f->Get("EdepPrimaryVsSteps");
    if (!hEdepVsSteps) {
        cout << "Warning: Histogram EdepPrimaryVsSteps not found" << endl;
    }

    TH2* hResVsEnd = (TH2*)f->Get("ResidualEnergyVsEndVolume");
    if (!hResVsEnd) {
        cout << "Warning: Histogram ResidualEnergyVsEndVolume not found" << endl;
    }

    TH2* hResVsProc = (TH2*)f->Get("ResidualEnergyVsLastProcess");
    if (!hResVsProc) {
        cout << "Warning: Histogram ResidualEnergyVsLastProcess not found" << endl;
    }

    TH2* hResVsStop = (TH2*)f->Get("ResidualEnergyVsStopStatus");
    if (!hResVsStop) {
        cout << "Warning: Histogram ResidualEnergyVsStopStatus not found" << endl;
    }
    
    // Read run metadata (if present)
    double primaryEnergyMeV = -1.0;
    double sampleThicknessNm = -1.0;
    double maxPrimaryEnergyMeV = -1.0;
    char primaryParticle[64] = "";
    char emModel[64] = "";
    TTree* meta = (TTree*)f->Get("RunMeta");
    if (meta) {
        meta->SetBranchAddress("primaryEnergyMeV", &primaryEnergyMeV);
        meta->SetBranchAddress("sampleThicknessNm", &sampleThicknessNm);
        if (meta->GetBranch("maxPrimaryEnergyMeV")) {
            meta->SetBranchAddress("maxPrimaryEnergyMeV", &maxPrimaryEnergyMeV);
        }
        if (meta->GetBranch("primaryParticle")) {
            meta->SetBranchAddress("primaryParticle", primaryParticle);
        }
        if (meta->GetBranch("emModel")) {
            meta->SetBranchAddress("emModel", emModel);
        }
        meta->GetEntry(0);
    }

    auto particleLabel = [](const std::string& name) {
        if (name == "e-" || name == "e") return std::string("e^{-}");
        if (name == "e+") return std::string("e^{+}");
        if (name == "mu-" || name == "mu") return std::string("#mu^{-}");
        if (name == "mu+") return std::string("#mu^{+}");
        if (!name.empty()) return name;
        return std::string("e^{-}");
    };
    auto modelLabel = [](std::string name) {
        std::transform(name.begin(), name.end(), name.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        if (name == "pai") return std::string("PAI");
        if (name == "livermore" || name == "livermorephysics" ||
            name == "g4emlivermorephysics") return std::string("Livermore");
        if (name == "penelope" || name == "penelopephysics" ||
            name == "g4empenelopephysics") return std::string("Penelope");
        if (!name.empty()) return name;
        return std::string("unknown");
    };
    auto energyLabel = [](double energyMeV) -> std::string {
        if (energyMeV <= 0.) return std::string("n/a");
        if (energyMeV < 1e-3) {
            return std::string(Form("%.0f eV", energyMeV * 1.0e6));
        }
        if (energyMeV < 1.0) {
            return std::string(Form("%.3g keV", energyMeV * 1.0e3));
        }
        return std::string(Form("%.3g MeV", energyMeV));
    };
    std::string particleText = particleLabel(primaryParticle);
    std::string modelText = modelLabel(emModel);
    auto trimXAxis = [](TH1* hist) {
        if (!hist) return;
        const int lastBin = hist->FindLastBinAbove(0.0);
        if (lastBin > 0) {
            hist->GetXaxis()->SetRangeUser(0.0, hist->GetXaxis()->GetBinUpEdge(lastBin));
        }
    };
    auto trimAxes = [](TH2* hist) {
        if (!hist) return;
        int maxX = 0;
        int maxY = 0;
        const int nx = hist->GetNbinsX();
        const int ny = hist->GetNbinsY();
        for (int ix = 1; ix <= nx; ++ix) {
            for (int iy = 1; iy <= ny; ++iy) {
                if (hist->GetBinContent(ix, iy) > 0.0) {
                    if (ix > maxX) maxX = ix;
                    if (iy > maxY) maxY = iy;
                }
            }
        }
        if (maxX > 0) {
            hist->GetXaxis()->SetRangeUser(0.0, hist->GetXaxis()->GetBinUpEdge(maxX));
        }
        if (maxY > 0) {
            hist->GetYaxis()->SetRangeUser(0.0, hist->GetYaxis()->GetBinUpEdge(maxY));
        }
    };
    std::string plotDir;
    {
        std::string path = fileName ? fileName : "";
        std::string baseDir;
        std::string subDir;
        const std::string marker = "/results/";
        const auto pos = path.rfind(marker);
        if (pos != std::string::npos) {
            baseDir = path.substr(0, pos);
            std::string rest = path.substr(pos + marker.size());
            const auto slash = rest.rfind('/');
            if (slash != std::string::npos) {
                subDir = rest.substr(0, slash);
            }
        } else {
            std::string dir = gSystem->DirName(path.c_str());
            baseDir = gSystem->DirName(dir.c_str());
            subDir = gSystem->BaseName(dir.c_str());
            if (baseDir == "results") {
                baseDir = ".";
            }
        }
        std::string fileBase = gSystem->BaseName(path.c_str());
        if (fileBase.size() >= 5 && fileBase.substr(fileBase.size() - 5) == ".root") {
            fileBase = fileBase.substr(0, fileBase.size() - 5);
        }
        if (baseDir.empty() || baseDir == ".") {
            plotDir = "plots";
        } else {
            plotDir = baseDir + "/plots";
        }
        if (!subDir.empty()) {
            plotDir += "/" + subDir;
        }
        plotDir += "/" + fileBase;
        gSystem->mkdir(plotDir.c_str(), true);
    }
    auto saveCanvas = [&](TCanvas* canvas, const char* name) {
        if (!canvas || !name || plotDir.empty()) {
            return;
        }
        std::string rootPath = plotDir + "/" + name + ".root";
        std::string pdfPath = plotDir + "/" + name + ".pdf";
        canvas->SaveAs(rootPath.c_str());
        canvas->SaveAs(pdfPath.c_str());
    };

    // Create a canvas
    TCanvas* c1 = new TCanvas("c1", "Energy Deposition", 800, 600);
    c1->SetGrid();
    c1->SetLogy();
    c1->SetLeftMargin(0.15);
    
    // Format Y axis to avoid label overlap and use scientific notation
    h->GetYaxis()->SetTitleOffset(1.35);
    h->GetYaxis()->SetLabelSize(0.035);

    // Style: magenta line with thicker width
    h->SetLineColor(kMagenta + 1);
    h->SetLineWidth(3);

    // Draw histogram with the normal default style
    h->SetTitle(Form("Primary %s energy deposition in Al_{2}O_{3} (%s)",
                     particleText.c_str(), modelText.c_str()));
    trimXAxis(h);
    h->Draw("HIST");

    // Add info text
    TLatex info1;
    info1.SetNDC(true);
    info1.SetTextSize(0.03);
    if (primaryEnergyMeV > 0.) {
        info1.DrawLatex(0.45, 0.85,
                        Form("Primary %s energy: %s", particleText.c_str(),
                             energyLabel(primaryEnergyMeV).c_str()));
    } else {
        info1.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
    }
    if (sampleThicknessNm > 0.) {
        info1.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
    } else {
        info1.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
    }
    info1.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));
    
    // Update and write canvas into the ROOT file (so TBrowser can show it)
    c1->Update();
    c1->Draw();
    f->cd();
    c1->Write("EdepPrimaryCanvas", TObject::kOverwrite);
    saveCanvas(c1, "EdepPrimaryCanvas");

    // Create a canvas for the interaction count histogram
    TCanvas* c2 = new TCanvas("c2", "Energy Deposition Steps", 800, 600);
    c2->SetGrid();
    c2->SetLogy();
    c2->SetLeftMargin(0.15);

    // Format Y axis to avoid label overlap and use scientific notation
    hSteps->GetYaxis()->SetTitleOffset(1.35);
    hSteps->GetYaxis()->SetLabelSize(0.035);

    // Style: magenta line with thicker width
    hSteps->SetLineColor(kMagenta + 1);
    hSteps->SetLineWidth(3);

    hSteps->SetTitle(Form("Energy-depositing steps per event (%s, %s)",
                          particleText.c_str(), modelText.c_str()));
    trimXAxis(hSteps);
    hSteps->Draw("HIST");

    // Add info text
    TLatex info2;
    info2.SetNDC(true);
    info2.SetTextSize(0.03);
    if (primaryEnergyMeV > 0.) {
        info2.DrawLatex(0.45, 0.85,
                        Form("Primary %s energy: %s", particleText.c_str(),
                             energyLabel(primaryEnergyMeV).c_str()));
    } else {
        info2.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
    }
    if (sampleThicknessNm > 0.) {
        info2.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
    } else {
        info2.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
    }
    info2.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));

    c2->Update();
    c2->Draw();
    f->cd();
    c2->Write("EdepInteractionsCanvas", TObject::kOverwrite);
    saveCanvas(c2, "EdepInteractionsCanvas");

    if (hStep) {
        TCanvas* c3 = new TCanvas("c3", "Energy Deposition per Step", 800, 600);
        c3->SetGrid();
        c3->SetLeftMargin(0.15);

        hStep->GetYaxis()->SetTitleOffset(1.35);
        hStep->GetYaxis()->SetLabelSize(0.035);
        hStep->SetLineColor(kMagenta + 1);
        hStep->SetLineWidth(3);
        hStep->SetTitle(Form("Energy deposition per step (%s, %s)",
                             particleText.c_str(), modelText.c_str()));
        trimXAxis(hStep);
        hStep->Draw("HIST");

        TLatex info3;
        info3.SetNDC(true);
        info3.SetTextSize(0.03);
        if (primaryEnergyMeV > 0.) {
            info3.DrawLatex(0.45, 0.85,
                            Form("Primary %s energy: %s", particleText.c_str(),
                                 energyLabel(primaryEnergyMeV).c_str()));
        } else {
            info3.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
        }
        if (sampleThicknessNm > 0.) {
            info3.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
        } else {
            info3.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
        }
        info3.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));

        c3->Update();
        c3->Draw();
        f->cd();
        c3->Write("EdepStepCanvas", TObject::kOverwrite);
        saveCanvas(c3, "EdepStepCanvas");
    }

    if (hPai) {
        TCanvas* cPai = new TCanvas("cPai", "PAI Energy Transfer per Step", 800, 600);
        cPai->SetGrid();
        cPai->SetLeftMargin(0.15);

        hPai->GetYaxis()->SetTitleOffset(1.35);
        hPai->GetYaxis()->SetLabelSize(0.035);
        hPai->SetLineColor(kMagenta + 1);
        hPai->SetLineWidth(3);
        hPai->SetTitle(Form("PAI energy transfer per step (%s, %s)",
                            particleText.c_str(), modelText.c_str()));
        hPai->Draw("HIST");

        TLatex infoPai;
        infoPai.SetNDC(true);
        infoPai.SetTextSize(0.03);
        if (primaryEnergyMeV > 0.) {
            infoPai.DrawLatex(0.45, 0.85,
                              Form("Primary %s energy: %s", particleText.c_str(),
                                   energyLabel(primaryEnergyMeV).c_str()));
        } else {
            infoPai.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
        }
        if (sampleThicknessNm > 0.) {
            infoPai.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
        } else {
            infoPai.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
        }
        infoPai.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));

        cPai->Update();
        cPai->Draw();
        f->cd();
        cPai->Write("PAITransferCanvas", TObject::kOverwrite);
        saveCanvas(cPai, "PAITransferCanvas");
    }

    if (hResidual) {
        TCanvas* cRes = new TCanvas("cRes", "Primary Residual Energy", 800, 600);
        cRes->SetGrid();
        cRes->SetLogy();
        cRes->SetLeftMargin(0.15);

        hResidual->GetYaxis()->SetTitleOffset(1.35);
        hResidual->GetYaxis()->SetLabelSize(0.035);
        hResidual->SetLineColor(kMagenta + 1);
        hResidual->SetLineWidth(3);
        hResidual->SetTitle(Form("Primary residual energy (%s, %s)",
                                 particleText.c_str(), modelText.c_str()));
        trimXAxis(hResidual);
        hResidual->Draw("HIST");

        TLatex infoRes;
        infoRes.SetNDC(true);
        infoRes.SetTextSize(0.03);
        if (primaryEnergyMeV > 0.) {
            infoRes.DrawLatex(0.45, 0.85,
                              Form("Primary %s energy: %s", particleText.c_str(),
                                   energyLabel(primaryEnergyMeV).c_str()));
        } else {
            infoRes.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
        }
        if (sampleThicknessNm > 0.) {
            infoRes.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
        } else {
            infoRes.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
        }
        infoRes.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));

        cRes->Update();
        cRes->Draw();
        f->cd();
        cRes->Write("PrimaryResidualEnergyCanvas", TObject::kOverwrite);
        saveCanvas(cRes, "PrimaryResidualEnergyCanvas");
    }

    if (hEndVol) {
        TCanvas* cEnd = new TCanvas("cEnd", "Primary End Volume", 800, 600);
        cEnd->SetGrid();
        cEnd->SetLogy();
        cEnd->SetLeftMargin(0.15);

        hEndVol->GetYaxis()->SetTitleOffset(1.35);
        hEndVol->GetYaxis()->SetLabelSize(0.035);
        hEndVol->SetLineColor(kMagenta + 1);
        hEndVol->SetLineWidth(3);
        hEndVol->SetTitle(Form("Primary end volume (%s, %s)",
                               particleText.c_str(), modelText.c_str()));
        hEndVol->GetXaxis()->SetBinLabel(1, "Unknown");
        hEndVol->GetXaxis()->SetBinLabel(2, "Al2O3");
        hEndVol->GetXaxis()->SetBinLabel(3, "World");
        hEndVol->GetXaxis()->SetBinLabel(4, "OutOfWorld");
        hEndVol->GetXaxis()->SetBinLabel(5, "Other");
        hEndVol->LabelsOption("v");
        hEndVol->Draw("HIST");

        TLatex infoEnd;
        infoEnd.SetNDC(true);
        infoEnd.SetTextSize(0.03);
        if (primaryEnergyMeV > 0.) {
            infoEnd.DrawLatex(0.45, 0.85,
                              Form("Primary %s energy: %s", particleText.c_str(),
                                   energyLabel(primaryEnergyMeV).c_str()));
        } else {
            infoEnd.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
        }
        if (sampleThicknessNm > 0.) {
            infoEnd.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
        } else {
            infoEnd.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
        }
        infoEnd.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));

        cEnd->Update();
        cEnd->Draw();
        f->cd();
        cEnd->Write("PrimaryEndVolumeCanvas", TObject::kOverwrite);
        saveCanvas(cEnd, "PrimaryEndVolumeCanvas");
    }

    if (hStepLen) {
        TCanvas* cStep = new TCanvas("cStep", "Step Length in Al2O3", 800, 600);
        cStep->SetGrid();
        cStep->SetLogy();
        cStep->SetLeftMargin(0.15);

        hStepLen->GetYaxis()->SetTitleOffset(1.35);
        hStepLen->GetYaxis()->SetLabelSize(0.035);
        hStepLen->SetLineColor(kMagenta + 1);
        hStepLen->SetLineWidth(3);
        hStepLen->SetTitle(Form("Step length in Al_{2}O_{3} (%s, %s)",
                                particleText.c_str(), modelText.c_str()));
        trimXAxis(hStepLen);
        hStepLen->Draw("HIST");

        TLatex infoStep;
        infoStep.SetNDC(true);
        infoStep.SetTextSize(0.03);
        if (primaryEnergyMeV > 0.) {
            infoStep.DrawLatex(0.45, 0.85,
                               Form("Primary %s energy: %s", particleText.c_str(),
                                    energyLabel(primaryEnergyMeV).c_str()));
        } else {
            infoStep.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
        }
        if (sampleThicknessNm > 0.) {
            infoStep.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
        } else {
            infoStep.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
        }
        infoStep.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));

        cStep->Update();
        cStep->Draw();
        f->cd();
        cStep->Write("StepLengthAl2O3Canvas", TObject::kOverwrite);
        saveCanvas(cStep, "StepLengthAl2O3Canvas");
    }

    if (hEdepVsSteps) {
        TCanvas* c4 = new TCanvas("c4", "EdepPrimary vs Steps", 800, 600);
        c4->SetGrid();
        c4->SetLeftMargin(0.15);
        c4->SetRightMargin(0.15);

        hEdepVsSteps->SetTitle(Form("Primary edep vs steps (%s, %s)",
                                    particleText.c_str(), modelText.c_str()));
        hEdepVsSteps->GetYaxis()->SetTitleOffset(1.35);
        hEdepVsSteps->GetYaxis()->SetLabelSize(0.035);
        trimAxes(hEdepVsSteps);
        hEdepVsSteps->Draw("COLZ");

        TLatex info4;
        info4.SetNDC(true);
        info4.SetTextSize(0.03);
        if (primaryEnergyMeV > 0.) {
            info4.DrawLatex(0.45, 0.85,
                            Form("Primary %s energy: %s", particleText.c_str(),
                                 energyLabel(primaryEnergyMeV).c_str()));
        } else {
            info4.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
        }
        if (sampleThicknessNm > 0.) {
            info4.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
        } else {
            info4.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
        }
        info4.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));

        c4->Update();
        c4->Draw();
        f->cd();
        c4->Write("EdepPrimaryVsStepsCanvas", TObject::kOverwrite);
        saveCanvas(c4, "EdepPrimaryVsStepsCanvas");
    }

    if (hResVsEnd) {
        TCanvas* c5 = new TCanvas("c5", "Residual Energy vs End Volume", 800, 600);
        c5->SetGrid();
        c5->SetLeftMargin(0.15);
        c5->SetRightMargin(0.18);

        hResVsEnd->SetTitle(Form("Residual energy vs end volume (%s, %s)",
                                 particleText.c_str(), modelText.c_str()));
        hResVsEnd->GetYaxis()->SetTitleOffset(1.35);
        hResVsEnd->GetYaxis()->SetLabelSize(0.035);
        hResVsEnd->GetYaxis()->SetBinLabel(1, "Unknown");
        hResVsEnd->GetYaxis()->SetBinLabel(2, "Al2O3");
        hResVsEnd->GetYaxis()->SetBinLabel(3, "World");
        hResVsEnd->GetYaxis()->SetBinLabel(4, "OutOfWorld");
        hResVsEnd->GetYaxis()->SetBinLabel(5, "Other");
        hResVsEnd->GetYaxis()->LabelsOption("v");
        hResVsEnd->Draw("COLZ");

        TLatex info5;
        info5.SetNDC(true);
        info5.SetTextSize(0.03);
        if (primaryEnergyMeV > 0.) {
            info5.DrawLatex(0.45, 0.85,
                            Form("Primary %s energy: %s", particleText.c_str(),
                                 energyLabel(primaryEnergyMeV).c_str()));
        } else {
            info5.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
        }
        if (sampleThicknessNm > 0.) {
            info5.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
        } else {
            info5.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
        }
        info5.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));

        c5->Update();
        c5->Draw();
        f->cd();
        c5->Write("ResidualEnergyVsEndVolumeCanvas", TObject::kOverwrite);
        saveCanvas(c5, "ResidualEnergyVsEndVolumeCanvas");
    }

    if (hResVsProc) {
        TCanvas* c6 = new TCanvas("c6", "Residual Energy vs Last Process", 800, 600);
        c6->SetGrid();
        c6->SetLeftMargin(0.15);
        c6->SetRightMargin(0.18);

        hResVsProc->SetTitle(Form("Residual energy vs last process (%s, %s)",
                                  particleText.c_str(), modelText.c_str()));
        hResVsProc->GetYaxis()->SetTitleOffset(1.35);
        hResVsProc->GetYaxis()->SetLabelSize(0.035);
        hResVsProc->GetYaxis()->SetBinLabel(1, "Unknown");
        hResVsProc->GetYaxis()->SetBinLabel(2, "Transportation");
        hResVsProc->GetYaxis()->SetBinLabel(3, "eIoni");
        hResVsProc->GetYaxis()->SetBinLabel(4, "msc");
        hResVsProc->GetYaxis()->SetBinLabel(5, "eBrem");
        hResVsProc->GetYaxis()->SetBinLabel(6, "CoulombScat");
        hResVsProc->GetYaxis()->SetBinLabel(7, "Other");
        hResVsProc->GetYaxis()->LabelsOption("v");
        hResVsProc->Draw("COLZ");

        TLatex info6;
        info6.SetNDC(true);
        info6.SetTextSize(0.03);
        if (primaryEnergyMeV > 0.) {
            info6.DrawLatex(0.45, 0.85,
                            Form("Primary %s energy: %s", particleText.c_str(),
                                 energyLabel(primaryEnergyMeV).c_str()));
        } else {
            info6.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
        }
        if (sampleThicknessNm > 0.) {
            info6.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
        } else {
            info6.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
        }
        info6.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));

        c6->Update();
        c6->Draw();
        f->cd();
        c6->Write("ResidualEnergyVsLastProcessCanvas", TObject::kOverwrite);
        saveCanvas(c6, "ResidualEnergyVsLastProcessCanvas");
    }

    if (hResVsStop) {
        TCanvas* c7 = new TCanvas("c7", "Residual Energy vs Stop Status", 800, 600);
        c7->SetGrid();
        c7->SetLeftMargin(0.15);
        c7->SetRightMargin(0.18);

        hResVsStop->SetTitle(Form("Residual energy vs stop status (%s, %s)",
                                  particleText.c_str(), modelText.c_str()));
        hResVsStop->GetYaxis()->SetTitleOffset(1.35);
        hResVsStop->GetYaxis()->SetLabelSize(0.035);
        hResVsStop->GetYaxis()->SetBinLabel(1, "Unknown");
        hResVsStop->GetYaxis()->SetBinLabel(2, "fAlive");
        hResVsStop->GetYaxis()->SetBinLabel(3, "fStopButAlive");
        hResVsStop->GetYaxis()->SetBinLabel(4, "fStopAndKill");
        hResVsStop->GetYaxis()->SetBinLabel(5, "fKillTrackAndSecondaries");
        hResVsStop->GetYaxis()->SetBinLabel(6, "fSuspend");
        hResVsStop->GetYaxis()->SetBinLabel(7, "fPostponeToNextEvent");
        hResVsStop->GetYaxis()->LabelsOption("v");
        hResVsStop->Draw("COLZ");

        TLatex info7;
        info7.SetNDC(true);
        info7.SetTextSize(0.03);
        if (primaryEnergyMeV > 0.) {
            info7.DrawLatex(0.45, 0.85,
                            Form("Primary %s energy: %s", particleText.c_str(),
                                 energyLabel(primaryEnergyMeV).c_str()));
        } else {
            info7.DrawLatex(0.45, 0.85, Form("Primary %s energy: n/a", particleText.c_str()));
        }
        if (sampleThicknessNm > 0.) {
            info7.DrawLatex(0.45, 0.79, Form("Sample thickness: %.2f nm", sampleThicknessNm));
        } else {
            info7.DrawLatex(0.45, 0.79, "Sample thickness: n/a");
        }
        info7.DrawLatex(0.45, 0.73, Form("EM model: %s", modelText.c_str()));

        c7->Update();
        c7->Draw();
        f->cd();
        c7->Write("ResidualEnergyVsStopStatusCanvas", TObject::kOverwrite);
        saveCanvas(c7, "ResidualEnergyVsStopStatusCanvas");
    }

    f->Close();
}
