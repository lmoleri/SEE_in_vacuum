// ROOT macro to draw the energy deposition histogram
// Usage: root -l draw_histo.C
// Or from build directory: root -l ../draw_histo.C
// Optional: root -l 'draw_histo.C("file.root")'

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

    TH2* hEdepVsSteps = (TH2*)f->Get("EdepPrimaryVsSteps");
    if (!hEdepVsSteps) {
        cout << "Warning: Histogram EdepPrimaryVsSteps not found" << endl;
    }
    
    // Read run metadata (if present)
    double primaryEnergyMeV = -1.0;
    double sampleThicknessNm = -1.0;
    double maxPrimaryEnergyMeV = -1.0;
    char primaryParticle[64] = "";
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
    h->SetTitle(Form("Primary %s energy deposition in Al_{2}O_{3}", particleText.c_str()));
    const double xMaxEv = h->GetXaxis()->GetXmax();
    if (xMaxEv > 0.) {
        h->GetXaxis()->SetRangeUser(0.0, xMaxEv);
    }
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
    
    // Update and write canvas into the ROOT file (so TBrowser can show it)
    c1->Update();
    c1->Draw();
    f->cd();
    c1->Write("EdepPrimaryCanvas", TObject::kOverwrite);

    // Create a canvas for the interaction count histogram
    TCanvas* c2 = new TCanvas("c2", "Energy Deposition Steps", 800, 600);
    c2->SetGrid();
    c2->SetLeftMargin(0.15);

    // Format Y axis to avoid label overlap and use scientific notation
    hSteps->GetYaxis()->SetTitleOffset(1.35);
    hSteps->GetYaxis()->SetLabelSize(0.035);

    // Style: magenta line with thicker width
    hSteps->SetLineColor(kMagenta + 1);
    hSteps->SetLineWidth(3);

    hSteps->SetTitle(Form("Energy-depositing steps per event (%s)", particleText.c_str()));
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

    c2->Update();
    c2->Draw();
    f->cd();
    c2->Write("EdepInteractionsCanvas", TObject::kOverwrite);

    if (hStep) {
        TCanvas* c3 = new TCanvas("c3", "Energy Deposition per Step", 800, 600);
        c3->SetGrid();
        c3->SetLeftMargin(0.15);

        hStep->GetYaxis()->SetTitleOffset(1.35);
        hStep->GetYaxis()->SetLabelSize(0.035);
        hStep->SetLineColor(kMagenta + 1);
        hStep->SetLineWidth(3);
        hStep->SetTitle(Form("Energy deposition per step (%s)", particleText.c_str()));
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

        c3->Update();
        c3->Draw();
        f->cd();
        c3->Write("EdepStepCanvas", TObject::kOverwrite);
    }

    if (hPai) {
        TCanvas* cPai = new TCanvas("cPai", "PAI Energy Transfer per Step", 800, 600);
        cPai->SetGrid();
        cPai->SetLeftMargin(0.15);

        hPai->GetYaxis()->SetTitleOffset(1.35);
        hPai->GetYaxis()->SetLabelSize(0.035);
        hPai->SetLineColor(kMagenta + 1);
        hPai->SetLineWidth(3);
        hPai->SetTitle(Form("PAI energy transfer per step (%s)", particleText.c_str()));
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

        cPai->Update();
        cPai->Draw();
        f->cd();
        cPai->Write("PAITransferCanvas", TObject::kOverwrite);
    }

    if (hEdepVsSteps) {
        TCanvas* c4 = new TCanvas("c4", "EdepPrimary vs Steps", 800, 600);
        c4->SetGrid();
        c4->SetLeftMargin(0.15);
        c4->SetRightMargin(0.15);

        hEdepVsSteps->SetTitle(Form("Primary edep vs steps (%s)", particleText.c_str()));
        hEdepVsSteps->GetYaxis()->SetTitleOffset(1.35);
        hEdepVsSteps->GetYaxis()->SetLabelSize(0.035);
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

        c4->Update();
        c4->Draw();
        f->cd();
        c4->Write("EdepPrimaryVsStepsCanvas", TObject::kOverwrite);
    }

    f->Close();
}
