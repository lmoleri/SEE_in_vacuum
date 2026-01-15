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
    
    // Read run metadata (if present)
    double primaryEnergyMeV = -1.0;
    double sampleThicknessNm = -1.0;
    char primaryParticle[64] = "";
    TTree* meta = (TTree*)f->Get("RunMeta");
    if (meta) {
        meta->SetBranchAddress("primaryEnergyMeV", &primaryEnergyMeV);
        meta->SetBranchAddress("sampleThicknessNm", &sampleThicknessNm);
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
    std::string particleText = particleLabel(primaryParticle);

    // Create a canvas
    TCanvas* c1 = new TCanvas("c1", "Energy Deposition", 800, 600);
    c1->SetGrid();
    c1->SetLogy();
    c1->SetLeftMargin(0.15);
    
    // Fold overflow into the last visible bin so the tail is visible
    int nbins = h->GetNbinsX();
    h->SetBinContent(nbins, h->GetBinContent(nbins) + h->GetBinContent(nbins + 1));
    h->SetBinError(nbins, std::sqrt(std::pow(h->GetBinError(nbins), 2) +
                                    std::pow(h->GetBinError(nbins + 1), 2)));
    h->SetBinContent(nbins + 1, 0.0);
    h->SetBinError(nbins + 1, 0.0);

    // Format Y axis to avoid label overlap and use scientific notation
    h->GetYaxis()->SetTitleOffset(1.35);
    h->GetYaxis()->SetLabelSize(0.035);

    // Style: magenta line with thicker width
    h->SetLineColor(kMagenta + 1);
    h->SetLineWidth(3);

    // Draw histogram with the normal default style
    h->SetTitle(Form("Primary %s energy deposition in Al_{2}O_{3}", particleText.c_str()));
    h->Draw("HIST");

    // Add info text
    TLatex info1;
    info1.SetNDC(true);
    info1.SetTextSize(0.03);
    if (primaryEnergyMeV > 0.) {
        info1.DrawLatex(0.45, 0.85,
                        Form("Primary %s energy: %.3f MeV", particleText.c_str(), primaryEnergyMeV));
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
                        Form("Primary %s energy: %.3f MeV", particleText.c_str(), primaryEnergyMeV));
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

    f->Write("", TObject::kOverwrite);
    f->Close();
}
