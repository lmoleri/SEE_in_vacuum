// ROOT macro to draw the energy deposition histogram
// Usage: root -l draw_histo.C
// Or from build directory: root -l ../draw_histo.C

void draw_histo() {
    // Open the ROOT file
    TFile* f = TFile::Open("SEE_in_vacuum.root");
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
    
    // Create a canvas
    TCanvas* c1 = new TCanvas("c1", "Energy Deposition", 800, 600);
    c1->SetGrid();
    
    // Draw histogram with the normal default style
    h->Draw("HIST");
    
    // Update the canvas
    c1->Update();
    c1->Draw();
}
