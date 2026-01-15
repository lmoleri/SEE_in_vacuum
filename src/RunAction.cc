#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"

#include "DetectorConstruction.hh"

RunAction::RunAction()
    : G4UserRunAction(),
      fNPrimaryElectrons(0),
      fNSecondaryElectrons(0),
      fMinNonZeroEdep(-1.),
      fPrimaryEnergy(0.),
      fPrimaryParticleName("e-"),
      fSampleThickness(0.),
      fOutputTag("SEE_in_vacuum")
{
}

RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    fNPrimaryElectrons   = 0;
    fNSecondaryElectrons = 0;
    fMinNonZeroEdep = -1.;

    auto* analysisManager = G4AnalysisManager::Instance();

    // Basic configuration
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFirstHistoId(0);

    // Capture sample thickness from detector if not set
    if (fSampleThickness <= 0.) {
        auto* det = dynamic_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        if (det) {
            fSampleThickness = det->GetSampleThickness();
        }
    }

    static G4bool metaCreated = false;
    static G4bool histosCreated = false;
    if (!metaCreated) {
        analysisManager->CreateNtuple("RunMeta", "Run metadata");
        analysisManager->CreateNtupleDColumn("primaryEnergyMeV");
        analysisManager->CreateNtupleDColumn("sampleThicknessNm");
        analysisManager->CreateNtupleSColumn("primaryParticle");
        analysisManager->FinishNtuple();
        metaCreated = true;
    }

    if (histosCreated) {
        analysisManager->Reset();
    } else {
        // Create a 1D histogram for primary e- energy deposition in Al2O3
        // ID 0: EdepPrimary
        // Based on typical energy deposition: mean ~6 eV, RMS ~29 eV
        // Use range 0-200 eV with fine binning for better resolution
        // Note: Values will be filled in eV units (converted in EventAction)
        G4int histoId = analysisManager->CreateH1(
            "EdepPrimary",
            "Primary e^{-} energy deposition in Al_{2}O_{3}",
            200,   // number of bins (1 eV per bin for good resolution)
            0.,    // Edep min (eV)
            200.   // Edep max (eV, covers mean + ~7*RMS)
        );

        // Set axis labels explicitly (X-axis in eV units)
        analysisManager->SetH1XAxisTitle(histoId, "Energy deposition (eV)");
        analysisManager->SetH1YAxisTitle(histoId, "Number of events");

        // Create a 1D histogram for the number of microscopic energy-depositing
        // steps in Al2O3 per event
        G4int stepsHistoId = analysisManager->CreateH1(
            "EdepInteractions",
            "Microscopic energy-depositing steps in Al_{2}O_{3} per event",
            10,   // number of bins (1 step per bin)
            0.,   // min count
            10.   // max count
        );

        analysisManager->SetH1XAxisTitle(stepsHistoId, "Energy-depositing steps per event");
        analysisManager->SetH1YAxisTitle(stepsHistoId, "Number of events");

        histosCreated = true;
    }

    // Open output file for this run
    G4String fileName = fOutputTag;
    if (fileName.empty()) {
        fileName = "SEE_in_vacuum";
    }
    if (fileName.size() < 5 || fileName.substr(fileName.size() - 5) != ".root") {
        fileName += ".root";
    }
    analysisManager->OpenFile(fileName);
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    const auto nEvents = run->GetNumberOfEvent();

    // Define "primary electrons" as the number of incident beam particles,
    // i.e. one primary electron per event.
    fNPrimaryElectrons = nEvents;

    G4cout << "\n=== Run summary (SEY) ===" << G4endl;
    G4cout << "  Events processed       : " << nEvents << G4endl;
    G4cout << "  Primary electrons      : " << fNPrimaryElectrons << G4endl;
    G4cout << "  Secondary electrons out: " << fNSecondaryElectrons << G4endl;

    if (fNPrimaryElectrons > 0) {
        const G4double sey = static_cast<G4double>(fNSecondaryElectrons) /
                             static_cast<G4double>(fNPrimaryElectrons);
        G4cout << "  Secondary electron yield (SEY) = "
               << sey << G4endl;
    } else {
        G4cout << "  Secondary electron yield (SEY) not defined (no primaries counted)"
               << G4endl;
    }
    if (fMinNonZeroEdep > 0.) {
        G4cout << "  Min non-zero primary edep : " << fMinNonZeroEdep / eV << " eV"
               << G4endl;
    } else {
        G4cout << "  Min non-zero primary edep : none" << G4endl;
    }
    G4cout << "==========================\n" << G4endl;

    // Write and close analysis output
    auto* analysisManager = G4AnalysisManager::Instance();
    if (analysisManager) {
        analysisManager->FillNtupleDColumn(0, fPrimaryEnergy / MeV);
        analysisManager->FillNtupleDColumn(1, fSampleThickness / nm);
        analysisManager->FillNtupleSColumn(2, fPrimaryParticleName);
        analysisManager->AddNtupleRow();

        analysisManager->Write();
        analysisManager->CloseFile();
    }
}

void RunAction::AddPrimaryElectron()
{
    // No-op for now; primary count is taken from the number of events.
}

void RunAction::AddSecondaryElectron()
{
    ++fNSecondaryElectrons;
}

void RunAction::UpdateMinNonZeroEdep(G4double edep)
{
    if (edep <= 0.) {
        return;
    }
    if (fMinNonZeroEdep < 0. || edep < fMinNonZeroEdep) {
        fMinNonZeroEdep = edep;
    }
}

void RunAction::SetPrimaryEnergy(G4double energy)
{
    fPrimaryEnergy = energy;
}

void RunAction::SetPrimaryParticleName(const G4String& name)
{
    fPrimaryParticleName = name;
}

void RunAction::SetSampleThickness(G4double thickness)
{
    fSampleThickness = thickness;
}

void RunAction::SetOutputTag(const G4String& tag)
{
    fOutputTag = tag;
}

