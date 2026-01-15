#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"

RunAction::RunAction()
    : G4UserRunAction(),
      fNPrimaryElectrons(0),
      fNSecondaryElectrons(0)
{
}

RunAction::~RunAction()
{
    delete G4AnalysisManager::Instance();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    fNPrimaryElectrons   = 0;
    fNSecondaryElectrons = 0;

    auto* analysisManager = G4AnalysisManager::Instance();

    // Basic configuration
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFirstHistoId(0);

    // Open output file
    analysisManager->OpenFile("SEE_in_vacuum.root");

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
    G4cout << "==========================\n" << G4endl;

    // Write and close analysis output
    auto* analysisManager = G4AnalysisManager::Instance();
    if (analysisManager) {
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

