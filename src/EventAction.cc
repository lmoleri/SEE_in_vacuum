#include "EventAction.hh"

#include "RunAction.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"

EventAction::EventAction(RunAction* runAction)
    : G4UserEventAction(),
      fRunAction(runAction),
      fEdepPrimary(0.),
      fNMicroscopicEdep(0)
{
}

EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction(const G4Event*)
{
    // Reset per-event accumulator
    fEdepPrimary = 0.;
    fNMicroscopicEdep = 0;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    auto* analysisManager = G4AnalysisManager::Instance();
    if (!analysisManager) return;

    // Histogram ID 0: primary e- energy deposition in Al2O3
    // Convert from internal units (keV) to eV for display
    analysisManager->FillH1(0, fEdepPrimary / eV);

    // Histogram ID 1: number of microscopic energy-depositing steps per event
    analysisManager->FillH1(1, static_cast<G4double>(fNMicroscopicEdep));

    if (fRunAction && fEdepPrimary > 0.) {
        fRunAction->UpdateMinNonZeroEdep(fEdepPrimary);
    }
}

void EventAction::AddPrimaryEdep(G4double edep)
{
    fEdepPrimary += edep;
}

void EventAction::AddEdepInteraction()
{
    ++fNMicroscopicEdep;
}

