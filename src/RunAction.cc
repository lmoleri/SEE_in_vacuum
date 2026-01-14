#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

RunAction::RunAction()
    : G4UserRunAction(),
      fNPrimaryElectrons(0),
      fNSecondaryElectrons(0)
{
}

RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    fNPrimaryElectrons = 0;
    fNSecondaryElectrons = 0;
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
}

void RunAction::AddPrimaryElectron()
{
    // No-op for now; primary count is taken from the number of events.
}

void RunAction::AddSecondaryElectron()
{
    ++fNSecondaryElectrons;
}

