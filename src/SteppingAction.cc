#include "SteppingAction.hh"

#include "RunAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction(RunAction* runAction)
    : G4UserSteppingAction(),
      fRunAction(runAction)
{
}

SteppingAction::~SteppingAction()
{
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    if (!fRunAction) return;

    const G4Track* track = step->GetTrack();
    const auto* particle = track->GetParticleDefinition();

    // We are interested only in electrons
    if (particle->GetParticleName() != "e-") {
        return;
    }

    // Count secondary electrons that EXIT the Al2O3 layer into vacuum (world)
    const auto* prePoint  = step->GetPreStepPoint();
    const auto* postPoint = step->GetPostStepPoint();

    const auto* preVol  = prePoint->GetPhysicalVolume();
    const auto* postVol = postPoint->GetPhysicalVolume();

    if (!preVol || !postVol) {
        return;
    }

    const G4String preName  = preVol->GetName();
    const G4String postName = postVol->GetName();

    // Only consider transitions from the Al2O3 target to the world
    // and exclude the primary beam electron (parentID == 0).
    if (preName == "Al2O3" && postName == "World") {
        const G4int parentId = track->GetParentID();
        if (parentId > 0) {
            // Apply a small kinetic energy threshold to avoid counting
            // essentially "stopped" electrons as secondaries.
            const G4double eKin = postPoint->GetKineticEnergy();
            const G4double eThreshold = 1.0 * eV; // adjustable

            if (eKin > eThreshold) {
                fRunAction->AddSecondaryElectron();
            }
        }
    }
}

