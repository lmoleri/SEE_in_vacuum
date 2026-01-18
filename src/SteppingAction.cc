#include "SteppingAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4VProcess.hh"

SteppingAction::SteppingAction(RunAction* runAction, EventAction* eventAction)
    : G4UserSteppingAction(),
      fRunAction(runAction),
      fEventAction(eventAction)
{
}

SteppingAction::~SteppingAction()
{
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    if (!fRunAction) return;

    const auto* prePoint  = step->GetPreStepPoint();
    const auto* postPoint = step->GetPostStepPoint();

    const auto* preVol  = prePoint->GetPhysicalVolume();
    const auto* postVol = postPoint->GetPhysicalVolume();

    if (!preVol || !postVol) {
        return;
    }

    const G4String preName  = preVol->GetName();
    const G4String postName = postVol->GetName();

    // Count any energy-depositing step in Al2O3 (all particles)
    if (fEventAction && preName == "Al2O3") {
        const G4double edep = step->GetTotalEnergyDeposit();
        if (edep > 0.) {
            fEventAction->AddEdepInteraction();
            auto* analysisManager = G4AnalysisManager::Instance();
            if (analysisManager) {
                analysisManager->FillH1(1, edep / eV);
            }
        }
    }

    const G4Track* track    = step->GetTrack();
    const auto*    particle = track->GetParticleDefinition();

    // We are interested only in electrons for the metrics below
    if (particle->GetParticleName() != "e-") {
        return;
    }

    // 1) Accumulate primary electron energy deposition in Al2O3
    //    (trackID==1, parentID==0, pre-step in Al2O3)
    if (fEventAction && track->GetTrackID() == 1 && track->GetParentID() == 0 &&
        preName == "Al2O3") {
        const G4double edep = step->GetTotalEnergyDeposit();
        if (edep > 0.) {
            fEventAction->AddPrimaryEdep(edep);
        }
    }

    // 1b) Capture per-step energy transfer in Al2O3 when PAI is enabled.
    //     Use energy deposition on electron steps as proxy for microscopic transfers.
    if (preName == "Al2O3" && fRunAction && fRunAction->IsPaiEnabled()) {
        const G4String& pname = particle->GetParticleName();
        if (pname == "e-") {
            const G4double edep = step->GetTotalEnergyDeposit();
            if (edep > 0.) {
                auto* analysisManager = G4AnalysisManager::Instance();
                if (analysisManager) {
                    analysisManager->FillH1(3, edep / eV);
                }
            }
        }
    }

    // 2) Count secondary electrons that EXIT the Al2O3 layer into vacuum (world)
    //    and exclude the primary beam electron (parentID == 0).
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

