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

    // Record step length in Al2O3 for diagnostics
    if (preName == "Al2O3") {
        const G4double stepLen = step->GetStepLength();
        if (stepLen > 0.) {
            auto* analysisManager = G4AnalysisManager::Instance();
            if (analysisManager) {
                analysisManager->FillH1(6, stepLen / nm);
            }
        }
    }

    const G4Track* track    = step->GetTrack();
    const auto*    particle = track->GetParticleDefinition();
    const G4String particleName = particle->GetParticleName();
    const G4bool isPrimary = (track->GetTrackID() == 1 && track->GetParentID() == 0);

    // Accumulate primary particle energy deposition in Al2O3 (for all primary particles)
    // (trackID==1, parentID==0, pre-step in Al2O3)
    if (fEventAction && isPrimary && preName == "Al2O3") {
        const G4double edep = step->GetTotalEnergyDeposit();
        if (edep > 0.) {
            fEventAction->AddPrimaryEdep(edep);
        }
    }

    // Track primary particle residual energy and other metrics (for all primary particles)
    if (fEventAction && isPrimary) {
        fEventAction->UpdatePrimaryResidualEnergy(postPoint->GetKineticEnergy());
        const auto* postVol = postPoint->GetPhysicalVolume();
        if (postVol) {
            fEventAction->UpdatePrimaryLastVolume(postVol->GetName());
        } else {
            fEventAction->UpdatePrimaryLastVolume("OutOfWorld");
        }
        const auto* proc = postPoint->GetProcessDefinedStep();
        if (proc) {
            fEventAction->UpdatePrimaryLastProcess(proc->GetProcessName());
        } else {
            fEventAction->UpdatePrimaryLastProcess("Unknown");
        }
        fEventAction->UpdatePrimaryStopStatus(static_cast<G4int>(track->GetTrackStatus()));
        if (track->GetTrackStatus() == fStopAndKill) {
            const auto* endVol = postPoint->GetPhysicalVolume();
            if (endVol) {
                fEventAction->UpdatePrimaryEndVolume(endVol->GetName());
            } else {
                fEventAction->UpdatePrimaryEndVolume("OutOfWorld");
            }
        }
    }

    // Electron-specific tracking below
    if (particleName != "e-") {
        return;
    }

    // 1b) Capture per-step energy transfer in Al2O3 when PAI is enabled.
    //     Use energy deposition on electron steps as proxy for microscopic transfers.
    if (preName == "Al2O3" && fRunAction && fRunAction->IsPaiEnabled()) {
        const G4double edep = step->GetTotalEnergyDeposit();
        if (edep > 0.) {
            auto* analysisManager = G4AnalysisManager::Instance();
            if (analysisManager) {
                analysisManager->FillH1(3, edep / eV);
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

