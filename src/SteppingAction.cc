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
#include <cmath>

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
        const G4double stepLen = step->GetStepLength();
        if (stepLen > 0.) {
            fEventAction->AddPrimaryTrackLength(stepLen);
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

    // Depth-weighted primary electron energy deposition in Al2O3.
    // Use only the entrance side (direction of the primary particle).
    if (fEventAction && isPrimary && preName == "Al2O3") {
        const G4double edep = step->GetTotalEnergyDeposit();
        if (edep > 0. && fRunAction) {
            const G4double alphaInvNm = fRunAction->GetSeyAlphaInvNm();
            if (alphaInvNm > 0.) {
                const G4double thickness = fRunAction->GetSampleThickness();
                const G4double zmin = -0.5 * thickness;
                const G4double zmax = 0.5 * thickness;
                const G4double z = prePoint->GetPosition().z();
                const G4double dirZ = fRunAction->GetPrimaryDirectionZ();
                G4double depth = (dirZ >= 0.) ? (z - zmin) : (zmax - z);
                if (depth < 0.) depth = 0.;
                if (depth > thickness) depth = thickness;
                const G4double depthNm = depth / nm;
                const G4double weight = std::exp(-alphaInvNm * depthNm);
                fEventAction->AddPrimaryEdepWeighted(edep * weight);
                auto* analysisManager = G4AnalysisManager::Instance();
                if (analysisManager) {
                    const G4int depthId = fRunAction->GetEdepDepthPrimaryId();
                    if (depthId >= 0) {
                        analysisManager->FillH1(depthId, depthNm, edep / eV);
                    }
                    const G4int depthWeightedId = fRunAction->GetEdepDepthPrimaryWeightedId();
                    if (depthWeightedId >= 0) {
                        analysisManager->FillH1(depthWeightedId, depthNm, edep * weight / eV);
                    }
                    const G4int depthCountsId = fRunAction->GetEdepDepthPrimaryCountsId();
                    if (depthCountsId >= 0) {
                        analysisManager->FillH1(depthCountsId, depthNm, 1.0);
                    }
                }
            }
        }
    }

    // Record primary electron track length vs depth in Al2O3 (path-length weighting).
    if (isPrimary && preName == "Al2O3" && fRunAction) {
        const G4double stepLen = step->GetStepLength();
        if (stepLen > 0.) {
            const G4double thickness = fRunAction->GetSampleThickness();
            const G4double zmin = -0.5 * thickness;
            const G4double zmax = 0.5 * thickness;
            const G4double z = prePoint->GetPosition().z();
            const G4double dirZ = fRunAction->GetPrimaryDirectionZ();
            G4double depth = (dirZ >= 0.) ? (z - zmin) : (zmax - z);
            if (depth < 0.) depth = 0.;
            if (depth > thickness) depth = thickness;
            const G4double depthNm = depth / nm;
            auto* analysisManager = G4AnalysisManager::Instance();
            if (analysisManager) {
                const G4int trackDepthId = fRunAction->GetPrimaryTrackLengthDepthId();
                if (trackDepthId >= 0) {
                    analysisManager->FillH1(trackDepthId, depthNm, stepLen / nm);
                }
            }
        }
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

    // 2) Count electrons that EXIT the Al2O3 layer into vacuum (world).
    //    Track both total emitted electrons (including primaries) and true secondaries.
    if (preName == "Al2O3" && postName == "World") {
        const G4int parentId = track->GetParentID();
        // Apply a small kinetic energy threshold to avoid counting
        // essentially "stopped" electrons as emitted.
        const G4double eKin = postPoint->GetKineticEnergy();
        const G4double eThreshold = 1.0 * eV; // adjustable

        if (eKin > eThreshold) {
            fRunAction->AddEmittedElectron();
            if (parentId > 0) {
                fRunAction->AddSecondaryElectron();
            }
        }
    }
}
