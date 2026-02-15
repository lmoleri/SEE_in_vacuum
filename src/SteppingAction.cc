#include "SteppingAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
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

    if (!preVol) {
        return;
    }

    const G4String preName  = preVol->GetName();
    const G4String postName = postVol ? postVol->GetName() : "OutOfWorld";

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
            const G4int stepLenId = fRunAction->GetStepLengthAl2O3Id();
            if (analysisManager && stepLenId >= 0) {
                analysisManager->FillH1(stepLenId, stepLen / nm);
            }
        }
    }

    const G4Track* track    = step->GetTrack();
    const auto*    particle = track->GetParticleDefinition();
    const G4String particleName = particle->GetParticleName();
    const G4bool isPrimary = (track->GetTrackID() == 1 && track->GetParentID() == 0);
    const auto* process = postPoint->GetProcessDefinedStep();

    const auto computeDepthNm = [&](G4double zCoord) -> G4double {
        const G4double thickness = fRunAction->GetSampleThickness();
        if (thickness <= 0.) {
            return 0.0;
        }
        const G4double zmin = -0.5 * thickness;
        const G4double zmax = 0.5 * thickness;
        const G4double dirZ = fRunAction->GetPrimaryDirectionZ();
        G4double depth = (dirZ >= 0.) ? (zCoord - zmin) : (zmax - zCoord);
        if (depth < 0.) depth = 0.;
        if (depth > thickness) depth = thickness;
        return depth / nm;
    };

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
        if (particleName == "e-") {
            const G4double preDepthNm = computeDepthNm(prePoint->GetPosition().z());
            const G4double postDepthNm = computeDepthNm(postPoint->GetPosition().z());
            const G4double preE = prePoint->GetKineticEnergy();
            const G4double postE = postPoint->GetKineticEnergy();
            const G4double stepLenNm = step->GetStepLength() / nm;
            const G4int stepStatus = static_cast<G4int>(postPoint->GetStepStatus());
            const G4String procName = process ? process->GetProcessName() : "None";
            const auto preDir = prePoint->GetMomentumDirection();
            const auto postDir = postPoint->GetMomentumDirection();
            auto signZ = [](G4double z) -> G4int {
                const G4double threshold = 1e-9;
                if (z > threshold) return 1;
                if (z < -threshold) return -1;
                return 0;
            };
            G4double dot = preDir.x() * postDir.x() + preDir.y() * postDir.y() + preDir.z() * postDir.z();
            if (dot > 1.0) dot = 1.0;
            if (dot < -1.0) dot = -1.0;
            const G4double deltaThetaDeg = std::acos(dot) / deg;
            const G4int reversalOnStep =
                (signZ(preDir.z()) != 0 && signZ(postDir.z()) != 0 &&
                 signZ(preDir.z()) != signZ(postDir.z()))
                    ? 1
                    : 0;

            fEventAction->UpdatePrimaryMaxDepthNm(preDepthNm);
            fEventAction->UpdatePrimaryMaxDepthNm(postDepthNm);
            fEventAction->AddPrimaryStepAudit(procName, stepStatus);
            fEventAction->UpdatePrimaryDirectionSignZ(track->GetMomentumDirection().z(),
                                                      track->GetCurrentStepNumber(), preDepthNm, postE,
                                                      procName, stepStatus, stepLenNm, preE, postE,
                                                      deltaThetaDeg);
            fEventAction->RecordPrimaryTrajectoryStep(track->GetCurrentStepNumber(),
                                                      preDepthNm, postDepthNm, stepLenNm, preE,
                                                      postE, edep, preDir.z(), postDir.z(),
                                                      deltaThetaDeg, reversalOnStep, procName,
                                                      stepStatus, preName, postName);
            fEventAction->UpdatePrimaryFirstProcessInAl2O3(
                procName);
            if (edep > 0.) {
                fEventAction->AddPrimaryEdepByProcess(
                    procName, edep, preDepthNm);
            }
        }
    }

    if (fEventAction && isPrimary && particleName == "e-") {
        const G4bool preInAl2O3 = (preName == "Al2O3");
        const G4bool postInAl2O3 = (postName == "Al2O3");
        if (preInAl2O3 != postInAl2O3) {
            G4double boundaryDepthNm = 0.0;
            if (preInAl2O3) {
                boundaryDepthNm = computeDepthNm(prePoint->GetPosition().z());
            } else if (postInAl2O3) {
                boundaryDepthNm = computeDepthNm(postPoint->GetPosition().z());
            }
            fEventAction->RecordPrimaryBoundaryCrossing(track->GetCurrentStepNumber(),
                                                        boundaryDepthNm,
                                                        postPoint->GetKineticEnergy(), preName,
                                                        postName);
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
        if (process) {
            fEventAction->UpdatePrimaryLastProcess(process->GetProcessName());
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
            const G4double depthNm = computeDepthNm(prePoint->GetPosition().z());
            const G4double alphaInvNm = fRunAction->GetSeyAlphaInvNm();
            auto* analysisManager = G4AnalysisManager::Instance();
            if (alphaInvNm > 0.) {
                const G4double weight = std::exp(-alphaInvNm * depthNm);
                fEventAction->AddPrimaryEdepWeighted(edep * weight);
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
                    const G4int stepDepthId = fRunAction->GetEdepStepDepthPrimaryId();
                    if (stepDepthId >= 0) {
                        analysisManager->FillH2(stepDepthId, depthNm, edep / eV);
                    }
                }
            }
            if (fRunAction->IsVerboseStepDiagnostics()) {
                const G4double preE = prePoint->GetKineticEnergy();
                const G4double postE = postPoint->GetKineticEnergy();
                const G4double frac = (preE > 0.) ? (edep / preE) : 0.0;
                if (frac >= fRunAction->GetVerboseStepThresholdFrac()) {
                    const auto slot = fRunAction->ConsumeVerboseStepSlot();
                    if (slot >= 0) {
                        const auto* proc = postPoint->GetProcessDefinedStep();
                        const auto stepStatus = postPoint->GetStepStatus();
                        const auto trackStatus = track->GetTrackStatus();
                        const auto ntupleId = fRunAction->GetVerboseStepNtupleId();
                        auto volumeCode = [](const G4String& name) -> G4int {
                            if (name == "Al2O3") return 1;
                            if (name == "World") return 2;
                            if (name == "Si") return 3;
                            if (name == "OutOfWorld") return 4;
                            return 0;
                        };
                        auto processCode = [](const G4String& name) -> G4int {
                            if (name == "msc") return 1;
                            if (name == "eIoni") return 2;
                            if (name == "eBrem") return 3;
                            if (name == "Transportation") return 4;
                            if (name == "eCoulombScattering") return 5;
                            return 0;
                        };
                        if (analysisManager && ntupleId >= 0) {
                            const auto* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
                            const auto eventId = evt ? evt->GetEventID() : -1;
                            analysisManager->FillNtupleIColumn(ntupleId, 0, eventId);
                            analysisManager->FillNtupleIColumn(ntupleId, 1, track->GetTrackID());
                            analysisManager->FillNtupleIColumn(ntupleId, 2, track->GetCurrentStepNumber());
                            analysisManager->FillNtupleDColumn(ntupleId, 3, depthNm);
                            analysisManager->FillNtupleDColumn(ntupleId, 4, step->GetStepLength() / nm);
                            analysisManager->FillNtupleDColumn(ntupleId, 5, edep / eV);
                            analysisManager->FillNtupleDColumn(ntupleId, 6, preE / eV);
                            analysisManager->FillNtupleDColumn(ntupleId, 7, postE / eV);
                            analysisManager->FillNtupleDColumn(ntupleId, 8, frac);
                            analysisManager->FillNtupleIColumn(ntupleId, 9, static_cast<G4int>(stepStatus));
                            analysisManager->FillNtupleIColumn(ntupleId, 10, static_cast<G4int>(trackStatus));
                            analysisManager->FillNtupleIColumn(ntupleId, 11, volumeCode(preName));
                            analysisManager->FillNtupleIColumn(ntupleId, 12, volumeCode(postName));
                            analysisManager->FillNtupleIColumn(ntupleId, 13,
                                                               proc ? processCode(proc->GetProcessName())
                                                                    : 0);
                            analysisManager->FillNtupleSColumn(ntupleId, 14, preName);
                            analysisManager->FillNtupleSColumn(ntupleId, 15, postName);
                            analysisManager->FillNtupleSColumn(ntupleId, 16,
                                                               proc ? proc->GetProcessName() : "None");
                            analysisManager->AddNtupleRow(ntupleId);
                        }
                        G4cout << "[VerboseStep] slot=" << slot
                               << " edep=" << edep / eV << " eV"
                               << " preE=" << preE / eV << " eV"
                               << " frac=" << frac
                               << " depth=" << depthNm << " nm"
                               << " stepLen=" << step->GetStepLength() / nm << " nm"
                               << " preVol=" << preName
                               << " postVol=" << postName
                               << " proc=" << (proc ? proc->GetProcessName() : "None")
                               << " stepStatus=" << stepStatus
                               << " trackStatus=" << trackStatus
                               << G4endl;
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

        // Diagnostics for primary electron exits: split by entrance/opposite/lateral side.
        if (isPrimary && particleName == "e-") {
            if (fEventAction && fRunAction) {
                const G4double thickness = fRunAction->GetSampleThickness();
                const G4double zMin = -0.5 * thickness;
                const G4double zMax = 0.5 * thickness;
                const G4double zExit = postPoint->GetPosition().z();
                const G4double tol = 0.2 * nm;
                const G4bool nearMin = (std::abs(zExit - zMin) <= tol);
                const G4bool nearMax = (std::abs(zExit - zMax) <= tol);
                const G4bool entranceAtMin = (fRunAction->GetPrimaryDirectionZ() >= 0.);

                // 1=entrance-side exit, 2=opposite-side exit, 3=lateral/edge exit
                G4int exitClass = 3;
                if (nearMin || nearMax) {
                    const G4bool entranceExit =
                        (nearMin && entranceAtMin) || (nearMax && !entranceAtMin);
                    exitClass = entranceExit ? 1 : 2;
                }

                // Keep only the most recent Al2O3->World crossing for this event.
                // Final per-event filling is done in EventAction::EndOfEventAction.
                fEventAction->UpdatePrimaryExitCandidate(
                    exitClass, eKin, postPoint->GetMomentumDirection());
            }
        }

        if (eKin > eThreshold) {
            fRunAction->AddEmittedElectron();
            if (parentId > 0) {
                fRunAction->AddSecondaryElectron();
            }
        }
    }
}
