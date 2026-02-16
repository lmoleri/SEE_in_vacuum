#include "EventAction.hh"

#include "RunAction.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4StepStatus.hh"

#include <cmath>

EventAction::EventAction(RunAction* runAction)
    : G4UserEventAction(),
      fRunAction(runAction),
      fEdepPrimary(0.),
      fEdepPrimaryWeighted(0.),
      fPrimaryTrackLength(0.),
      fNMicroscopicEdep(0),
      fPrimaryResidualEnergy(0.),
      fPrimaryEndLocation(0),
      fPrimaryLastLocation(0),
      fPrimaryLastProcess(""),
      fPrimaryStopStatus(-1),
      fHasPrimaryExitCandidate(false),
      fPrimaryExitClassCandidate(0),
      fPrimaryExitEnergyCandidate(0.),
      fPrimaryExitDirectionCandidate(0., 0., 0.),
      fPrimaryEdepByEIoni(0.),
      fPrimaryEdepByMsc(0.),
      fPrimaryEdepByOther(0.),
      fPrimaryFirstStepEdep(0.),
      fPrimaryMaxStepEdep(0.),
      fDepthFirstEdepNm(-1.),
      fPrimaryMaxDepthNm(0.),
      fPrimaryBoundaryCrossings(0),
      fPrimaryDirectionReversals(0),
      fLastPrimaryDirectionSign(0),
      fHasLastPrimaryDirectionSign(false),
      fPrimaryFirstProcessInAl2O3(""),
      fHasFirstDirectionReversal(false),
      fFirstDirectionReversalStep(-1),
      fFirstDirectionReversalDepthNm(-1.),
      fFirstDirectionReversalEnergy(-1.),
      fFirstDirectionReversalProcess(""),
      fFirstDirectionReversalStepStatus(-1),
      fFirstDirectionReversalStepLenNm(-1.),
      fFirstDirectionReversalPreEnergy(-1.),
      fFirstDirectionReversalPostEnergy(-1.),
      fFirstDirectionReversalDeltaThetaDeg(-1.),
      fHasFirstBoundaryCrossing(false),
      fFirstBoundaryStep(-1),
      fFirstBoundaryDepthNm(-1.),
      fFirstBoundaryEnergy(-1.),
      fFirstBoundaryType(0),
      fNStepStatusGeomBoundary(0),
      fNStepStatusPostStepProc(0),
      fNStepStatusAlongStepProc(0),
      fNStepStatusUserLimit(0),
      fNStepStatusOther(0),
      fNProcMsc(0),
      fNProcStepLimiter(0),
      fNProcTransportation(0),
      fNProcEIoni(0),
      fNProcOther(0)
{
}

EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction(const G4Event*)
{
    // Reset per-event accumulator
    fEdepPrimary = 0.;
    fEdepPrimaryWeighted = 0.;
    fPrimaryTrackLength = 0.;
    fNMicroscopicEdep = 0;
    fPrimaryResidualEnergy = 0.;
    fPrimaryEndLocation = 0;
    fPrimaryLastLocation = 0;
    fPrimaryLastProcess = "";
    fPrimaryStopStatus = -1;
    fHasPrimaryExitCandidate = false;
    fPrimaryExitClassCandidate = 0;
    fPrimaryExitEnergyCandidate = 0.;
    fPrimaryExitDirectionCandidate = G4ThreeVector();
    fPrimaryEdepByEIoni = 0.;
    fPrimaryEdepByMsc = 0.;
    fPrimaryEdepByOther = 0.;
    fPrimaryFirstStepEdep = 0.;
    fPrimaryMaxStepEdep = 0.;
    fDepthFirstEdepNm = -1.;
    fPrimaryMaxDepthNm = 0.;
    fPrimaryBoundaryCrossings = 0;
    fPrimaryDirectionReversals = 0;
    fLastPrimaryDirectionSign = 0;
    fHasLastPrimaryDirectionSign = false;
    fPrimaryFirstProcessInAl2O3 = "";
    fHasFirstDirectionReversal = false;
    fFirstDirectionReversalStep = -1;
    fFirstDirectionReversalDepthNm = -1.;
    fFirstDirectionReversalEnergy = -1.;
    fFirstDirectionReversalProcess = "";
    fFirstDirectionReversalStepStatus = -1;
    fFirstDirectionReversalStepLenNm = -1.;
    fFirstDirectionReversalPreEnergy = -1.;
    fFirstDirectionReversalPostEnergy = -1.;
    fFirstDirectionReversalDeltaThetaDeg = -1.;
    fHasFirstBoundaryCrossing = false;
    fFirstBoundaryStep = -1;
    fFirstBoundaryDepthNm = -1.;
    fFirstBoundaryEnergy = -1.;
    fFirstBoundaryType = 0;
    fNStepStatusGeomBoundary = 0;
    fNStepStatusPostStepProc = 0;
    fNStepStatusAlongStepProc = 0;
    fNStepStatusUserLimit = 0;
    fNStepStatusOther = 0;
    fNProcMsc = 0;
    fNProcStepLimiter = 0;
    fNProcTransportation = 0;
    fNProcEIoni = 0;
    fNProcOther = 0;
    fPrimaryTrajectorySteps.clear();
    if (fRunAction && fRunAction->IsTrajectoryDiagnostics()) {
        const auto reserveSize = fRunAction->GetTrajectoryMaxStepsPerEvent();
        if (reserveSize > 0) {
            fPrimaryTrajectorySteps.reserve(static_cast<std::size_t>(reserveSize));
        }
    }
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    auto* analysisManager = G4AnalysisManager::Instance();
    if (!analysisManager) return;

    // Histogram ID 0: primary particle energy deposition in active scoring material
    // Convert from internal units (keV) to eV for display
    analysisManager->FillH1(0, fEdepPrimary / eV);
    if (fRunAction && fRunAction->GetEdepPrimaryWeightedId() >= 0) {
        analysisManager->FillH1(fRunAction->GetEdepPrimaryWeightedId(), fEdepPrimaryWeighted / eV);
    }
    if (fRunAction && fRunAction->GetPrimaryTrackLengthId() >= 0) {
        analysisManager->FillH1(fRunAction->GetPrimaryTrackLengthId(), fPrimaryTrackLength / nm);
    }

    // Histogram ID 2: number of microscopic energy-depositing steps per event
    analysisManager->FillH1(2, static_cast<G4double>(fNMicroscopicEdep));

    // 2D: primary edep vs step count per event
    if (fRunAction && fRunAction->GetEdepVsStepsId() >= 0) {
        analysisManager->FillH2(
            fRunAction->GetEdepVsStepsId(),
            fEdepPrimary / eV,
            static_cast<G4double>(fNMicroscopicEdep)
        );
    }

    // Histogram ID 4: primary residual kinetic energy at end of event
    analysisManager->FillH1(4, fPrimaryResidualEnergy / eV);

    // Histogram ID 5: primary end volume (categorical)
    if (fPrimaryEndLocation == 0 && fPrimaryLastLocation != 0) {
        fPrimaryEndLocation = fPrimaryLastLocation;
    }
    analysisManager->FillH1(5, static_cast<G4double>(fPrimaryEndLocation));

    // Fill exactly one primary-exit classification per event (if any exit occurred).
    // Accept exits only when the final location is outside active scoring material.
    const G4int finalLocation = (fPrimaryEndLocation != 0) ? fPrimaryEndLocation : fPrimaryLastLocation;
    G4int eventExitClass = 4; // 1=entrance, 2=opposite, 3=lateral, 4=stop/no-valid-exit
    G4bool entranceExitSpecularAccepted = false;
    if (fRunAction && fHasPrimaryExitCandidate && finalLocation != 1) {
        eventExitClass = fPrimaryExitClassCandidate;
        const G4int classId = fRunAction->GetPrimaryExitClassId();
        if (classId >= 0) {
            analysisManager->FillH1(classId, static_cast<G4double>(fPrimaryExitClassCandidate));
        }
        if (fPrimaryExitEnergyCandidate > 0.) {
            if (fPrimaryExitClassCandidate == 1) {
                entranceExitSpecularAccepted = true;
                if (fRunAction->IsSpecularAcceptanceEnabled()) {
                    const auto incidentDir = fRunAction->GetPrimaryDirection();
                    if (incidentDir.mag2() > 0. && fPrimaryExitDirectionCandidate.mag2() > 0.) {
                        const G4double sign = (incidentDir.z() >= 0.) ? -1.0 : 1.0;
                        const G4ThreeVector normal(0., 0., sign); // outward normal of entrance side
                        const G4ThreeVector specDir =
                            (incidentDir - 2.0 * incidentDir.dot(normal) * normal).unit();
                        const G4ThreeVector exitDir = fPrimaryExitDirectionCandidate.unit();
                        G4double cosAng = specDir.dot(exitDir);
                        if (cosAng > 1.0) cosAng = 1.0;
                        if (cosAng < -1.0) cosAng = -1.0;
                        const G4double angleDeg = std::acos(cosAng) / deg;
                        entranceExitSpecularAccepted =
                            (angleDeg <= fRunAction->GetSpecularAcceptanceHalfAngleDeg());
                    } else {
                        entranceExitSpecularAccepted = false;
                    }
                }
                const G4int id = fRunAction->GetPrimaryExitEnergyEntranceId();
                if (id >= 0) analysisManager->FillH1(id, fPrimaryExitEnergyCandidate / eV);
                const G4int specId = fRunAction->GetPrimaryExitEnergyEntranceSpecularId();
                if (entranceExitSpecularAccepted && specId >= 0) {
                    analysisManager->FillH1(specId, fPrimaryExitEnergyCandidate / eV);
                }
            } else if (fPrimaryExitClassCandidate == 2) {
                const G4int id = fRunAction->GetPrimaryExitEnergyOppositeId();
                if (id >= 0) analysisManager->FillH1(id, fPrimaryExitEnergyCandidate / eV);
            } else if (fPrimaryExitClassCandidate == 3) {
                const G4int id = fRunAction->GetPrimaryExitEnergyLateralId();
                if (id >= 0) analysisManager->FillH1(id, fPrimaryExitEnergyCandidate / eV);
            }
        }
    }

    // Fill class-conditioned EdepPrimary histograms.
    if (fRunAction) {
        G4int classEdepId = -1;
        if (eventExitClass == 1) {
            classEdepId = fRunAction->GetEdepPrimaryExitEntranceId();
        } else if (eventExitClass == 2) {
            classEdepId = fRunAction->GetEdepPrimaryExitOppositeId();
        } else if (eventExitClass == 3) {
            classEdepId = fRunAction->GetEdepPrimaryExitLateralId();
        } else {
            classEdepId = fRunAction->GetEdepPrimaryStopId();
        }
        if (classEdepId >= 0) {
            analysisManager->FillH1(classEdepId, fEdepPrimary / eV);
        }
    }

    // 2D: residual energy vs end volume category
    if (fRunAction && fRunAction->GetResidualVsEndVolumeId() >= 0) {
        analysisManager->FillH2(
            fRunAction->GetResidualVsEndVolumeId(),
            fPrimaryResidualEnergy / eV,
            static_cast<G4double>(fPrimaryEndLocation)
        );
    }

    auto processCategory = [](const G4String& name) -> G4int {
        if (name == "Transportation") return 1;
        if (name == "eIoni") return 2;
        if (name == "msc") return 3;
        if (name == "eBrem") return 4;
        if (name == "CoulombScat") return 5;
        return 6; // Other/unknown
    };
    auto stopStatusCategory = [](G4int status) -> G4int {
        switch (status) {
            case 0: return 1; // fAlive
            case 1: return 2; // fStopButAlive
            case 2: return 3; // fStopAndKill
            case 3: return 4; // fKillTrackAndSecondaries
            case 4: return 5; // fSuspend
            case 5: return 6; // fPostponeToNextEvent
            default: return 0; // Unknown
        }
    };

    // 2D: residual energy vs last process category
    if (fRunAction && fRunAction->GetResidualVsLastProcessId() >= 0) {
        analysisManager->FillH2(
            fRunAction->GetResidualVsLastProcessId(),
            fPrimaryResidualEnergy / eV,
            static_cast<G4double>(processCategory(fPrimaryLastProcess))
        );
    }

    // 2D: residual energy vs stop status category
    if (fRunAction && fRunAction->GetResidualVsStopStatusId() >= 0) {
        analysisManager->FillH2(
            fRunAction->GetResidualVsStopStatusId(),
            fPrimaryResidualEnergy / eV,
            static_cast<G4double>(stopStatusCategory(fPrimaryStopStatus))
        );
    }

    // One-row-per-event diagnostics ntuple for primary electron transport debugging.
    if (fRunAction && fRunAction->GetEventDiagnosticsNtupleId() >= 0) {
        const G4int ntupleId = fRunAction->GetEventDiagnosticsNtupleId();
        analysisManager->FillNtupleIColumn(ntupleId, 0, event ? event->GetEventID() : -1);
        analysisManager->FillNtupleDColumn(ntupleId, 1, fRunAction->GetPrimaryEnergy() / eV);
        analysisManager->FillNtupleDColumn(ntupleId, 2, fEdepPrimary / eV);
        analysisManager->FillNtupleDColumn(ntupleId, 3, fPrimaryResidualEnergy / eV);
        analysisManager->FillNtupleIColumn(ntupleId, 4, eventExitClass);
        analysisManager->FillNtupleDColumn(ntupleId, 5, fPrimaryExitEnergyCandidate / eV);
        analysisManager->FillNtupleIColumn(ntupleId, 6, fPrimaryStopStatus);
        analysisManager->FillNtupleIColumn(ntupleId, 7, finalLocation);
        analysisManager->FillNtupleIColumn(ntupleId, 8, fNMicroscopicEdep);
        analysisManager->FillNtupleDColumn(ntupleId, 9, fPrimaryTrackLength / nm);
        analysisManager->FillNtupleDColumn(ntupleId, 10, fPrimaryMaxDepthNm);
        analysisManager->FillNtupleIColumn(ntupleId, 11, fPrimaryBoundaryCrossings);
        analysisManager->FillNtupleIColumn(ntupleId, 12, fPrimaryDirectionReversals);
        analysisManager->FillNtupleSColumn(ntupleId, 13, fPrimaryFirstProcessInAl2O3);
        analysisManager->FillNtupleSColumn(ntupleId, 14, fPrimaryLastProcess);
        analysisManager->FillNtupleDColumn(ntupleId, 15, fPrimaryEdepByEIoni / eV);
        analysisManager->FillNtupleDColumn(ntupleId, 16, fPrimaryEdepByMsc / eV);
        analysisManager->FillNtupleDColumn(ntupleId, 17, fPrimaryEdepByOther / eV);
        analysisManager->FillNtupleDColumn(ntupleId, 18, fPrimaryFirstStepEdep / eV);
        analysisManager->FillNtupleDColumn(ntupleId, 19, fPrimaryMaxStepEdep / eV);
        analysisManager->FillNtupleDColumn(ntupleId, 20,
                                           (fDepthFirstEdepNm >= 0.) ? fDepthFirstEdepNm : -1.0);
        analysisManager->FillNtupleIColumn(ntupleId, 21, fFirstDirectionReversalStep);
        analysisManager->FillNtupleDColumn(ntupleId, 22, fFirstDirectionReversalDepthNm);
        analysisManager->FillNtupleDColumn(ntupleId, 23,
                                           (fFirstDirectionReversalEnergy >= 0.)
                                               ? (fFirstDirectionReversalEnergy / eV)
                                               : -1.0);
        analysisManager->FillNtupleIColumn(ntupleId, 24, fFirstBoundaryStep);
        analysisManager->FillNtupleDColumn(ntupleId, 25, fFirstBoundaryDepthNm);
        analysisManager->FillNtupleDColumn(ntupleId, 26,
                                           (fFirstBoundaryEnergy >= 0.)
                                               ? (fFirstBoundaryEnergy / eV)
                                               : -1.0);
        analysisManager->FillNtupleIColumn(ntupleId, 27, fFirstBoundaryType);
        analysisManager->FillNtupleSColumn(ntupleId, 28, fFirstDirectionReversalProcess);
        analysisManager->FillNtupleIColumn(ntupleId, 29, fFirstDirectionReversalStepStatus);
        analysisManager->FillNtupleDColumn(ntupleId, 30, fFirstDirectionReversalStepLenNm);
        analysisManager->FillNtupleDColumn(ntupleId, 31,
                                           (fFirstDirectionReversalPreEnergy >= 0.)
                                               ? (fFirstDirectionReversalPreEnergy / eV)
                                               : -1.0);
        analysisManager->FillNtupleDColumn(ntupleId, 32,
                                           (fFirstDirectionReversalPostEnergy >= 0.)
                                               ? (fFirstDirectionReversalPostEnergy / eV)
                                               : -1.0);
        analysisManager->FillNtupleDColumn(ntupleId, 33, fFirstDirectionReversalDeltaThetaDeg);
        analysisManager->FillNtupleIColumn(ntupleId, 34, fNStepStatusGeomBoundary);
        analysisManager->FillNtupleIColumn(ntupleId, 35, fNStepStatusPostStepProc);
        analysisManager->FillNtupleIColumn(ntupleId, 36, fNStepStatusAlongStepProc);
        analysisManager->FillNtupleIColumn(ntupleId, 37, fNStepStatusUserLimit);
        analysisManager->FillNtupleIColumn(ntupleId, 38, fNStepStatusOther);
        analysisManager->FillNtupleIColumn(ntupleId, 39, fNProcMsc);
        analysisManager->FillNtupleIColumn(ntupleId, 40, fNProcStepLimiter);
        analysisManager->FillNtupleIColumn(ntupleId, 41, fNProcTransportation);
        analysisManager->FillNtupleIColumn(ntupleId, 42, fNProcEIoni);
        analysisManager->FillNtupleIColumn(ntupleId, 43, fNProcOther);
        analysisManager->AddNtupleRow(ntupleId);
    }

    if (fRunAction && fRunAction->IsTrajectoryDiagnostics() &&
        fRunAction->GetTrajectoryDiagnosticsNtupleId() >= 0 &&
        !fPrimaryTrajectorySteps.empty()) {
        G4int sampleIndex = -1;
        if (fRunAction->AcquireTrajectorySampleSlot(eventExitClass, sampleIndex)) {
            const G4int ntupleId = fRunAction->GetTrajectoryDiagnosticsNtupleId();
            const G4int eventId = event ? event->GetEventID() : -1;
            const G4String activeMaterial = fRunAction->GetActiveScoringMaterial();
            for (const auto& step : fPrimaryTrajectorySteps) {
                const G4int isBoundary = (step.preVolume != step.postVolume) ? 1 : 0;
                const G4int isOutwardBoundary =
                    (step.preVolume == activeMaterial && step.postVolume != activeMaterial) ? 1 : 0;
                const G4int isFirstReversalStep =
                    (fFirstDirectionReversalStep > 0 &&
                     step.stepNumber == fFirstDirectionReversalStep)
                        ? 1
                        : 0;
                analysisManager->FillNtupleIColumn(ntupleId, 0, eventId);
                analysisManager->FillNtupleIColumn(ntupleId, 1, sampleIndex);
                analysisManager->FillNtupleIColumn(ntupleId, 2, eventExitClass);
                analysisManager->FillNtupleIColumn(ntupleId, 3, step.stepNumber);
                analysisManager->FillNtupleDColumn(ntupleId, 4, step.preDepthNm);
                analysisManager->FillNtupleDColumn(ntupleId, 5, step.postDepthNm);
                analysisManager->FillNtupleDColumn(ntupleId, 6, step.stepLenNm);
                analysisManager->FillNtupleDColumn(ntupleId, 7, step.preEnergy / eV);
                analysisManager->FillNtupleDColumn(ntupleId, 8, step.postEnergy / eV);
                analysisManager->FillNtupleDColumn(ntupleId, 9, step.edep / eV);
                analysisManager->FillNtupleDColumn(ntupleId, 10, step.dirZPre);
                analysisManager->FillNtupleDColumn(ntupleId, 11, step.dirZPost);
                analysisManager->FillNtupleDColumn(ntupleId, 12, step.deltaThetaDeg);
                analysisManager->FillNtupleIColumn(ntupleId, 13, step.reversalOnStep);
                analysisManager->FillNtupleSColumn(ntupleId, 14, step.processName);
                analysisManager->FillNtupleIColumn(ntupleId, 15, step.stepStatus);
                analysisManager->FillNtupleSColumn(ntupleId, 16, step.preVolume);
                analysisManager->FillNtupleSColumn(ntupleId, 17, step.postVolume);
                analysisManager->FillNtupleIColumn(ntupleId, 18, isBoundary);
                analysisManager->FillNtupleIColumn(ntupleId, 19, isOutwardBoundary);
                analysisManager->FillNtupleIColumn(ntupleId, 20, isFirstReversalStep);
                analysisManager->AddNtupleRow(ntupleId);
            }
        }
    }

    if (fRunAction && fEdepPrimary > 0.) {
        fRunAction->UpdateMinNonZeroEdep(fEdepPrimary);
    }
}

void EventAction::AddPrimaryEdep(G4double edep)
{
    fEdepPrimary += edep;
}

void EventAction::AddPrimaryEdepWeighted(G4double edepWeighted)
{
    fEdepPrimaryWeighted += edepWeighted;
}

void EventAction::AddPrimaryTrackLength(G4double stepLength)
{
    if (stepLength > 0.) {
        fPrimaryTrackLength += stepLength;
    }
}

void EventAction::AddEdepInteraction()
{
    ++fNMicroscopicEdep;
}

void EventAction::UpdatePrimaryResidualEnergy(G4double energy)
{
    if (energy >= 0.) {
        fPrimaryResidualEnergy = energy;
    }
}

void EventAction::UpdatePrimaryEndVolume(const G4String& volumeName)
{
    if (fPrimaryEndLocation != 0) {
        return;
    }
    const G4String activeMaterial = fRunAction ? fRunAction->GetActiveScoringMaterial() : "Al2O3";
    if (volumeName == activeMaterial) {
        fPrimaryEndLocation = 1;
    } else if (volumeName == "World") {
        fPrimaryEndLocation = 2;
    } else if (volumeName == "OutOfWorld") {
        fPrimaryEndLocation = 3;
    } else {
        fPrimaryEndLocation = 4;
    }
}

void EventAction::UpdatePrimaryLastVolume(const G4String& volumeName)
{
    const G4String activeMaterial = fRunAction ? fRunAction->GetActiveScoringMaterial() : "Al2O3";
    if (volumeName == activeMaterial) {
        fPrimaryLastLocation = 1;
    } else if (volumeName == "World") {
        fPrimaryLastLocation = 2;
    } else if (volumeName == "OutOfWorld") {
        fPrimaryLastLocation = 3;
    } else {
        fPrimaryLastLocation = 4;
    }
}

void EventAction::UpdatePrimaryLastProcess(const G4String& processName)
{
    fPrimaryLastProcess = processName;
}

void EventAction::UpdatePrimaryStopStatus(G4int status)
{
    fPrimaryStopStatus = status;
}

void EventAction::UpdatePrimaryExitCandidate(G4int exitClass, G4double kineticEnergy,
                                             const G4ThreeVector& exitDirection)
{
    if (exitClass < 1 || exitClass > 3) {
        return;
    }
    fHasPrimaryExitCandidate = true;
    fPrimaryExitClassCandidate = exitClass;
    fPrimaryExitEnergyCandidate = (kineticEnergy > 0.) ? kineticEnergy : 0.;
    fPrimaryExitDirectionCandidate = exitDirection;
}

void EventAction::AddPrimaryEdepByProcess(const G4String& processName, G4double edep, G4double depthNm)
{
    if (edep <= 0.) {
        return;
    }
    if (fPrimaryFirstStepEdep <= 0.) {
        fPrimaryFirstStepEdep = edep;
        fDepthFirstEdepNm = depthNm;
    }
    if (edep > fPrimaryMaxStepEdep) {
        fPrimaryMaxStepEdep = edep;
    }
    if (processName == "eIoni") {
        fPrimaryEdepByEIoni += edep;
    } else if (processName == "msc") {
        fPrimaryEdepByMsc += edep;
    } else {
        fPrimaryEdepByOther += edep;
    }
}

void EventAction::UpdatePrimaryMaxDepthNm(G4double depthNm)
{
    if (depthNm > fPrimaryMaxDepthNm) {
        fPrimaryMaxDepthNm = depthNm;
    }
}

void EventAction::AddPrimaryStepAudit(const G4String& processName, G4int stepStatus)
{
    switch (stepStatus) {
        case fGeomBoundary:
            ++fNStepStatusGeomBoundary;
            break;
        case fPostStepDoItProc:
            ++fNStepStatusPostStepProc;
            break;
        case fAlongStepDoItProc:
            ++fNStepStatusAlongStepProc;
            break;
        case fUserDefinedLimit:
            ++fNStepStatusUserLimit;
            break;
        default:
            ++fNStepStatusOther;
            break;
    }

    if (processName == "msc") {
        ++fNProcMsc;
    } else if (processName == "StepLimiter") {
        ++fNProcStepLimiter;
    } else if (processName == "Transportation") {
        ++fNProcTransportation;
    } else if (processName == "eIoni") {
        ++fNProcEIoni;
    } else {
        ++fNProcOther;
    }
}

void EventAction::RecordPrimaryBoundaryCrossing(G4int stepNumber, G4double depthNm,
                                                 G4double kineticEnergy, const G4String& preVolume,
                                                 const G4String& postVolume)
{
    ++fPrimaryBoundaryCrossings;
    const G4String activeMaterial = fRunAction ? fRunAction->GetActiveScoringMaterial() : "Al2O3";
    const G4bool isOutwardFromFilm =
        (preVolume == activeMaterial && postVolume != activeMaterial);
    if (isOutwardFromFilm && !fHasFirstBoundaryCrossing) {
        auto boundaryType = [&activeMaterial](const G4String& preName,
                                              const G4String& postName) -> G4int {
            if (preName == activeMaterial && postName == "World") return 1; // out of layer to world
            if (preName == "World" && postName == activeMaterial) return 2; // back into layer
            if (preName == activeMaterial && postName != activeMaterial) return 3; // out to other
            if (preName != activeMaterial && postName == activeMaterial) return 4; // into layer from other
            return 0;
        };
        fHasFirstBoundaryCrossing = true;
        fFirstBoundaryStep = stepNumber;
        fFirstBoundaryDepthNm = depthNm;
        fFirstBoundaryEnergy = kineticEnergy;
        fFirstBoundaryType = boundaryType(preVolume, postVolume);
    }
}

void EventAction::UpdatePrimaryDirectionSignZ(G4double dirZ, G4int stepNumber, G4double depthNm,
                                               G4double kineticEnergy, const G4String& processName,
                                               G4int stepStatus, G4double stepLenNm,
                                               G4double preEnergy, G4double postEnergy,
                                               G4double deltaThetaDeg)
{
    const G4double threshold = 1e-9;
    G4int sign = 0;
    if (dirZ > threshold) {
        sign = 1;
    } else if (dirZ < -threshold) {
        sign = -1;
    }
    if (sign == 0) {
        return;
    }
    if (fHasLastPrimaryDirectionSign && fLastPrimaryDirectionSign != sign) {
        ++fPrimaryDirectionReversals;
        if (!fHasFirstDirectionReversal) {
            fHasFirstDirectionReversal = true;
            fFirstDirectionReversalStep = stepNumber;
            fFirstDirectionReversalDepthNm = depthNm;
            fFirstDirectionReversalEnergy = kineticEnergy;
            fFirstDirectionReversalProcess = processName;
            fFirstDirectionReversalStepStatus = stepStatus;
            fFirstDirectionReversalStepLenNm = stepLenNm;
            fFirstDirectionReversalPreEnergy = preEnergy;
            fFirstDirectionReversalPostEnergy = postEnergy;
            fFirstDirectionReversalDeltaThetaDeg = deltaThetaDeg;
        }
    }
    fLastPrimaryDirectionSign = sign;
    fHasLastPrimaryDirectionSign = true;
}

void EventAction::UpdatePrimaryFirstProcessInAl2O3(const G4String& processName)
{
    if (fPrimaryFirstProcessInAl2O3.empty() && !processName.empty() && processName != "None") {
        fPrimaryFirstProcessInAl2O3 = processName;
    }
}

void EventAction::RecordPrimaryTrajectoryStep(G4int stepNumber, G4double preDepthNm,
                                              G4double postDepthNm, G4double stepLenNm,
                                              G4double preEnergy, G4double postEnergy,
                                              G4double edep, G4double dirZPre, G4double dirZPost,
                                              G4double deltaThetaDeg, G4int reversalOnStep,
                                              const G4String& processName, G4int stepStatus,
                                              const G4String& preVolume,
                                              const G4String& postVolume)
{
    if (!fRunAction || !fRunAction->IsTrajectoryDiagnostics()) {
        return;
    }
    const auto maxSteps = fRunAction->GetTrajectoryMaxStepsPerEvent();
    if (maxSteps > 0 &&
        static_cast<G4int>(fPrimaryTrajectorySteps.size()) >= maxSteps) {
        return;
    }
    PrimaryTrajectoryStep rec;
    rec.stepNumber = stepNumber;
    rec.preDepthNm = preDepthNm;
    rec.postDepthNm = postDepthNm;
    rec.stepLenNm = stepLenNm;
    rec.preEnergy = preEnergy;
    rec.postEnergy = postEnergy;
    rec.edep = edep;
    rec.dirZPre = dirZPre;
    rec.dirZPost = dirZPost;
    rec.deltaThetaDeg = deltaThetaDeg;
    rec.reversalOnStep = reversalOnStep;
    rec.processName = processName;
    rec.stepStatus = stepStatus;
    rec.preVolume = preVolume;
    rec.postVolume = postVolume;
    fPrimaryTrajectorySteps.push_back(rec);
}
