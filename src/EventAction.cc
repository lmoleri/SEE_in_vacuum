#include "EventAction.hh"

#include "RunAction.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"

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
      fPrimaryFirstProcessInAl2O3("")
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
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    auto* analysisManager = G4AnalysisManager::Instance();
    if (!analysisManager) return;

    // Histogram ID 0: primary particle energy deposition in Al2O3
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
    // Accept exits only when the final location is outside Al2O3.
    const G4int finalLocation = (fPrimaryEndLocation != 0) ? fPrimaryEndLocation : fPrimaryLastLocation;
    G4int eventExitClass = 4; // 1=entrance, 2=opposite, 3=lateral, 4=stop/no-valid-exit
    if (fRunAction && fHasPrimaryExitCandidate && finalLocation != 1) {
        eventExitClass = fPrimaryExitClassCandidate;
        const G4int classId = fRunAction->GetPrimaryExitClassId();
        if (classId >= 0) {
            analysisManager->FillH1(classId, static_cast<G4double>(fPrimaryExitClassCandidate));
        }
        if (fPrimaryExitEnergyCandidate > 0.) {
            if (fPrimaryExitClassCandidate == 1) {
                const G4int id = fRunAction->GetPrimaryExitEnergyEntranceId();
                if (id >= 0) analysisManager->FillH1(id, fPrimaryExitEnergyCandidate / eV);
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
        analysisManager->AddNtupleRow(ntupleId);
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
    if (volumeName == "Al2O3") {
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
    if (volumeName == "Al2O3") {
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

void EventAction::UpdatePrimaryExitCandidate(G4int exitClass, G4double kineticEnergy)
{
    if (exitClass < 1 || exitClass > 3) {
        return;
    }
    fHasPrimaryExitCandidate = true;
    fPrimaryExitClassCandidate = exitClass;
    fPrimaryExitEnergyCandidate = (kineticEnergy > 0.) ? kineticEnergy : 0.;
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

void EventAction::AddPrimaryBoundaryCrossing()
{
    ++fPrimaryBoundaryCrossings;
}

void EventAction::UpdatePrimaryDirectionSignZ(G4double dirZ)
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
