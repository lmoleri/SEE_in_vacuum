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
      fPrimaryExitEnergyCandidate(0.)
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
    if (fRunAction && fHasPrimaryExitCandidate && finalLocation != 1) {
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
