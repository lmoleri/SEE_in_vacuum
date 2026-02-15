#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <vector>

class G4Event;
class RunAction;

// Collect per-event observables, e.g. primary e- energy deposition in Al2O3
class EventAction : public G4UserEventAction
{
public:
    explicit EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event) override;
    virtual void EndOfEventAction(const G4Event* event) override;

    // Called from SteppingAction to accumulate primary e- energy deposition
    void AddPrimaryEdep(G4double edep);
    void AddPrimaryEdepWeighted(G4double edepWeighted);
    void AddPrimaryTrackLength(G4double stepLength);

    // Count microscopic energy deposition interactions in Al2O3 per event
    void AddEdepInteraction();
    void UpdatePrimaryResidualEnergy(G4double energy);
    void UpdatePrimaryEndVolume(const G4String& volumeName);
    void UpdatePrimaryLastVolume(const G4String& volumeName);
    void UpdatePrimaryLastProcess(const G4String& processName);
    void UpdatePrimaryStopStatus(G4int status);
    void UpdatePrimaryExitCandidate(G4int exitClass, G4double kineticEnergy,
                                    const G4ThreeVector& exitDirection);
    void AddPrimaryEdepByProcess(const G4String& processName, G4double edep, G4double depthNm);
    void UpdatePrimaryMaxDepthNm(G4double depthNm);
    void AddPrimaryStepAudit(const G4String& processName, G4int stepStatus);
    void RecordPrimaryBoundaryCrossing(G4int stepNumber, G4double depthNm, G4double kineticEnergy,
                                       const G4String& preVolume, const G4String& postVolume);
    void UpdatePrimaryDirectionSignZ(G4double dirZ, G4int stepNumber, G4double depthNm,
                                     G4double kineticEnergy, const G4String& processName,
                                     G4int stepStatus, G4double stepLenNm, G4double preEnergyEv,
                                     G4double postEnergyEv, G4double deltaThetaDeg);
    void UpdatePrimaryFirstProcessInAl2O3(const G4String& processName);
    void RecordPrimaryTrajectoryStep(G4int stepNumber, G4double preDepthNm, G4double postDepthNm,
                                     G4double stepLenNm, G4double preEnergy, G4double postEnergy,
                                     G4double edep, G4double dirZPre, G4double dirZPost,
                                     G4double deltaThetaDeg, G4int reversalOnStep,
                                     const G4String& processName, G4int stepStatus,
                                     const G4String& preVolume, const G4String& postVolume);

private:
    struct PrimaryTrajectoryStep {
        G4int stepNumber = 0;
        G4double preDepthNm = 0.;
        G4double postDepthNm = 0.;
        G4double stepLenNm = 0.;
        G4double preEnergy = 0.;
        G4double postEnergy = 0.;
        G4double edep = 0.;
        G4double dirZPre = 0.;
        G4double dirZPost = 0.;
        G4double deltaThetaDeg = 0.;
        G4int reversalOnStep = 0;
        G4String processName;
        G4int stepStatus = 0;
        G4String preVolume;
        G4String postVolume;
    };

    RunAction* fRunAction;
    G4double fEdepPrimary; // total primary e- energy deposited in Al2O3 for this event
    G4double fEdepPrimaryWeighted; // depth-weighted primary e- energy deposition
    G4double fPrimaryTrackLength; // primary track length in Al2O3 for this event
    G4int fNMicroscopicEdep; // number of energy-depositing steps in Al2O3 for this event
    G4double fPrimaryResidualEnergy;
    G4int fPrimaryEndLocation;
    G4int fPrimaryLastLocation;
    G4String fPrimaryLastProcess;
    G4int fPrimaryStopStatus;
    G4bool fHasPrimaryExitCandidate;
    G4int fPrimaryExitClassCandidate;
    G4double fPrimaryExitEnergyCandidate;
    G4ThreeVector fPrimaryExitDirectionCandidate;
    G4double fPrimaryEdepByEIoni;
    G4double fPrimaryEdepByMsc;
    G4double fPrimaryEdepByOther;
    G4double fPrimaryFirstStepEdep;
    G4double fPrimaryMaxStepEdep;
    G4double fDepthFirstEdepNm;
    G4double fPrimaryMaxDepthNm;
    G4int fPrimaryBoundaryCrossings;
    G4int fPrimaryDirectionReversals;
    G4int fLastPrimaryDirectionSign;
    G4bool fHasLastPrimaryDirectionSign;
    G4String fPrimaryFirstProcessInAl2O3;
    G4bool fHasFirstDirectionReversal;
    G4int fFirstDirectionReversalStep;
    G4double fFirstDirectionReversalDepthNm;
    G4double fFirstDirectionReversalEnergy;
    G4String fFirstDirectionReversalProcess;
    G4int fFirstDirectionReversalStepStatus;
    G4double fFirstDirectionReversalStepLenNm;
    G4double fFirstDirectionReversalPreEnergy;
    G4double fFirstDirectionReversalPostEnergy;
    G4double fFirstDirectionReversalDeltaThetaDeg;
    G4bool fHasFirstBoundaryCrossing;
    G4int fFirstBoundaryStep;
    G4double fFirstBoundaryDepthNm;
    G4double fFirstBoundaryEnergy;
    G4int fFirstBoundaryType;
    G4int fNStepStatusGeomBoundary;
    G4int fNStepStatusPostStepProc;
    G4int fNStepStatusAlongStepProc;
    G4int fNStepStatusUserLimit;
    G4int fNStepStatusOther;
    G4int fNProcMsc;
    G4int fNProcStepLimiter;
    G4int fNProcTransportation;
    G4int fNProcEIoni;
    G4int fNProcOther;
    std::vector<PrimaryTrajectoryStep> fPrimaryTrajectorySteps;
};

#endif
