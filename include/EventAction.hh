#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

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
    void UpdatePrimaryExitCandidate(G4int exitClass, G4double kineticEnergy);
    void AddPrimaryEdepByProcess(const G4String& processName, G4double edep, G4double depthNm);
    void UpdatePrimaryMaxDepthNm(G4double depthNm);
    void AddPrimaryBoundaryCrossing();
    void UpdatePrimaryDirectionSignZ(G4double dirZ);
    void UpdatePrimaryFirstProcessInAl2O3(const G4String& processName);

private:
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
};

#endif
