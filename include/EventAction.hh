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

    // Count microscopic energy deposition interactions in Al2O3 per event
    void AddEdepInteraction();

private:
    RunAction* fRunAction;
    G4double fEdepPrimary; // total primary e- energy deposited in Al2O3 for this event
    G4int fNMicroscopicEdep; // number of energy-depositing steps in Al2O3 for this event
};

#endif
