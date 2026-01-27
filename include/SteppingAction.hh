#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class RunAction;
class EventAction;
class G4Step;

// SteppingAction:
//  - counts electrons leaving the Al2O3 volume into vacuum (world)
//  - accumulates primary particle energy deposition in Al2O3 (electrons, muons, etc.)
class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(RunAction* runAction, EventAction* eventAction);
    virtual ~SteppingAction();

    virtual void UserSteppingAction(const G4Step*) override;

private:
    RunAction*   fRunAction;   // not owned
    EventAction* fEventAction; // not owned
};

#endif
