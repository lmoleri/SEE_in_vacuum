#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class RunAction;
class G4Step;

// SteppingAction:
//  - counts primary electrons
//  - counts electrons leaving the Al2O3 volume into vacuum (world)
class SteppingAction : public G4UserSteppingAction
{
public:
    explicit SteppingAction(RunAction* runAction);
    virtual ~SteppingAction();

    virtual void UserSteppingAction(const G4Step*);

private:
    RunAction* fRunAction;  // not owned
};

#endif

