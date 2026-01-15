#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization()
{
}

ActionInitialization::~ActionInitialization()
{
}

void ActionInitialization::BuildForMaster() const
{
    // For this simple example we don't need actions on the master thread.
}

void ActionInitialization::Build() const
{
    // Run action (collect SEY statistics and manage analysis)
    fRunAction = new RunAction();
    SetUserAction(fRunAction);

    // Primary generator
    fPrimaryGenerator = new PrimaryGeneratorAction(fRunAction);
    SetUserAction(fPrimaryGenerator);

    // Event action (per-event accumulation)
    fEventAction = new EventAction(fRunAction);
    SetUserAction(fEventAction);

    // Stepping action (detect primaries, secondaries, and edep)
    SetUserAction(new SteppingAction(fRunAction, fEventAction));
}
