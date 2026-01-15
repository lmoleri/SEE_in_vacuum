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
    // Primary generator
    SetUserAction(new PrimaryGeneratorAction());

    // Run action (collect SEY statistics and manage analysis)
    auto* runAction = new RunAction();
    SetUserAction(runAction);

    // Event action (per-event accumulation)
    auto* eventAction = new EventAction();
    SetUserAction(eventAction);

    // Stepping action (detect primaries, secondaries, and edep)
    SetUserAction(new SteppingAction(runAction, eventAction));
}
