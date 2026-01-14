#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
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

    // Run action (collect SEY statistics)
    auto* runAction = new RunAction();
    SetUserAction(runAction);

    // Stepping action (detect primaries and secondaries)
    SetUserAction(new SteppingAction(runAction));
}
