#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"

ActionInitialization::ActionInitialization()
{
}

ActionInitialization::~ActionInitialization()
{
}

void ActionInitialization::BuildForMaster() const
{
}

void ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction());
}
