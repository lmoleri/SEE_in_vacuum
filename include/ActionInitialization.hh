#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class RunAction;
class EventAction;
class PrimaryGeneratorAction;

class ActionInitialization : public G4VUserActionInitialization
{
public:
    ActionInitialization();
    virtual ~ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

    RunAction* GetRunAction() const { return fRunAction; }
    EventAction* GetEventAction() const { return fEventAction; }
    PrimaryGeneratorAction* GetPrimaryGenerator() const { return fPrimaryGenerator; }

private:
    mutable RunAction* fRunAction = nullptr;
    mutable EventAction* fEventAction = nullptr;
    mutable PrimaryGeneratorAction* fPrimaryGenerator = nullptr;
};

#endif
