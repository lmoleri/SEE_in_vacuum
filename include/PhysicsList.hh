#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class PhysicsList : public G4VModularPhysicsList
{
public:
    PhysicsList();
    virtual ~PhysicsList();

    // Configure physics processes
    virtual void ConstructProcess() override;

    // Set production cuts
    virtual void SetCuts() override;
};

#endif
