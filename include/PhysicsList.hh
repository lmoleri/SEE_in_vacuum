#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class PhysicsList : public G4VModularPhysicsList
{
public:
    explicit PhysicsList(const G4String& emModel = "PAI");
    virtual ~PhysicsList();

    void SetPaiEnabledOverride(G4bool enabled);
    void SetLivermoreAtomicDeexcitationOverride(G4bool enabled);

    // Configure physics processes
    virtual void ConstructProcess() override;

    // Set production cuts
    virtual void SetCuts() override;

private:
    void ConfigureEmPhysics();
    G4String fEmModel;
    G4int fPaiEnabledOverride;
    G4int fLivermoreAtomicDeexcitationOverride;
};

#endif
