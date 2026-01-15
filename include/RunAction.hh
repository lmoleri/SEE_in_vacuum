#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

// Simple RunAction to accumulate secondary electron yield (SEY)
class RunAction : public G4UserRunAction
{
public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    // Called from SteppingAction
    void AddPrimaryElectron();
    void AddSecondaryElectron();

    // Track minimum non-zero primary energy deposition per event
    void UpdateMinNonZeroEdep(G4double edep);

    void SetPrimaryEnergy(G4double energy);
    void SetSampleThickness(G4double thickness);
    void SetOutputTag(const G4String& tag);

private:
    G4int fNPrimaryElectrons;
    G4int fNSecondaryElectrons;
    G4double fMinNonZeroEdep;
    G4double fPrimaryEnergy;
    G4double fSampleThickness;
    G4String fOutputTag;
};

#endif

