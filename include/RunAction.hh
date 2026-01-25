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
    void SetMaxPrimaryEnergy(G4double energy);
    void SetPrimaryParticleName(const G4String& name);
    void SetEmModel(const G4String& model);
    void SetSampleThickness(G4double thickness);
    void SetOutputTag(const G4String& tag);
    void SetPaiEnabled(G4bool enabled);
    G4bool IsPaiEnabled() const;

private:
    G4int fNPrimaryElectrons;
    G4int fNSecondaryElectrons;
    G4double fMinNonZeroEdep;
    G4double fPrimaryEnergy;
    G4double fMaxPrimaryEnergy;
    G4String fPrimaryParticleName;
    G4String fEmModel;
    G4double fSampleThickness;
    G4String fOutputTag;
    G4bool fPaiEnabled;
};

#endif

