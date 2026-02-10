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
    void AddEmittedElectron();

    // Track minimum non-zero primary energy deposition per event
    void UpdateMinNonZeroEdep(G4double edep);

    void SetPrimaryEnergy(G4double energy);
    void SetMaxPrimaryEnergy(G4double energy);
    void SetPrimaryParticleName(const G4String& name);
    void SetEmModel(const G4String& model);
    void SetSampleThickness(G4double thickness);
    void SetSubstrateThickness(G4double thickness);
    void SetSampleRadius(G4double radius);
    void SetMaxStep(G4double maxStep);
    void SetSeyAlphaInvNm(G4double alphaInvNm);
    void SetPrimaryDirectionZ(G4double dirZ);
    void SetOutputTag(const G4String& tag);
    void SetPaiEnabled(G4bool enabled);
    void SetLivermoreAtomicDeexcitation(G4int value);
    G4bool IsPaiEnabled() const;
    G4double GetSampleThickness() const;
    G4double GetSubstrateThickness() const;
    G4double GetSampleRadius() const;
    G4double GetMaxStep() const;
    G4double GetSeyAlphaInvNm() const;
    G4double GetPrimaryDirectionZ() const;
    G4int GetEdepPrimaryWeightedId() const;
    G4int GetEdepDepthPrimaryId() const;
    G4int GetEdepDepthPrimaryWeightedId() const;
    G4int GetEdepDepthPrimaryCountsId() const;
    G4int GetPrimaryTrackLengthDepthId() const;
    G4int GetPrimaryTrackLengthId() const;

private:
    void OptimizeHistogramInFile(const G4String& fileName);
    G4int fNPrimaryElectrons;
    G4int fNSecondaryElectrons;
    G4int fNEmittedElectrons;
    G4double fMinNonZeroEdep;
    G4double fPrimaryEnergy;
    G4double fMaxPrimaryEnergy;
    G4String fPrimaryParticleName;
    G4String fEmModel;
    G4double fSampleThickness;
    G4double fSubstrateThickness;
    G4double fSampleRadius;
    G4double fMaxStep;
    G4double fSeyAlphaInvNm;
    G4double fPrimaryDirectionZ;
    G4int fEdepPrimaryWeightedId;
    G4int fEdepDepthPrimaryId;
    G4int fEdepDepthPrimaryWeightedId;
    G4int fEdepDepthPrimaryCountsId;
    G4int fPrimaryTrackLengthDepthId;
    G4int fPrimaryTrackLengthId;
    G4String fOutputTag;
    G4bool fPaiEnabled;
    G4int fLivermoreAtomicDeexcitation;
};

#endif
