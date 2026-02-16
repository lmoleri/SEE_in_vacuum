#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4ThreeVector.hh"
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
    void SetPrimaryDirection(const G4ThreeVector& dir);
    G4ThreeVector GetPrimaryDirection() const;
    void SetPrimaryDirectionZ(G4double dirZ);
    void SetSpecularAcceptance(G4bool enabled, G4double halfAngleDeg);
    void SetActiveScoringMaterial(const G4String& material);
    G4bool IsSpecularAcceptanceEnabled() const;
    G4double GetSpecularAcceptanceHalfAngleDeg() const;
    G4String GetActiveScoringMaterial() const;
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
    G4int GetEdepStepDepthPrimaryId() const;
    G4int GetPrimaryTrackLengthDepthId() const;
    G4int GetPrimaryTrackLengthId() const;
    G4int GetPrimaryExitClassId() const;
    G4int GetPrimaryExitEnergyEntranceId() const;
    G4int GetPrimaryExitEnergyEntranceSpecularId() const;
    G4int GetPrimaryExitEnergyOppositeId() const;
    G4int GetPrimaryExitEnergyLateralId() const;
    G4int GetStepLengthAl2O3Id() const;
    G4int GetEdepVsStepsId() const;
    G4int GetResidualVsEndVolumeId() const;
    G4int GetResidualVsLastProcessId() const;
    G4int GetResidualVsStopStatusId() const;
    G4int GetEdepPrimaryStopId() const;
    G4int GetEdepPrimaryExitEntranceId() const;
    G4int GetEdepPrimaryExitOppositeId() const;
    G4int GetEdepPrimaryExitLateralId() const;
    G4int GetEventDiagnosticsNtupleId() const;
    G4double GetPrimaryEnergy() const;
    G4int GetVerboseStepNtupleId() const;
    G4bool IsVerboseStepDiagnostics() const;
    G4double GetVerboseStepThresholdFrac() const;
    G4int GetVerboseStepMaxCount() const;
    G4bool IsTrajectoryDiagnostics() const;
    G4int GetTrajectoryDiagnosticsNtupleId() const;
    G4int GetTrajectorySamplePerClass() const;
    G4int GetTrajectoryMaxStepsPerEvent() const;
    void SetVerboseStepDiagnostics(G4bool enabled);
    void SetVerboseStepThresholdFrac(G4double frac);
    void SetVerboseStepMaxCount(G4int maxCount);
    void SetTrajectoryDiagnostics(G4bool enabled);
    void SetTrajectorySamplePerClass(G4int maxCount);
    void SetTrajectoryMaxStepsPerEvent(G4int maxCount);
    G4int ConsumeVerboseStepSlot();
    G4bool AcquireTrajectorySampleSlot(G4int exitClass, G4int& sampleIndex);

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
    G4ThreeVector fPrimaryDirection;
    G4double fPrimaryDirectionZ;
    G4bool fSpecularAcceptanceEnabled;
    G4double fSpecularAcceptanceHalfAngleDeg;
    G4String fActiveScoringMaterial;
    G4int fEdepPrimaryWeightedId;
    G4int fEdepDepthPrimaryId;
    G4int fEdepDepthPrimaryWeightedId;
    G4int fEdepDepthPrimaryCountsId;
    G4int fEdepStepDepthPrimaryId;
    G4int fPrimaryTrackLengthDepthId;
    G4int fPrimaryTrackLengthId;
    G4int fPrimaryExitClassId;
    G4int fPrimaryExitEnergyEntranceId;
    G4int fPrimaryExitEnergyEntranceSpecularId;
    G4int fPrimaryExitEnergyOppositeId;
    G4int fPrimaryExitEnergyLateralId;
    G4int fStepLengthAl2O3Id;
    G4int fEdepVsStepsId;
    G4int fResidualVsEndVolumeId;
    G4int fResidualVsLastProcessId;
    G4int fResidualVsStopStatusId;
    G4int fEdepPrimaryStopId;
    G4int fEdepPrimaryExitEntranceId;
    G4int fEdepPrimaryExitOppositeId;
    G4int fEdepPrimaryExitLateralId;
    G4int fEventDiagnosticsNtupleId;
    G4int fVerboseStepNtupleId;
    G4String fOutputTag;
    G4bool fPaiEnabled;
    G4int fLivermoreAtomicDeexcitation;
    G4bool fVerboseStepDiagnostics = false;
    G4double fVerboseStepThresholdFrac = 0.9;
    G4int fVerboseStepMaxCount = 1000;
    G4int fVerboseStepUsed = 0;
    G4bool fTrajectoryDiagnostics = false;
    G4int fTrajectorySamplePerClass = 300;
    G4int fTrajectoryMaxStepsPerEvent = 3000;
    G4int fTrajectoryDiagnosticsNtupleId = -1;
    G4int fTrajectoryClass2Used = 0;
    G4int fTrajectoryClass4Used = 0;
};

#endif
