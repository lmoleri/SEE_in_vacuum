#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class RunAction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    explicit PrimaryGeneratorAction(RunAction* runAction);
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    void SetEnergy(G4double energy) { fParticleGun->SetParticleEnergy(energy); }
    void SetPosition(G4ThreeVector pos) { fParticleGun->SetParticlePosition(pos); }
    void SetDirection(G4ThreeVector dir) { fParticleGun->SetParticleMomentumDirection(dir); }

private:
    G4ParticleGun* fParticleGun;
    RunAction* fRunAction;
};

#endif
