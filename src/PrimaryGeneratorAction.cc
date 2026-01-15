#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction* runAction)
    : fRunAction(runAction)
{
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);

    // Default particle type: electron
    SetParticleName("e-");

    // Default energy: 1 MeV
    fParticleGun->SetParticleEnergy(1.0 * MeV);

    // Default position: 1 micron above the Al2O3 layer
    // Al2O3 layer is at z=0, thickness is 20nm, so center is at z=0
    // Place gun at z = -1 micron (above the layer, shooting downward)
    fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, -1.0 * um));

    // Default direction: along +z axis (downward toward the Al2O3 layer)
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if (fRunAction) {
        fRunAction->SetPrimaryEnergy(fParticleGun->GetParticleEnergy());
        if (auto* def = fParticleGun->GetParticleDefinition()) {
            fRunAction->SetPrimaryParticleName(def->GetParticleName());
        }
    }
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::SetParticleName(const G4String& name)
{
    auto* particleTable = G4ParticleTable::GetParticleTable();
    auto* particle = particleTable->FindParticle(name);
    if (!particle) {
        G4cerr << "PrimaryGeneratorAction: Unknown particle '" << name
               << "', falling back to e-" << G4endl;
        particle = particleTable->FindParticle("e-");
    }
    fParticleGun->SetParticleDefinition(particle);
    if (fRunAction && particle) {
        fRunAction->SetPrimaryParticleName(particle->GetParticleName());
    }
}
