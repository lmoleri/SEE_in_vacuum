#include "PhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4ProductionCutsTable.hh"

#include "G4RegionStore.hh"
#include "G4PAIModel.hh"
#include "G4eIonisation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

PhysicsList::PhysicsList()
{
    // Register physics constructors
    RegisterPhysics(new G4EmStandardPhysics_option4());  // Electromagnetic physics
    RegisterPhysics(new G4DecayPhysics());               // Decay physics
    RegisterPhysics(new G4RadioactiveDecayPhysics());    // Radioactive decay
}

PhysicsList::~PhysicsList()
{
}

void PhysicsList::ConstructProcess()
{
    // Let the modular physics list construct all standard processes first
    G4VModularPhysicsList::ConstructProcess();

    // Attach PAI model for electron ionisation in the Al2O3 region
    G4Region* targetRegion =
        G4RegionStore::GetInstance()->GetRegion("Al2O3Region", false);

    if (!targetRegion) {
        // Region not defined yet; fall back to default models
        return;
    }

    auto particleIterator = GetParticleIterator();
    particleIterator->reset();
    while ((*particleIterator)()) {
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        if (!pmanager) continue;

        const G4String& name = particle->GetParticleName();

        if (name == "e-") {
            // Find existing eIonisation process added by EM physics
            G4eIonisation* eIoni = nullptr;
            G4int nProc = pmanager->GetProcessListLength();
            for (G4int i = 0; i < nProc; ++i) {
                auto* proc = (*pmanager->GetProcessList())[i];
                eIoni = dynamic_cast<G4eIonisation*>(proc);
                if (eIoni) break;
            }

            if (eIoni) {
                auto* paiModel = new G4PAIModel(particle, "PAI_e-_Al2O3");
                // Highest priority (0) for this region
                eIoni->AddEmModel(0, paiModel, paiModel, targetRegion);
            }
        }
    }
}

void PhysicsList::SetCuts()
{
    // Set default cut value for all particle types
    SetCutsWithDefault();
    
    // Set production thresholds
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);
}
