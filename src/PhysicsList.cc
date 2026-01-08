#include "PhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList()
{
    // Register physics constructors
    RegisterPhysics(new G4EmStandardPhysics_option4());  // Electromagnetic physics
    RegisterPhysics(new G4DecayPhysics());                // Decay physics
    RegisterPhysics(new G4RadioactiveDecayPhysics());     // Radioactive decay
}

PhysicsList::~PhysicsList()
{
}

void PhysicsList::SetCuts()
{
    // Set default cut value for all particle types
    SetCutsWithDefault();
    
    // Set production thresholds
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);
}
