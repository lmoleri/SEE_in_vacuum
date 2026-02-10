#include "PhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmParameters.hh"

#include "G4RegionStore.hh"
#include "G4PAIModel.hh"
#include "G4eIonisation.hh"
#include "G4MuIonisation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include <algorithm>
#include <cstdlib>
#include <string>

PhysicsList::PhysicsList(const G4String& emModel)
    : fEmModel(emModel),
      fPaiEnabledOverride(-1),
      fLivermoreAtomicDeexcitationOverride(-1),
      fAtomicDeexcitationOverride(-1),
      fDeexcitationIgnoreCutOverride(-1)
{
    ConfigureEmPhysics();
    RegisterPhysics(new G4DecayPhysics());               // Decay physics
    RegisterPhysics(new G4RadioactiveDecayPhysics());    // Radioactive decay
}

PhysicsList::~PhysicsList()
{
}

void PhysicsList::ConfigureEmPhysics()
{
    G4String model = fEmModel;
    model.toLower();
    if (model == "pai") {
        RegisterPhysics(new G4EmStandardPhysics_option4());
        fPaiEnabledOverride = 1;
        return;
    }
    if (model == "g4emlivermorephysics" || model == "livermore" || model == "livermorephysics") {
        RegisterPhysics(new G4EmLivermorePhysics());
        fPaiEnabledOverride = 0;
        return;
    }
    if (model == "g4empenelopephysics" || model == "penelope" || model == "penelopephysics") {
        RegisterPhysics(new G4EmPenelopePhysics());
        fPaiEnabledOverride = 0;
        return;
    }

    G4cout << "[PhysicsList] Unknown EM model '" << fEmModel
           << "'. Falling back to PAI (Option4 + PAI)." << G4endl;
    RegisterPhysics(new G4EmStandardPhysics_option4());
    fPaiEnabledOverride = 1;
}

void PhysicsList::SetPaiEnabledOverride(G4bool enabled)
{
    fPaiEnabledOverride = enabled ? 1 : 0;
}

void PhysicsList::SetLivermoreAtomicDeexcitationOverride(G4bool enabled)
{
    fLivermoreAtomicDeexcitationOverride = enabled ? 1 : 0;
}

void PhysicsList::SetAtomicDeexcitationOverride(G4bool enabled)
{
    fAtomicDeexcitationOverride = enabled ? 1 : 0;
}

void PhysicsList::SetDeexcitationIgnoreCutOverride(G4bool enabled)
{
    fDeexcitationIgnoreCutOverride = enabled ? 1 : 0;
}

void PhysicsList::ConstructProcess()
{
    // Enforce low-energy tracking thresholds before building EM tables.
    auto* emParams = G4EmParameters::Instance();
    emParams->SetMinEnergy(0.1 * eV);
    emParams->SetLowestElectronEnergy(0.1 * eV);
    emParams->SetLowestMuHadEnergy(0.1 * eV);

    if (fAtomicDeexcitationOverride >= 0) {
        const G4bool enabled = (fAtomicDeexcitationOverride > 0);
        emParams->SetFluo(enabled);
        emParams->SetAuger(enabled);
        // Respect production cuts when deexcitation is disabled
        if (!enabled) {
            emParams->SetDeexcitationIgnoreCut(false);
        }
        G4cout << "[EM] Atomic deexcitation (Fluo/Auger) "
               << (enabled ? "enabled" : "disabled") << " via config override" << G4endl;
    }
    if (fDeexcitationIgnoreCutOverride >= 0) {
        const G4bool ignoreCuts = (fDeexcitationIgnoreCutOverride > 0);
        emParams->SetDeexcitationIgnoreCut(ignoreCuts);
        G4cout << "[EM] DeexcitationIgnoreCut set to "
               << (ignoreCuts ? "true" : "false") << " via config override" << G4endl;
    }

    if (fAtomicDeexcitationOverride < 0 && fLivermoreAtomicDeexcitationOverride >= 0) {
        G4String model = fEmModel;
        model.toLower();
        if (model == "g4emlivermorephysics" || model == "livermore" || model == "livermorephysics") {
            const G4bool enabled = (fLivermoreAtomicDeexcitationOverride > 0);
            emParams->SetFluo(enabled);
            emParams->SetAuger(enabled);
            // Note: PIXE control may not be available in all GEANT4 versions
            // Fluo and Auger are the main atomic deexcitation processes
            G4cout << "[Livermore] Atomic deexcitation (Fluo/Auger) "
                   << (enabled ? "enabled" : "disabled") << " via config override" << G4endl;
        }
    }

    // Let the modular physics list construct all standard processes first
    G4VModularPhysicsList::ConstructProcess();

    if (fPaiEnabledOverride == 0) {
        G4cout << "[PAI] Disabled via config override" << G4endl;
        return;
    }
    if (fPaiEnabledOverride < 0) {
        const char* disablePaiEnv = std::getenv("SEE_DISABLE_PAI");
        if (disablePaiEnv) {
            std::string flag(disablePaiEnv);
            if (flag == "1" || flag == "true" || flag == "TRUE" || flag == "yes" || flag == "YES") {
                G4cout << "[PAI] Disabled via SEE_DISABLE_PAI=" << flag << G4endl;
                return;
            }
        }
    }

    // Attach PAI model for electron ionisation in the Al2O3 region
    G4Region* targetRegion =
        G4RegionStore::GetInstance()->GetRegion("Al2O3Region", false);

    if (!targetRegion) {
        // Region not defined yet; fall back to default models
        G4cout << "[PAI] Al2O3Region not found; PAI not applied" << G4endl;
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
                if (fPaiEnabledOverride > 0) {
                    G4cout << "[PAI] Enabled via config override for e- in Al2O3Region" << G4endl;
                } else {
                    G4cout << "[PAI] Enabled for e- in Al2O3Region" << G4endl;
                }
            }
        } else if (name == "mu-") {
            // Find existing muIonisation process added by EM physics
            G4MuIonisation* muIoni = nullptr;
            G4int nProc = pmanager->GetProcessListLength();
            for (G4int i = 0; i < nProc; ++i) {
                auto* proc = (*pmanager->GetProcessList())[i];
                muIoni = dynamic_cast<G4MuIonisation*>(proc);
                if (muIoni) break;
            }

            if (muIoni) {
                auto* paiModel = new G4PAIModel(particle, "PAI_mu-_Al2O3");
                // Highest priority (0) for this region
                muIoni->AddEmModel(0, paiModel, paiModel, targetRegion);
                if (fPaiEnabledOverride > 0) {
                    G4cout << "[PAI] Enabled via config override for mu- in Al2O3Region" << G4endl;
                } else {
                    G4cout << "[PAI] Enabled for mu- in Al2O3Region" << G4endl;
                }
            }
        }
    }
}

void PhysicsList::SetCuts()
{
    // Set default cut value for all particle types
    SetCutsWithDefault();
    
    // Set production thresholds
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(0.1*eV, 1*GeV);
}
