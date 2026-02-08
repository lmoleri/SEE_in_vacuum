#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4EmParameters.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "DetectorConstruction.hh"

#include <cmath>
#include <iomanip>

// ROOT headers for histogram optimization
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

namespace {
void PrintRegionCuts(const G4String& regionName)
{
    auto* region = G4RegionStore::GetInstance()->GetRegion(regionName, false);
    if (!region) {
        G4cout << "  Region not found: " << regionName << G4endl;
        return;
    }
    auto* cuts = region->GetProductionCuts();
    if (!cuts) {
        G4cout << "  Region " << regionName << ": no production cuts assigned." << G4endl;
        return;
    }
    auto* table = G4ProductionCutsTable::GetProductionCutsTable();
    if (!table) {
        G4cout << "  ProductionCutsTable not available." << G4endl;
        return;
    }

    struct ParticleCutInfo {
        const char* name;
        const G4ParticleDefinition* particle;
        G4int index;
    };
    const ParticleCutInfo particles[] = {
        {"gamma", G4Gamma::Gamma(), idxG4GammaCut},
        {"e-", G4Electron::Electron(), idxG4ElectronCut},
        {"e+", G4Positron::Positron(), idxG4PositronCut},
        {"proton", G4Proton::Proton(), idxG4ProtonCut},
    };

    const auto nMaterials = region->GetNumberOfMaterials();
    if (nMaterials == 0) {
        G4cout << "  Region " << regionName << ": no materials registered." << G4endl;
        return;
    }

    G4cout << "  Region: " << regionName << G4endl;
    auto matBegin = region->GetMaterialIterator();
    auto matEnd = matBegin + nMaterials;
    for (auto it = matBegin; it != matEnd; ++it) {
        const auto* material = *it;
        if (!material) continue;
        G4cout << "    Material: " << material->GetName() << G4endl;
        for (const auto& p : particles) {
            const G4double range = cuts->GetProductionCut(p.index);
            const G4double energy = table->ConvertRangeToEnergy(p.particle, material, range);
            G4cout << "      " << std::setw(6) << p.name
                   << " cut: " << std::setw(10) << range / nm << " nm"
                   << " (" << range / mm << " mm)";
            if (energy > 0.) {
                G4cout << " -> " << energy / eV << " eV";
            } else {
                G4cout << " -> n/a";
            }
            G4cout << G4endl;
        }
    }
}
}  // namespace

RunAction::RunAction()
    : G4UserRunAction(),
      fNPrimaryElectrons(0),
      fNSecondaryElectrons(0),
      fNEmittedElectrons(0),
      fMinNonZeroEdep(-1.),
      fPrimaryEnergy(0.),
      fMaxPrimaryEnergy(0.),
      fPrimaryParticleName("e-"),
      fEmModel("PAI"),
      fSampleThickness(0.),
      fSubstrateThickness(0.),
      fOutputTag("SEE_in_vacuum"),
      fPaiEnabled(false),
      fLivermoreAtomicDeexcitation(-1),
      fSeyAlphaInvNm(0.0),
      fPrimaryDirectionZ(1.0),
      fEdepPrimaryWeightedId(-1),
      fEdepDepthPrimaryId(-1),
      fEdepDepthPrimaryWeightedId(-1),
      fPrimaryTrackLengthId(-1)
{
}

RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    fNPrimaryElectrons   = 0;
    fNSecondaryElectrons = 0;
    fNEmittedElectrons = 0;
    fMinNonZeroEdep = -1.;

    auto* analysisManager = G4AnalysisManager::Instance();

    // Basic configuration
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFirstHistoId(0);

    // Capture sample thickness from detector if not set
    if (fSampleThickness <= 0.) {
        auto* det = dynamic_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        if (det) {
            fSampleThickness = det->GetSampleThickness();
        }
    }
    // Capture substrate thickness from detector if not set
    if (fSubstrateThickness < 0.) {
        fSubstrateThickness = 0.;
    }
    if (fSubstrateThickness == 0.) {
        auto* det = dynamic_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        if (det) {
            fSubstrateThickness = det->GetSubstrateThickness();
        }
    }

    static G4bool metaCreated = false;
    static G4bool histosCreated = false;
    if (!metaCreated) {
        analysisManager->CreateNtuple("RunMeta", "Run metadata");
        analysisManager->CreateNtupleDColumn("primaryEnergyMeV");
        analysisManager->CreateNtupleDColumn("sampleThicknessNm");
        analysisManager->CreateNtupleDColumn("substrateThicknessNm");
        analysisManager->CreateNtupleDColumn("maxPrimaryEnergyMeV");
        analysisManager->CreateNtupleIColumn("paiEnabled");
        analysisManager->CreateNtupleSColumn("primaryParticle");
        analysisManager->CreateNtupleSColumn("emModel");
        analysisManager->CreateNtupleIColumn("livermoreAtomicDeexcitation");
        analysisManager->CreateNtupleIColumn("primaryElectrons");
        analysisManager->CreateNtupleIColumn("secondaryElectrons");
        analysisManager->CreateNtupleIColumn("emittedElectrons");
        analysisManager->CreateNtupleDColumn("sey");
        analysisManager->CreateNtupleDColumn("minNonZeroPrimaryEdepEv");
        analysisManager->FinishNtuple();
        metaCreated = true;
    }

    if (histosCreated) {
        analysisManager->Reset();
    } else {
        // Create a 1D histogram for total energy deposition in Al2O3
        // ID 0: EdepPrimary (includes primary particle + all its descendants/secondaries)
        // Use a reasonable range for energy deposition, not the full primary energy
        // Energy deposition is typically much smaller than primary energy
        // For electrons: typically 0-1000 eV, for muons: typically 0-10000 eV
        // Use adaptive range based on primary particle type and energy
        G4double maxEdepRange = 10000. * eV; // Default: 10 keV max (covers most cases)
        
        // For high-energy particles, use a larger but still reasonable range
        G4double primaryEnergy = fPrimaryEnergy;
        if (primaryEnergy <= 0.) {
            primaryEnergy = fMaxPrimaryEnergy;
        }
        
        // Estimate reasonable max energy deposition based on particle type and energy
        // For electrons: energy deposition scales roughly with sqrt(energy) up to ~1 MeV, then more slowly
        // For muons: energy deposition is typically much smaller
        if (fPrimaryParticleName == "e-" || fPrimaryParticleName == "e+") {
            // Electrons: use min(primaryEnergy/100, 10000 eV) to avoid huge ranges
            maxEdepRange = std::min(primaryEnergy / 100.0, 10000. * eV);
            // But ensure at least 1000 eV range for reasonable statistics
            if (maxEdepRange < 1000. * eV) {
                maxEdepRange = 1000. * eV;
            }
        } else if (fPrimaryParticleName == "mu-" || fPrimaryParticleName == "mu+") {
            // Muons: energy deposition is typically very small, use 10000 eV max
            maxEdepRange = 10000. * eV;
        }
        
        // Use fine binning: 1 eV bins for good resolution
        // This ensures we can properly resolve the actual data distribution
        const G4double binWidth = 1.0 * eV;
        const G4int primaryBins = static_cast<G4int>(std::ceil(maxEdepRange / binWidth));
        // Cap at reasonable maximum to avoid memory issues, but ensure fine binning
        const G4int maxPrimaryBins = 100000;
        const G4int finalBins = std::min(primaryBins, maxPrimaryBins);

        G4int histoId = analysisManager->CreateH1(
            "EdepPrimary",
            "Total energy deposition in Al_{2}O_{3} (primary + descendants)",
            finalBins,   // 1 eV per bin for fine resolution
            0.,    // Edep min (eV)
            maxEdepRange / eV   // Edep max (eV) - reasonable range, not full primary energy
        );

        // Set axis labels explicitly (X-axis in eV units)
        analysisManager->SetH1XAxisTitle(histoId, "Energy deposition (eV)");
        analysisManager->SetH1YAxisTitle(histoId, "Number of events");


        // Get maxEnergy for other histograms (use primary energy for residual energy, etc.)
        G4double maxEnergy = fMaxPrimaryEnergy;
        if (maxEnergy <= 0.) {
            maxEnergy = fPrimaryEnergy;
        }
        if (maxEnergy <= 0.) {
            maxEnergy = 200. * eV;
        }

        // Create a 1D histogram for energy deposited per step in Al2O3
        // ID 2: EdepStep
        const G4double maxStepEdepEv =
            std::max(50.0, std::min(maxEnergy / eV, 50000.0));
        const G4int maxStepBins = 200000;
        const G4double idealStepBins = std::ceil(maxStepEdepEv);
        const G4int stepBins = std::max(
            1, static_cast<G4int>(std::min(idealStepBins, static_cast<G4double>(maxStepBins))));
        G4int stepEdepId = analysisManager->CreateH1(
            "EdepStep",
            "Energy deposition per step in Al_{2}O_{3}",
            stepBins,  // ~1 eV per bin up to capped max
            0.,
            maxStepEdepEv
        );
        analysisManager->SetH1XAxisTitle(stepEdepId, "Energy deposition per step (eV)");
        analysisManager->SetH1YAxisTitle(stepEdepId, "Number of steps");

        // Create a 1D histogram for the number of microscopic energy-depositing
        // steps in Al2O3 per event
        const G4int stepsMax = 200;
        G4int stepsHistoId = analysisManager->CreateH1(
            "EdepInteractions",
            "Microscopic energy-depositing steps in Al_{2}O_{3} per event",
            stepsMax,   // number of bins (1 step per bin)
            0.,   // min count
            stepsMax   // max count
        );

        analysisManager->SetH1XAxisTitle(stepsHistoId, "Energy-depositing steps per event");
        analysisManager->SetH1YAxisTitle(stepsHistoId, "Number of events");

        // Create a 1D histogram for PAI electron-ionisation energy transfers
        // ID 3: PAITransfer
        const G4double maxPaiTransferEv =
            std::max(50.0, std::min(maxEnergy / eV, 50000.0));
        const G4int maxPaiBins = 200000;
        const G4double idealPaiBins = std::ceil(maxPaiTransferEv);
        const G4int paiBins = std::max(
            1, static_cast<G4int>(std::min(idealPaiBins, static_cast<G4double>(maxPaiBins))));
        G4int paiTransferId = analysisManager->CreateH1(
            "PAITransfer",
            "PAI energy transfer per step in Al_{2}O_{3}",
            paiBins,  // ~1 eV per bin up to capped max
            0.,
            maxPaiTransferEv
        );
        analysisManager->SetH1XAxisTitle(paiTransferId, "Energy transfer (eV)");
        analysisManager->SetH1YAxisTitle(paiTransferId, "Number of PAI steps");

        // Create a 1D histogram for primary residual kinetic energy at end of event
        // ID 4: PrimaryResidualEnergy
        // Use primary energy range for residual energy (can be up to full primary energy)
        const G4int maxResidualBins = 200000;
        const G4double idealResidualBins = std::ceil(maxEnergy / eV);
        const G4int residualBins = std::max(
            1, static_cast<G4int>(std::min(idealResidualBins, static_cast<G4double>(maxResidualBins))));
        G4int residualId = analysisManager->CreateH1(
            "PrimaryResidualEnergy",
            "Primary residual kinetic energy at end of event",
            residualBins, // 1 eV per bin across scan max energy
            0.,
            maxEnergy / eV
        );
        analysisManager->SetH1XAxisTitle(residualId, "Residual kinetic energy (eV)");
        analysisManager->SetH1YAxisTitle(residualId, "Number of events");

        // Create a 1D histogram for primary end volume category
        // ID 5: PrimaryEndVolume
        G4int endVolId = analysisManager->CreateH1(
            "PrimaryEndVolume",
            "Primary end volume category",
            5,   // 0: unknown, 1: Al2O3, 2: World, 3: OutOfWorld, 4: Other
            0.,
            5.
        );
        analysisManager->SetH1XAxisTitle(endVolId, "End volume category");
        analysisManager->SetH1YAxisTitle(endVolId, "Number of events");

        // Create a 1D histogram for step length in Al2O3
        // ID 6: StepLengthAl2O3
        // Use a larger range to accommodate longer steps, especially for higher energy electrons
        // Minimum 100nm, maximum 10000nm (10 microns)
        const G4double thicknessNm = (fSampleThickness > 0.) ? fSampleThickness / nm : 0.0;
        const G4double maxStepLenNm = std::max(100.0, std::min(thicknessNm * 100.0, 10000.0));
        const G4double stepLenBinNm = 0.1;
        const G4int maxStepLenBins = 200000;
        const G4double idealStepLenBins = std::ceil(maxStepLenNm / stepLenBinNm);
        const G4int stepLenBins = std::max(
            1, static_cast<G4int>(std::min(idealStepLenBins, static_cast<G4double>(maxStepLenBins))));
        G4int stepLenId = analysisManager->CreateH1(
            "StepLengthAl2O3",
            "Step length in Al_{2}O_{3}",
            stepLenBins,  // ~0.1 nm bins up to capped max
            0.,
            maxStepLenNm
        );
        analysisManager->SetH1XAxisTitle(stepLenId, "Step length (nm)");
        analysisManager->SetH1YAxisTitle(stepLenId, "Number of steps");

        // Create depth profiles for primary electron energy deposition (unweighted + weighted).
        if (thicknessNm > 0.) {
            const G4double depthBinNm = 0.05;
            const G4int maxDepthBins = 4000;
            const G4double idealDepthBins = std::ceil(thicknessNm / depthBinNm);
            const G4int depthBins = std::max(
                1, static_cast<G4int>(std::min(idealDepthBins, static_cast<G4double>(maxDepthBins))));
            fEdepDepthPrimaryId = analysisManager->CreateH1(
                "EdepDepthPrimary",
                "Primary e- energy deposition vs depth in Al_{2}O_{3}",
                depthBins,
                0.,
                thicknessNm
            );
            analysisManager->SetH1XAxisTitle(fEdepDepthPrimaryId, "Depth from entrance (nm)");
            analysisManager->SetH1YAxisTitle(fEdepDepthPrimaryId, "Energy deposition (eV)");

            fEdepDepthPrimaryWeightedId = analysisManager->CreateH1(
                "EdepDepthPrimaryWeighted",
                "Depth-weighted primary e- energy deposition vs depth in Al_{2}O_{3}",
                depthBins,
                0.,
                thicknessNm
            );
            analysisManager->SetH1XAxisTitle(fEdepDepthPrimaryWeightedId, "Depth from entrance (nm)");
            analysisManager->SetH1YAxisTitle(fEdepDepthPrimaryWeightedId, "Weighted energy deposition (eV)");
        }

        // Create a 1D histogram for depth-weighted energy deposition in Al2O3 (primary electron only).
        // Append at the end to preserve existing H1 IDs used elsewhere.
        fEdepPrimaryWeightedId = analysisManager->CreateH1(
            "EdepPrimaryWeighted",
            "Depth-weighted energy deposition in Al_{2}O_{3} (primary e-)",
            finalBins,
            0.,
            maxEdepRange / eV
        );
        analysisManager->SetH1XAxisTitle(fEdepPrimaryWeightedId, "Weighted energy deposition (eV)");
        analysisManager->SetH1YAxisTitle(fEdepPrimaryWeightedId, "Number of events");

        // Create a 1D histogram for primary track length in Al2O3 (per event).
        // Append at the end to keep existing H1 IDs stable.
        const G4double maxTrackLenNm = std::max(100.0, std::min(thicknessNm * 200.0, 10000.0));
        const G4double trackLenBinNm = 0.1;
        const G4int maxTrackLenBins = 200000;
        const G4double idealTrackLenBins = std::ceil(maxTrackLenNm / trackLenBinNm);
        const G4int trackLenBins = std::max(
            1, static_cast<G4int>(std::min(idealTrackLenBins, static_cast<G4double>(maxTrackLenBins))));
        fPrimaryTrackLengthId = analysisManager->CreateH1(
            "PrimaryTrackLengthAl2O3",
            "Primary track length in Al_{2}O_{3} per event",
            trackLenBins,
            0.,
            maxTrackLenNm
        );
        analysisManager->SetH1XAxisTitle(fPrimaryTrackLengthId, "Track length (nm)");
        analysisManager->SetH1YAxisTitle(fPrimaryTrackLengthId, "Number of events");

        // Create a 2D histogram: event edep vs number of steps in Al2O3
        // ID 0 for H2: EdepPrimaryVsSteps
        // Use the same binning as EdepPrimary (maxEdepRange, not maxEnergy)
        G4int edepVsStepsId = analysisManager->CreateH2(
            "EdepPrimaryVsSteps",
            "Primary energy deposition vs steps in Al_{2}O_{3}",
            finalBins,  // Use same bins as EdepPrimary
            0.,
            maxEdepRange / eV,  // Use same range as EdepPrimary
            stepsMax,
            0.,
            stepsMax
        );
        analysisManager->SetH2XAxisTitle(edepVsStepsId, "Primary energy deposition (eV)");
        analysisManager->SetH2YAxisTitle(edepVsStepsId, "Energy-depositing steps per event");

        // Create a 2D histogram: primary residual energy vs end volume category
        // ID 1 for H2: ResidualEnergyVsEndVolume
        G4int resVsEndId = analysisManager->CreateH2(
            "ResidualEnergyVsEndVolume",
            "Primary residual energy vs end volume",
            primaryBins,
            0.,
            maxEnergy / eV,
            5,
            0.,
            5.
        );
        analysisManager->SetH2XAxisTitle(resVsEndId, "Residual kinetic energy (eV)");
        analysisManager->SetH2YAxisTitle(resVsEndId, "End volume category");

        // Create a 2D histogram: primary residual energy vs last process category
        // ID 2 for H2: ResidualEnergyVsLastProcess
        G4int resVsProcId = analysisManager->CreateH2(
            "ResidualEnergyVsLastProcess",
            "Primary residual energy vs last process",
            primaryBins,
            0.,
            maxEnergy / eV,
            7,
            0.,
            7.
        );
        analysisManager->SetH2XAxisTitle(resVsProcId, "Residual kinetic energy (eV)");
        analysisManager->SetH2YAxisTitle(resVsProcId, "Last process category");

        // Create a 2D histogram: primary residual energy vs stop status category
        // ID 3 for H2: ResidualEnergyVsStopStatus
        G4int resVsStopId = analysisManager->CreateH2(
            "ResidualEnergyVsStopStatus",
            "Primary residual energy vs stop status",
            primaryBins,
            0.,
            maxEnergy / eV,
            7,
            0.,
            7.
        );
        analysisManager->SetH2XAxisTitle(resVsStopId, "Residual kinetic energy (eV)");
        analysisManager->SetH2YAxisTitle(resVsStopId, "Stop status category");

        histosCreated = true;
    }

    // Dump effective cuts and EM step settings once for diagnostics
    static G4bool printedCuts = false;
    if (!printedCuts) {
        printedCuts = true;
        auto* cutsTable = G4ProductionCutsTable::GetProductionCutsTable();
        if (cutsTable) {
            G4cout << "\n--- Production cuts table ---" << G4endl;
            G4cout << "  Energy range: "
                   << cutsTable->GetLowEdgeEnergy() / eV << " eV -> "
                   << cutsTable->GetHighEdgeEnergy() / eV << " eV" << G4endl;
        }
        G4cout << "\n--- Effective production cuts ---" << G4endl;
        PrintRegionCuts("Al2O3Region");
        PrintRegionCuts("DefaultRegionForTheWorld");

        auto* emParams = G4EmParameters::Instance();
        if (emParams) {
            G4cout << "\n--- EM step settings (global) ---" << G4endl;
            G4cout << "  MscRangeFactor (e-/e+): " << emParams->MscRangeFactor() << G4endl;
            G4cout << "  MscGeomFactor: " << emParams->MscGeomFactor() << G4endl;
            G4cout << "  MscSafetyFactor: " << emParams->MscSafetyFactor() << G4endl;
            G4cout << "  MscLambdaLimit: " << emParams->MscLambdaLimit() / mm << " mm" << G4endl;
            G4cout << "  MscSkin: " << emParams->MscSkin() << G4endl;
            G4cout << "  LinearLossLimit: " << emParams->LinearLossLimit() << G4endl;
            G4cout << "  MinKinEnergy: " << emParams->MinKinEnergy() / eV << " eV" << G4endl;
            G4cout << "  LowestElectronEnergy: " << emParams->LowestElectronEnergy() / eV << " eV" << G4endl;
        }
    }

    // Open output file for this run
    G4String fileName = fOutputTag;
    if (fileName.empty()) {
        fileName = "SEE_in_vacuum";
    }
    if (fileName.size() < 5 || fileName.substr(fileName.size() - 5) != ".root") {
        fileName += ".root";
    }
    analysisManager->OpenFile(fileName);
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    const auto nEvents = run->GetNumberOfEvent();

    // Define "primary electrons" as the number of incident beam particles,
    // i.e. one primary electron per event.
    fNPrimaryElectrons = nEvents;

    G4cout << "\n=== Run summary (SEY) ===" << G4endl;
    G4cout << "  Events processed       : " << nEvents << G4endl;
    G4cout << "  Primary electrons      : " << fNPrimaryElectrons << G4endl;
    G4cout << "  Secondary electrons out: " << fNSecondaryElectrons << G4endl;
    G4cout << "  Emitted electrons out  : " << fNEmittedElectrons << G4endl;

    if (fNPrimaryElectrons > 0) {
        const G4double sey = static_cast<G4double>(fNEmittedElectrons) /
                             static_cast<G4double>(fNPrimaryElectrons);
        G4cout << "  Secondary electron yield (SEY) = "
               << sey << G4endl;
    } else {
        G4cout << "  Secondary electron yield (SEY) not defined (no primaries counted)"
               << G4endl;
    }
    if (fMinNonZeroEdep > 0.) {
        G4cout << "  Min non-zero primary edep : " << fMinNonZeroEdep / eV << " eV"
               << G4endl;
    } else {
        G4cout << "  Min non-zero primary edep : none" << G4endl;
    }
    G4cout << "==========================\n" << G4endl;

    // Write and close analysis output
    auto* analysisManager = G4AnalysisManager::Instance();
    if (analysisManager) {
        analysisManager->FillNtupleDColumn(0, fPrimaryEnergy / MeV);
        analysisManager->FillNtupleDColumn(1, fSampleThickness / nm);
        analysisManager->FillNtupleDColumn(2, fSubstrateThickness / nm);
        analysisManager->FillNtupleDColumn(3, fMaxPrimaryEnergy / MeV);
        analysisManager->FillNtupleIColumn(4, fPaiEnabled ? 1 : 0);
        analysisManager->FillNtupleSColumn(5, fPrimaryParticleName);
        analysisManager->FillNtupleSColumn(6, fEmModel);
        analysisManager->FillNtupleIColumn(7, fLivermoreAtomicDeexcitation);
        analysisManager->FillNtupleIColumn(8, fNPrimaryElectrons);
        analysisManager->FillNtupleIColumn(9, fNSecondaryElectrons);
        analysisManager->FillNtupleIColumn(10, fNEmittedElectrons);
        const G4double sey = (fNPrimaryElectrons > 0)
                                 ? static_cast<G4double>(fNEmittedElectrons) /
                                       static_cast<G4double>(fNPrimaryElectrons)
                                 : 0.0;
        analysisManager->FillNtupleDColumn(11, sey);
        const G4double minNonZeroEdepEv = (fMinNonZeroEdep > 0.) ? (fMinNonZeroEdep / eV) : 0.0;
        analysisManager->FillNtupleDColumn(12, minNonZeroEdepEv);
        analysisManager->AddNtupleRow();

        analysisManager->Write();
        
        // Get the file name to optimize histograms
        G4String fileName = fOutputTag;
        if (fileName.empty()) {
            fileName = "SEE_in_vacuum";
        }
        if (fileName.size() < 5 || fileName.substr(fileName.size() - 5) != ".root") {
            fileName += ".root";
        }
        
        analysisManager->CloseFile();
        
        // Optimize EdepPrimary histogram in the ROOT file
        OptimizeHistogramInFile(fileName);
    }
}

void RunAction::AddPrimaryElectron()
{
    // No-op for now; primary count is taken from the number of events.
}

void RunAction::AddSecondaryElectron()
{
    ++fNSecondaryElectrons;
}

void RunAction::AddEmittedElectron()
{
    ++fNEmittedElectrons;
}

void RunAction::UpdateMinNonZeroEdep(G4double edep)
{
    if (edep <= 0.) {
        return;
    }
    if (fMinNonZeroEdep < 0. || edep < fMinNonZeroEdep) {
        fMinNonZeroEdep = edep;
    }
}

void RunAction::SetPrimaryEnergy(G4double energy)
{
    fPrimaryEnergy = energy;
}

void RunAction::SetMaxPrimaryEnergy(G4double energy)
{
    fMaxPrimaryEnergy = energy;
}

void RunAction::SetPrimaryParticleName(const G4String& name)
{
    fPrimaryParticleName = name;
}

void RunAction::SetEmModel(const G4String& model)
{
    fEmModel = model;
}

void RunAction::SetSampleThickness(G4double thickness)
{
    fSampleThickness = thickness;
}

void RunAction::SetSubstrateThickness(G4double thickness)
{
    fSubstrateThickness = thickness;
}

G4double RunAction::GetSampleThickness() const
{
    return fSampleThickness;
}

G4double RunAction::GetSubstrateThickness() const
{
    return fSubstrateThickness;
}

void RunAction::SetSeyAlphaInvNm(G4double alphaInvNm)
{
    fSeyAlphaInvNm = alphaInvNm;
}

G4double RunAction::GetSeyAlphaInvNm() const
{
    return fSeyAlphaInvNm;
}

void RunAction::SetPrimaryDirectionZ(G4double dirZ)
{
    fPrimaryDirectionZ = dirZ;
}

G4double RunAction::GetPrimaryDirectionZ() const
{
    return fPrimaryDirectionZ;
}

G4int RunAction::GetEdepPrimaryWeightedId() const
{
    return fEdepPrimaryWeightedId;
}

G4int RunAction::GetEdepDepthPrimaryId() const
{
    return fEdepDepthPrimaryId;
}

G4int RunAction::GetEdepDepthPrimaryWeightedId() const
{
    return fEdepDepthPrimaryWeightedId;
}

G4int RunAction::GetPrimaryTrackLengthId() const
{
    return fPrimaryTrackLengthId;
}

void RunAction::SetOutputTag(const G4String& tag)
{
    fOutputTag = tag;
}

void RunAction::SetPaiEnabled(G4bool enabled)
{
    fPaiEnabled = enabled;
}

void RunAction::SetLivermoreAtomicDeexcitation(G4int value)
{
    fLivermoreAtomicDeexcitation = value;
}

G4bool RunAction::IsPaiEnabled() const
{
    return fPaiEnabled;
}

void RunAction::OptimizeHistogramInFile(const G4String& fileName)
{
    // Open the ROOT file in update mode
    TFile* file = TFile::Open(fileName.c_str(), "UPDATE");
    if (!file || file->IsZombie()) {
        G4cerr << "Warning: Could not open file for histogram optimization: " << fileName << G4endl;
        return;
    }
    
    // Get the EdepPrimary histogram
    TH1* hOriginal = dynamic_cast<TH1*>(file->Get("EdepPrimary"));
    if (!hOriginal) {
        file->Close();
        return;
    }
    
    // Find the last non-empty bin to trim the range
    int lastNonEmptyBin = hOriginal->FindLastBinAbove(0.0);
    if (lastNonEmptyBin <= 0) {
        file->Close();
        return;
    }
    
    // Get the upper edge of the last non-empty bin
    double lastBinEdge = hOriginal->GetXaxis()->GetBinUpEdge(lastNonEmptyBin);
    
    // Add small padding (5% or at least 2 bin widths)
    const double binWidth = hOriginal->GetBinWidth(1);
    const double padding = std::max(lastBinEdge * 0.05, binWidth * 2.0);
    const double newMax = lastBinEdge + padding;
    const double currentMax = hOriginal->GetXaxis()->GetXmax();
    
    // Only optimize if we can trim significantly (at least 10% reduction)
    if (newMax < currentMax * 0.9) {
        // Create new histogram with trimmed range but same binning
        const int currentBins = hOriginal->GetNbinsX();
        const double currentMin = hOriginal->GetXaxis()->GetXmin();
        // Calculate new number of bins to maintain same bin width
        const int newBins = static_cast<int>(std::ceil((newMax - currentMin) / binWidth));
        
        TH1* hOptimized = new TH1D("EdepPrimary_optimized",
                                   hOriginal->GetTitle(),
                                   newBins, currentMin, newMax);
        hOptimized->SetDirectory(nullptr);
        
        // Copy axis titles
        hOptimized->GetXaxis()->SetTitle(hOriginal->GetXaxis()->GetTitle());
        hOptimized->GetYaxis()->SetTitle(hOriginal->GetYaxis()->GetTitle());
        
        // Copy bin contents from original (only up to newMax)
        for (int i = 1; i <= hOriginal->GetNbinsX(); i++) {
            const double binCenter = hOriginal->GetXaxis()->GetBinCenter(i);
            if (binCenter <= newMax) {
                const double content = hOriginal->GetBinContent(i);
                if (content > 0) {
                    // Find corresponding bin in new histogram
                    int newBin = hOptimized->GetXaxis()->FindBin(binCenter);
                    if (newBin >= 1 && newBin <= hOptimized->GetNbinsX()) {
                        hOptimized->SetBinContent(newBin, hOptimized->GetBinContent(newBin) + content);
                    }
                }
            }
        }
        
        // Update entries
        double totalContent = 0;
        for (int i = 1; i <= hOptimized->GetNbinsX(); i++) {
            totalContent += hOptimized->GetBinContent(i);
        }
        hOptimized->SetEntries(totalContent);
        
        // Replace the original histogram
        file->cd();
        file->Delete("EdepPrimary;*");
        hOptimized->SetName("EdepPrimary");
        hOptimized->SetDirectory(file);
        hOptimized->Write("EdepPrimary", TObject::kOverwrite);
        file->Flush();
    }
    
    file->Close();
}
