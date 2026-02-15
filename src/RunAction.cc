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

#include <algorithm>
#include <cmath>
#include <iomanip>

// ROOT headers for histogram optimization
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
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
      fSampleRadius(0.),
      fMaxStep(0.),
      fOutputTag("SEE_in_vacuum"),
      fPaiEnabled(false),
      fLivermoreAtomicDeexcitation(-1),
      fSeyAlphaInvNm(0.0),
      fPrimaryDirection(0., 0., 1.),
      fPrimaryDirectionZ(1.0),
      fSpecularAcceptanceEnabled(false),
      fSpecularAcceptanceHalfAngleDeg(5.0),
      fEdepPrimaryWeightedId(-1),
      fEdepDepthPrimaryId(-1),
      fEdepDepthPrimaryWeightedId(-1),
      fEdepDepthPrimaryCountsId(-1),
      fEdepStepDepthPrimaryId(-1),
      fPrimaryTrackLengthDepthId(-1),
      fPrimaryTrackLengthId(-1),
      fPrimaryExitClassId(-1),
      fPrimaryExitEnergyEntranceId(-1),
      fPrimaryExitEnergyEntranceSpecularId(-1),
      fPrimaryExitEnergyOppositeId(-1),
      fPrimaryExitEnergyLateralId(-1),
      fStepLengthAl2O3Id(-1),
      fEdepVsStepsId(-1),
      fResidualVsEndVolumeId(-1),
      fResidualVsLastProcessId(-1),
      fResidualVsStopStatusId(-1),
      fEdepPrimaryStopId(-1),
      fEdepPrimaryExitEntranceId(-1),
      fEdepPrimaryExitOppositeId(-1),
      fEdepPrimaryExitLateralId(-1),
      fEventDiagnosticsNtupleId(-1),
      fVerboseStepNtupleId(-1),
      fTrajectoryDiagnostics(false),
      fTrajectorySamplePerClass(300),
      fTrajectoryMaxStepsPerEvent(3000),
      fTrajectoryDiagnosticsNtupleId(-1),
      fTrajectoryClass2Used(0),
      fTrajectoryClass4Used(0)
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
    fVerboseStepUsed = 0;
    fTrajectoryClass2Used = 0;
    fTrajectoryClass4Used = 0;

    auto* analysisManager = G4AnalysisManager::Instance();

    // Basic configuration
    analysisManager->SetVerboseLevel(1);

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
    if (fSampleRadius <= 0.) {
        auto* det = dynamic_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        if (det) {
            fSampleRadius = det->GetSampleRadius();
        }
    }
    if (fMaxStep <= 0.) {
        auto* det = dynamic_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        if (det) {
            fMaxStep = det->GetMaxStep();
        }
    }

    static G4bool metaCreated = false;
    static G4bool eventDiagCreated = false;
    static G4int eventDiagNtupleId = -1;
    static G4bool histosCreated = false;
    static G4bool verboseCreated = false;
    static G4int verboseNtupleId = -1;
    static G4bool trajectoryCreated = false;
    static G4int trajectoryNtupleId = -1;
    if (!metaCreated) {
        analysisManager->SetFirstHistoId(0);
        analysisManager->CreateNtuple("RunMeta", "Run metadata");
        analysisManager->CreateNtupleDColumn("primaryEnergyMeV");
        analysisManager->CreateNtupleDColumn("sampleThicknessNm");
        analysisManager->CreateNtupleDColumn("substrateThicknessNm");
        analysisManager->CreateNtupleDColumn("sampleRadiusNm");
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
        analysisManager->CreateNtupleDColumn("maxStepNm");
        analysisManager->CreateNtupleDColumn("primaryDirX");
        analysisManager->CreateNtupleDColumn("primaryDirY");
        analysisManager->CreateNtupleDColumn("primaryDirZ");
        analysisManager->CreateNtupleDColumn("incidenceAngleSurfaceDeg");
        analysisManager->CreateNtupleDColumn("incidenceAngleNormalDeg");
        analysisManager->CreateNtupleIColumn("specularAcceptanceEnabled");
        analysisManager->CreateNtupleDColumn("specularAcceptanceHalfAngleDeg");
        analysisManager->FinishNtuple();
        metaCreated = true;
    }

    if (!eventDiagCreated) {
        eventDiagNtupleId = analysisManager->CreateNtuple(
            "EventDiagnostics",
            "Primary electron event-level diagnostics in Al2O3");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "eventId");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "primaryEnergyEv");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "edepPrimaryEv");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "primaryResidualEv");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "primaryExitClass");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "primaryExitEnergyEv");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "primaryStopStatus");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "primaryEndLocation");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nEdepSteps");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "primaryTrackLengthNm");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "maxDepthNm");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nBoundaryCrossings");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nDirectionReversalsZ");
        analysisManager->CreateNtupleSColumn(eventDiagNtupleId, "firstProcessInAl2O3");
        analysisManager->CreateNtupleSColumn(eventDiagNtupleId, "lastProcess");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "edepByEIoniEv");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "edepByMscEv");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "edepByOtherEv");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "edepFirstStepEv");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "edepMaxStepEv");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "depthFirstEdepNm");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "firstDirectionReversalStep");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "firstDirectionReversalDepthNm");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "firstDirectionReversalEnergyEv");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "firstBoundaryStep");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "firstBoundaryDepthNm");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "firstBoundaryEnergyEv");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "firstBoundaryType");
        analysisManager->CreateNtupleSColumn(eventDiagNtupleId, "firstDirectionReversalProcess");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "firstDirectionReversalStepStatus");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "firstDirectionReversalStepLenNm");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "firstDirectionReversalPreEnergyEv");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "firstDirectionReversalPostEnergyEv");
        analysisManager->CreateNtupleDColumn(eventDiagNtupleId, "firstDirectionReversalDeltaThetaDeg");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nStepStatusGeomBoundary");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nStepStatusPostStepProc");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nStepStatusAlongStepProc");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nStepStatusUserLimit");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nStepStatusOther");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nProcMsc");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nProcStepLimiter");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nProcTransportation");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nProcEIoni");
        analysisManager->CreateNtupleIColumn(eventDiagNtupleId, "nProcOther");
        analysisManager->FinishNtuple(eventDiagNtupleId);
        eventDiagCreated = true;
    }
    fEventDiagnosticsNtupleId = eventDiagNtupleId;

    if (fVerboseStepDiagnostics) {
        if (!verboseCreated) {
            verboseNtupleId = analysisManager->CreateNtuple("VerboseStepDiagnostics",
                                                            "High-fraction energy-deposit steps");
            analysisManager->CreateNtupleIColumn(verboseNtupleId, "eventId");
            analysisManager->CreateNtupleIColumn(verboseNtupleId, "trackId");
            analysisManager->CreateNtupleIColumn(verboseNtupleId, "stepNumber");
            analysisManager->CreateNtupleDColumn(verboseNtupleId, "depthNm");
            analysisManager->CreateNtupleDColumn(verboseNtupleId, "stepLenNm");
            analysisManager->CreateNtupleDColumn(verboseNtupleId, "edepEv");
            analysisManager->CreateNtupleDColumn(verboseNtupleId, "preEv");
            analysisManager->CreateNtupleDColumn(verboseNtupleId, "postEv");
            analysisManager->CreateNtupleDColumn(verboseNtupleId, "frac");
            analysisManager->CreateNtupleIColumn(verboseNtupleId, "stepStatus");
            analysisManager->CreateNtupleIColumn(verboseNtupleId, "trackStatus");
            analysisManager->CreateNtupleIColumn(verboseNtupleId, "preVolCode");
            analysisManager->CreateNtupleIColumn(verboseNtupleId, "postVolCode");
            analysisManager->CreateNtupleIColumn(verboseNtupleId, "procCode");
            analysisManager->CreateNtupleSColumn(verboseNtupleId, "preVol");
            analysisManager->CreateNtupleSColumn(verboseNtupleId, "postVol");
            analysisManager->CreateNtupleSColumn(verboseNtupleId, "process");
            analysisManager->FinishNtuple(verboseNtupleId);
            verboseCreated = true;
        }
        fVerboseStepNtupleId = verboseNtupleId;
    } else {
        fVerboseStepNtupleId = -1;
    }

    if (fTrajectoryDiagnostics) {
        if (!trajectoryCreated) {
            trajectoryNtupleId = analysisManager->CreateNtuple(
                "PrimaryTrajectoryDiagnostics",
                "Sampled full step-by-step primary electron trajectories in Al2O3");
            analysisManager->CreateNtupleIColumn(trajectoryNtupleId, "eventId");
            analysisManager->CreateNtupleIColumn(trajectoryNtupleId, "sampleIndex");
            analysisManager->CreateNtupleIColumn(trajectoryNtupleId, "primaryExitClass");
            analysisManager->CreateNtupleIColumn(trajectoryNtupleId, "stepNumber");
            analysisManager->CreateNtupleDColumn(trajectoryNtupleId, "preDepthNm");
            analysisManager->CreateNtupleDColumn(trajectoryNtupleId, "postDepthNm");
            analysisManager->CreateNtupleDColumn(trajectoryNtupleId, "stepLenNm");
            analysisManager->CreateNtupleDColumn(trajectoryNtupleId, "preEnergyEv");
            analysisManager->CreateNtupleDColumn(trajectoryNtupleId, "postEnergyEv");
            analysisManager->CreateNtupleDColumn(trajectoryNtupleId, "edepEv");
            analysisManager->CreateNtupleDColumn(trajectoryNtupleId, "dirZPre");
            analysisManager->CreateNtupleDColumn(trajectoryNtupleId, "dirZPost");
            analysisManager->CreateNtupleDColumn(trajectoryNtupleId, "deltaThetaDeg");
            analysisManager->CreateNtupleIColumn(trajectoryNtupleId, "reversalOnStep");
            analysisManager->CreateNtupleSColumn(trajectoryNtupleId, "process");
            analysisManager->CreateNtupleIColumn(trajectoryNtupleId, "stepStatus");
            analysisManager->CreateNtupleSColumn(trajectoryNtupleId, "preVol");
            analysisManager->CreateNtupleSColumn(trajectoryNtupleId, "postVol");
            analysisManager->CreateNtupleIColumn(trajectoryNtupleId, "isBoundaryCrossing");
            analysisManager->CreateNtupleIColumn(trajectoryNtupleId, "isOutwardBoundary");
            analysisManager->CreateNtupleIColumn(trajectoryNtupleId, "isFirstReversalStep");
            analysisManager->FinishNtuple(trajectoryNtupleId);
            trajectoryCreated = true;
        }
        fTrajectoryDiagnosticsNtupleId = trajectoryNtupleId;
    } else {
        fTrajectoryDiagnosticsNtupleId = -1;
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


        // Get maxEnergy for other histograms (use scan max when available).
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

        // Primary e- exits from Al2O3 -> World classified by side:
        // 1: entrance-side exit (backscatter-like), 2: opposite-side exit (transmission-like),
        // 3: lateral/edge exit.
        fPrimaryExitClassId = analysisManager->CreateH1(
            "PrimaryExitClass",
            "Primary e- exit class at Al_{2}O_{3}#rightarrowWorld",
            4,
            0.,
            4.
        );
        analysisManager->SetH1XAxisTitle(fPrimaryExitClassId, "Exit class");
        analysisManager->SetH1YAxisTitle(fPrimaryExitClassId, "Number of exits");

        fPrimaryExitEnergyEntranceId = analysisManager->CreateH1(
            "PrimaryExitEnergyEntrance",
            "Primary e- exit kinetic energy (entrance-side exit)",
            residualBins,
            0.,
            maxEnergy / eV
        );
        analysisManager->SetH1XAxisTitle(fPrimaryExitEnergyEntranceId, "Exit kinetic energy (eV)");
        analysisManager->SetH1YAxisTitle(fPrimaryExitEnergyEntranceId, "Number of exits");

        fPrimaryExitEnergyEntranceSpecularId = analysisManager->CreateH1(
            "PrimaryExitEnergyEntranceSpecular",
            "Primary e- exit kinetic energy (entrance-side, specular-cone accepted)",
            residualBins,
            0.,
            maxEnergy / eV
        );
        analysisManager->SetH1XAxisTitle(
            fPrimaryExitEnergyEntranceSpecularId, "Exit kinetic energy (eV)");
        analysisManager->SetH1YAxisTitle(
            fPrimaryExitEnergyEntranceSpecularId, "Number of accepted exits");

        fPrimaryExitEnergyOppositeId = analysisManager->CreateH1(
            "PrimaryExitEnergyOpposite",
            "Primary e- exit kinetic energy (opposite-side exit)",
            residualBins,
            0.,
            maxEnergy / eV
        );
        analysisManager->SetH1XAxisTitle(fPrimaryExitEnergyOppositeId, "Exit kinetic energy (eV)");
        analysisManager->SetH1YAxisTitle(fPrimaryExitEnergyOppositeId, "Number of exits");

        fPrimaryExitEnergyLateralId = analysisManager->CreateH1(
            "PrimaryExitEnergyLateral",
            "Primary e- exit kinetic energy (lateral/edge exit)",
            residualBins,
            0.,
            maxEnergy / eV
        );
        analysisManager->SetH1XAxisTitle(fPrimaryExitEnergyLateralId, "Exit kinetic energy (eV)");
        analysisManager->SetH1YAxisTitle(fPrimaryExitEnergyLateralId, "Number of exits");

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
        fStepLengthAl2O3Id = analysisManager->CreateH1(
            "StepLengthAl2O3",
            "Step length in Al_{2}O_{3}",
            stepLenBins,  // ~0.1 nm bins up to capped max
            0.,
            maxStepLenNm
        );
        analysisManager->SetH1XAxisTitle(fStepLengthAl2O3Id, "Step length (nm)");
        analysisManager->SetH1YAxisTitle(fStepLengthAl2O3Id, "Number of steps");

        // Create depth profiles for primary electron energy deposition (unweighted + weighted).
        // Tie binning to the enforced max step when available so smaller steps yield finer depth resolution.
        if (thicknessNm > 0.) {
            const G4double maxStepNm = (fMaxStep > 0.) ? (fMaxStep / nm) : 0.5;
            const G4double depthBinNm = (maxStepNm > 0.)
                                            ? std::max(0.005, maxStepNm)
                                            : 0.05;
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

            fEdepDepthPrimaryCountsId = analysisManager->CreateH1(
                "EdepDepthPrimaryCounts",
                "Primary e- energy-depositing step count vs depth in Al_{2}O_{3}",
                depthBins,
                0.,
                thicknessNm
            );
            analysisManager->SetH1XAxisTitle(fEdepDepthPrimaryCountsId, "Depth from entrance (nm)");
            analysisManager->SetH1YAxisTitle(fEdepDepthPrimaryCountsId, "Energy-depositing steps");

            // 2D: energy deposited per step vs depth (primary e- only)
            const G4double maxStepEdepEv =
                std::max(50.0, std::min(maxEnergy / eV, 50000.0));
            const G4int maxStepBins2D = 2000;
            const G4double idealStepBins2D = std::ceil(maxStepEdepEv);
            const G4int stepBins2D = std::max(
                1, static_cast<G4int>(std::min(idealStepBins2D, static_cast<G4double>(maxStepBins2D))));
            const G4int stepDepthBins = (maxStepNm > 0.)
                                            ? std::max(1, static_cast<G4int>(std::ceil(thicknessNm / maxStepNm)))
                                            : depthBins;
            fEdepStepDepthPrimaryId = analysisManager->CreateH2(
                "EdepStepVsDepthPrimary",
                "Primary e- energy deposition per step vs depth in Al_{2}O_{3}",
                stepDepthBins,
                0.,
                thicknessNm,
                stepBins2D,
                0.,
                maxStepEdepEv
            );
            analysisManager->SetH2XAxisTitle(fEdepStepDepthPrimaryId, "Depth from entrance (nm)");
            analysisManager->SetH2YAxisTitle(fEdepStepDepthPrimaryId, "Energy deposition per step (eV)");

            fPrimaryTrackLengthDepthId = analysisManager->CreateH1(
                "PrimaryTrackLengthDepth",
                "Primary e- track length vs depth in Al_{2}O_{3}",
                depthBins,
                0.,
                thicknessNm
            );
            analysisManager->SetH1XAxisTitle(fPrimaryTrackLengthDepthId, "Depth from entrance (nm)");
            analysisManager->SetH1YAxisTitle(fPrimaryTrackLengthDepthId, "Track length (nm)");
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
        fEdepVsStepsId = analysisManager->CreateH2(
            "EdepPrimaryVsSteps",
            "Primary energy deposition vs steps in Al_{2}O_{3}",
            finalBins,  // Use same bins as EdepPrimary
            0.,
            maxEdepRange / eV,  // Use same range as EdepPrimary
            stepsMax,
            0.,
            stepsMax
        );
        analysisManager->SetH2XAxisTitle(fEdepVsStepsId, "Primary energy deposition (eV)");
        analysisManager->SetH2YAxisTitle(fEdepVsStepsId, "Energy-depositing steps per event");

        // Create a 2D histogram: primary residual energy vs end volume category
        // ID 1 for H2: ResidualEnergyVsEndVolume
        fResidualVsEndVolumeId = analysisManager->CreateH2(
            "ResidualEnergyVsEndVolume",
            "Primary residual energy vs end volume",
            primaryBins,
            0.,
            maxEnergy / eV,
            5,
            0.,
            5.
        );
        analysisManager->SetH2XAxisTitle(fResidualVsEndVolumeId, "Residual kinetic energy (eV)");
        analysisManager->SetH2YAxisTitle(fResidualVsEndVolumeId, "End volume category");

        // Create a 2D histogram: primary residual energy vs last process category
        // ID 2 for H2: ResidualEnergyVsLastProcess
        fResidualVsLastProcessId = analysisManager->CreateH2(
            "ResidualEnergyVsLastProcess",
            "Primary residual energy vs last process",
            primaryBins,
            0.,
            maxEnergy / eV,
            7,
            0.,
            7.
        );
        analysisManager->SetH2XAxisTitle(fResidualVsLastProcessId, "Residual kinetic energy (eV)");
        analysisManager->SetH2YAxisTitle(fResidualVsLastProcessId, "Last process category");

        // Create a 2D histogram: primary residual energy vs stop status category
        // ID 3 for H2: ResidualEnergyVsStopStatus
        fResidualVsStopStatusId = analysisManager->CreateH2(
            "ResidualEnergyVsStopStatus",
            "Primary residual energy vs stop status",
            primaryBins,
            0.,
            maxEnergy / eV,
            7,
            0.,
            7.
        );
        analysisManager->SetH2XAxisTitle(fResidualVsStopStatusId, "Residual kinetic energy (eV)");
        analysisManager->SetH2YAxisTitle(fResidualVsStopStatusId, "Stop status category");

        // Class-conditioned primary edep spectra for event topology diagnostics.
        fEdepPrimaryStopId = analysisManager->CreateH1(
            "EdepPrimaryStop",
            "Primary edep in Al_{2}O_{3} (stop/no-exit events)",
            finalBins,
            0.,
            maxEdepRange / eV
        );
        analysisManager->SetH1XAxisTitle(fEdepPrimaryStopId, "Primary energy deposition (eV)");
        analysisManager->SetH1YAxisTitle(fEdepPrimaryStopId, "Number of events");

        fEdepPrimaryExitEntranceId = analysisManager->CreateH1(
            "EdepPrimaryExitEntrance",
            "Primary edep in Al_{2}O_{3} (entrance-side exit events)",
            finalBins,
            0.,
            maxEdepRange / eV
        );
        analysisManager->SetH1XAxisTitle(fEdepPrimaryExitEntranceId, "Primary energy deposition (eV)");
        analysisManager->SetH1YAxisTitle(fEdepPrimaryExitEntranceId, "Number of events");

        fEdepPrimaryExitOppositeId = analysisManager->CreateH1(
            "EdepPrimaryExitOpposite",
            "Primary edep in Al_{2}O_{3} (opposite-side exit events)",
            finalBins,
            0.,
            maxEdepRange / eV
        );
        analysisManager->SetH1XAxisTitle(fEdepPrimaryExitOppositeId, "Primary energy deposition (eV)");
        analysisManager->SetH1YAxisTitle(fEdepPrimaryExitOppositeId, "Number of events");

        fEdepPrimaryExitLateralId = analysisManager->CreateH1(
            "EdepPrimaryExitLateral",
            "Primary edep in Al_{2}O_{3} (lateral-exit events)",
            finalBins,
            0.,
            maxEdepRange / eV
        );
        analysisManager->SetH1XAxisTitle(fEdepPrimaryExitLateralId, "Primary energy deposition (eV)");
        analysisManager->SetH1YAxisTitle(fEdepPrimaryExitLateralId, "Number of events");

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
    if (fTrajectoryDiagnostics) {
        G4cout << "  Trajectory diagnostics (sampled) class-2 opposite exits : "
               << fTrajectoryClass2Used << " / " << fTrajectorySamplePerClass << G4endl;
        G4cout << "  Trajectory diagnostics (sampled) class-4 trapped events : "
               << fTrajectoryClass4Used << " / " << fTrajectorySamplePerClass << G4endl;
    }
    G4cout << "==========================\n" << G4endl;

    // Write and close analysis output
    auto* analysisManager = G4AnalysisManager::Instance();
    if (analysisManager) {
        analysisManager->FillNtupleDColumn(0, fPrimaryEnergy / MeV);
        analysisManager->FillNtupleDColumn(1, fSampleThickness / nm);
        analysisManager->FillNtupleDColumn(2, fSubstrateThickness / nm);
        analysisManager->FillNtupleDColumn(3, fSampleRadius / nm);
        analysisManager->FillNtupleDColumn(4, fMaxPrimaryEnergy / MeV);
        analysisManager->FillNtupleIColumn(5, fPaiEnabled ? 1 : 0);
        analysisManager->FillNtupleSColumn(6, fPrimaryParticleName);
        analysisManager->FillNtupleSColumn(7, fEmModel);
        analysisManager->FillNtupleIColumn(8, fLivermoreAtomicDeexcitation);
        analysisManager->FillNtupleIColumn(9, fNPrimaryElectrons);
        analysisManager->FillNtupleIColumn(10, fNSecondaryElectrons);
        analysisManager->FillNtupleIColumn(11, fNEmittedElectrons);
        const G4double sey = (fNPrimaryElectrons > 0)
                                 ? static_cast<G4double>(fNEmittedElectrons) /
                                       static_cast<G4double>(fNPrimaryElectrons)
                                 : 0.0;
        analysisManager->FillNtupleDColumn(12, sey);
        const G4double minNonZeroEdepEv = (fMinNonZeroEdep > 0.) ? (fMinNonZeroEdep / eV) : 0.0;
        analysisManager->FillNtupleDColumn(13, minNonZeroEdepEv);
        analysisManager->FillNtupleDColumn(14, fMaxStep / nm);
        const G4ThreeVector dirUnit = (fPrimaryDirection.mag2() > 0.)
                                          ? fPrimaryDirection.unit()
                                          : G4ThreeVector(0., 0., 1.);
        G4double absDz = std::abs(dirUnit.z());
        if (absDz > 1.0) absDz = 1.0;
        if (absDz < 0.0) absDz = 0.0;
        const G4double incidenceSurfaceDeg = std::asin(absDz) / deg;
        const G4double incidenceNormalDeg = 90.0 - incidenceSurfaceDeg;
        analysisManager->FillNtupleDColumn(15, dirUnit.x());
        analysisManager->FillNtupleDColumn(16, dirUnit.y());
        analysisManager->FillNtupleDColumn(17, dirUnit.z());
        analysisManager->FillNtupleDColumn(18, incidenceSurfaceDeg);
        analysisManager->FillNtupleDColumn(19, incidenceNormalDeg);
        analysisManager->FillNtupleIColumn(20, fSpecularAcceptanceEnabled ? 1 : 0);
        analysisManager->FillNtupleDColumn(21, fSpecularAcceptanceHalfAngleDeg);
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

void RunAction::SetSampleRadius(G4double radius)
{
    fSampleRadius = radius;
}

void RunAction::SetMaxStep(G4double maxStep)
{
    fMaxStep = maxStep;
}

G4double RunAction::GetSampleThickness() const
{
    return fSampleThickness;
}

G4double RunAction::GetSubstrateThickness() const
{
    return fSubstrateThickness;
}

G4double RunAction::GetSampleRadius() const
{
    return fSampleRadius;
}

G4double RunAction::GetMaxStep() const
{
    return fMaxStep;
}

void RunAction::SetSeyAlphaInvNm(G4double alphaInvNm)
{
    fSeyAlphaInvNm = alphaInvNm;
}

G4double RunAction::GetSeyAlphaInvNm() const
{
    return fSeyAlphaInvNm;
}

void RunAction::SetPrimaryDirection(const G4ThreeVector& dir)
{
    if (dir.mag2() <= 0.) {
        return;
    }
    fPrimaryDirection = dir.unit();
    fPrimaryDirectionZ = fPrimaryDirection.z();
}

G4ThreeVector RunAction::GetPrimaryDirection() const
{
    return fPrimaryDirection;
}

void RunAction::SetPrimaryDirectionZ(G4double dirZ)
{
    fPrimaryDirectionZ = dirZ;
    if (fPrimaryDirection.mag2() <= 0.) {
        fPrimaryDirection = G4ThreeVector(0., 0., 1.);
    }
    const G4double clamped = std::max(-1.0, std::min(1.0, static_cast<double>(dirZ)));
    if (std::abs(clamped) > 1e-12) {
        fPrimaryDirection.setZ(clamped);
        fPrimaryDirection = fPrimaryDirection.unit();
    }
}

G4double RunAction::GetPrimaryDirectionZ() const
{
    return fPrimaryDirectionZ;
}

void RunAction::SetSpecularAcceptance(G4bool enabled, G4double halfAngleDeg)
{
    fSpecularAcceptanceEnabled = enabled;
    if (halfAngleDeg < 0.) {
        halfAngleDeg = 0.;
    }
    if (halfAngleDeg > 180.) {
        halfAngleDeg = 180.;
    }
    fSpecularAcceptanceHalfAngleDeg = halfAngleDeg;
}

G4bool RunAction::IsSpecularAcceptanceEnabled() const
{
    return fSpecularAcceptanceEnabled;
}

G4double RunAction::GetSpecularAcceptanceHalfAngleDeg() const
{
    return fSpecularAcceptanceHalfAngleDeg;
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

G4int RunAction::GetEdepDepthPrimaryCountsId() const
{
    return fEdepDepthPrimaryCountsId;
}

G4int RunAction::GetEdepStepDepthPrimaryId() const
{
    return fEdepStepDepthPrimaryId;
}

G4int RunAction::GetPrimaryTrackLengthDepthId() const
{
    return fPrimaryTrackLengthDepthId;
}

G4int RunAction::GetPrimaryTrackLengthId() const
{
    return fPrimaryTrackLengthId;
}

G4int RunAction::GetPrimaryExitClassId() const
{
    return fPrimaryExitClassId;
}

G4int RunAction::GetPrimaryExitEnergyEntranceId() const
{
    return fPrimaryExitEnergyEntranceId;
}

G4int RunAction::GetPrimaryExitEnergyEntranceSpecularId() const
{
    return fPrimaryExitEnergyEntranceSpecularId;
}

G4int RunAction::GetPrimaryExitEnergyOppositeId() const
{
    return fPrimaryExitEnergyOppositeId;
}

G4int RunAction::GetPrimaryExitEnergyLateralId() const
{
    return fPrimaryExitEnergyLateralId;
}

G4int RunAction::GetStepLengthAl2O3Id() const
{
    return fStepLengthAl2O3Id;
}

G4int RunAction::GetEdepVsStepsId() const
{
    return fEdepVsStepsId;
}

G4int RunAction::GetResidualVsEndVolumeId() const
{
    return fResidualVsEndVolumeId;
}

G4int RunAction::GetResidualVsLastProcessId() const
{
    return fResidualVsLastProcessId;
}

G4int RunAction::GetResidualVsStopStatusId() const
{
    return fResidualVsStopStatusId;
}

G4int RunAction::GetEdepPrimaryStopId() const
{
    return fEdepPrimaryStopId;
}

G4int RunAction::GetEdepPrimaryExitEntranceId() const
{
    return fEdepPrimaryExitEntranceId;
}

G4int RunAction::GetEdepPrimaryExitOppositeId() const
{
    return fEdepPrimaryExitOppositeId;
}

G4int RunAction::GetEdepPrimaryExitLateralId() const
{
    return fEdepPrimaryExitLateralId;
}

G4int RunAction::GetEventDiagnosticsNtupleId() const
{
    return fEventDiagnosticsNtupleId;
}

G4double RunAction::GetPrimaryEnergy() const
{
    return fPrimaryEnergy;
}

G4int RunAction::GetVerboseStepNtupleId() const
{
    return fVerboseStepNtupleId;
}

G4bool RunAction::IsTrajectoryDiagnostics() const
{
    return fTrajectoryDiagnostics;
}

G4int RunAction::GetTrajectoryDiagnosticsNtupleId() const
{
    return fTrajectoryDiagnosticsNtupleId;
}

G4int RunAction::GetTrajectorySamplePerClass() const
{
    return fTrajectorySamplePerClass;
}

G4int RunAction::GetTrajectoryMaxStepsPerEvent() const
{
    return fTrajectoryMaxStepsPerEvent;
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

G4bool RunAction::IsVerboseStepDiagnostics() const
{
    return fVerboseStepDiagnostics;
}

G4double RunAction::GetVerboseStepThresholdFrac() const
{
    return fVerboseStepThresholdFrac;
}

G4int RunAction::GetVerboseStepMaxCount() const
{
    return fVerboseStepMaxCount;
}

void RunAction::SetVerboseStepDiagnostics(G4bool enabled)
{
    fVerboseStepDiagnostics = enabled;
}

void RunAction::SetVerboseStepThresholdFrac(G4double frac)
{
    if (frac > 0.0 && frac <= 1.0) {
        fVerboseStepThresholdFrac = frac;
    }
}

void RunAction::SetVerboseStepMaxCount(G4int maxCount)
{
    if (maxCount > 0) {
        fVerboseStepMaxCount = maxCount;
    }
}

void RunAction::SetTrajectoryDiagnostics(G4bool enabled)
{
    fTrajectoryDiagnostics = enabled;
}

void RunAction::SetTrajectorySamplePerClass(G4int maxCount)
{
    if (maxCount > 0) {
        fTrajectorySamplePerClass = maxCount;
    }
}

void RunAction::SetTrajectoryMaxStepsPerEvent(G4int maxCount)
{
    if (maxCount > 0) {
        fTrajectoryMaxStepsPerEvent = maxCount;
    }
}

G4int RunAction::ConsumeVerboseStepSlot()
{
    if (!fVerboseStepDiagnostics) {
        return -1;
    }
    if (fVerboseStepUsed >= fVerboseStepMaxCount) {
        return -1;
    }
    return fVerboseStepUsed++;
}

G4bool RunAction::AcquireTrajectorySampleSlot(G4int exitClass, G4int& sampleIndex)
{
    sampleIndex = -1;
    if (!fTrajectoryDiagnostics) {
        return false;
    }
    if (exitClass == 2) {
        if (fTrajectoryClass2Used >= fTrajectorySamplePerClass) {
            return false;
        }
        sampleIndex = fTrajectoryClass2Used;
        ++fTrajectoryClass2Used;
        return true;
    }
    if (exitClass == 4) {
        if (fTrajectoryClass4Used >= fTrajectorySamplePerClass) {
            return false;
        }
        sampleIndex = fTrajectoryClass4Used;
        ++fTrajectoryClass4Used;
        return true;
    }
    return false;
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
    
    // Ensure EdepStepVsDepthPrimary Y-axis is limited to the primary energy for this run.
    TH2* hStepDepth = dynamic_cast<TH2*>(file->Get("EdepStepVsDepthPrimary"));
    if (hStepDepth) {
        const double maxEnergyEv = (fPrimaryEnergy > 0.)
                                       ? (fPrimaryEnergy / eV)
                                       : ((fMaxPrimaryEnergy > 0.) ? (fMaxPrimaryEnergy / eV) : 0.0);
        if (maxEnergyEv > 0.) {
            hStepDepth->GetYaxis()->SetRangeUser(0.0, maxEnergyEv);
            hStepDepth->Write("EdepStepVsDepthPrimary", TObject::kOverwrite);
            file->Flush();
        }
        if (fNPrimaryElectrons > 0) {
            auto* hStepDepthPerEvent =
                dynamic_cast<TH2*>(hStepDepth->Clone("EdepStepVsDepthPrimaryPerEvent"));
            if (hStepDepthPerEvent) {
                hStepDepthPerEvent->Scale(1.0 / static_cast<double>(fNPrimaryElectrons));
                hStepDepthPerEvent->Write("EdepStepVsDepthPrimaryPerEvent", TObject::kOverwrite);
                file->Flush();
            }
        }
    }

    if (fNPrimaryElectrons > 0) {
        const double norm = 1.0 / static_cast<double>(fNPrimaryElectrons);
        const char* histNames[] = {
            "EdepDepthPrimary",
            "EdepDepthPrimaryWeighted",
            "EdepDepthPrimaryCounts"
        };
        for (int i = 0; i < 3; ++i) {
            TH1* h = dynamic_cast<TH1*>(file->Get(histNames[i]));
            if (!h) {
                continue;
            }
            auto* hPerEvent = dynamic_cast<TH1*>(h->Clone(histNames[i]));
            if (!hPerEvent) {
                continue;
            }
            hPerEvent->Scale(norm);
            hPerEvent->Write(histNames[i], TObject::kOverwrite);
            file->Flush();
        }
        file->Delete("EdepDepthPrimaryPerEvent;*");
        file->Delete("EdepDepthPrimaryWeightedPerEvent;*");
        file->Delete("EdepDepthPrimaryCountsPerEvent;*");
    }

    file->Close();
}
