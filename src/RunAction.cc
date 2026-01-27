#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"

#include "DetectorConstruction.hh"

#include <cmath>

RunAction::RunAction()
    : G4UserRunAction(),
      fNPrimaryElectrons(0),
      fNSecondaryElectrons(0),
      fMinNonZeroEdep(-1.),
      fPrimaryEnergy(0.),
      fMaxPrimaryEnergy(0.),
      fPrimaryParticleName("e-"),
      fEmModel("PAI"),
      fSampleThickness(0.),
      fOutputTag("SEE_in_vacuum"),
      fPaiEnabled(false),
      fLivermoreAtomicDeexcitation(-1)
{
}

RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    fNPrimaryElectrons   = 0;
    fNSecondaryElectrons = 0;
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

    static G4bool metaCreated = false;
    static G4bool histosCreated = false;
    if (!metaCreated) {
        analysisManager->CreateNtuple("RunMeta", "Run metadata");
        analysisManager->CreateNtupleDColumn("primaryEnergyMeV");
        analysisManager->CreateNtupleDColumn("sampleThicknessNm");
        analysisManager->CreateNtupleDColumn("maxPrimaryEnergyMeV");
        analysisManager->CreateNtupleIColumn("paiEnabled");
        analysisManager->CreateNtupleSColumn("primaryParticle");
        analysisManager->CreateNtupleSColumn("emModel");
        analysisManager->CreateNtupleIColumn("livermoreAtomicDeexcitation");
        analysisManager->FinishNtuple();
        metaCreated = true;
    }

    if (histosCreated) {
        analysisManager->Reset();
    } else {
        // Create a 1D histogram for primary particle energy deposition in Al2O3
        // ID 0: EdepPrimary
        // Based on typical energy deposition: mean ~6 eV, RMS ~29 eV (for electrons)
        // Use range 0-200 eV with fine binning for better resolution
        // Note: Values will be filled in eV units (converted in EventAction)
        G4double maxEnergy = fMaxPrimaryEnergy;
        if (maxEnergy <= 0.) {
            maxEnergy = fPrimaryEnergy;
        }
        if (maxEnergy <= 0.) {
            maxEnergy = 200. * eV;
        }
        const G4int maxPrimaryBins = 200000;
        const G4double idealBins = std::ceil(maxEnergy / eV);
        const G4int primaryBins = std::max(
            1, static_cast<G4int>(std::min(idealBins, static_cast<G4double>(maxPrimaryBins))));

        G4int histoId = analysisManager->CreateH1(
            "EdepPrimary",
            "Primary particle energy deposition in Al_{2}O_{3}",
            primaryBins,   // ~1 eV per bin across scan max energy
            0.,    // Edep min (eV)
            maxEnergy / eV   // Edep max (eV)
        );

        // Set axis labels explicitly (X-axis in eV units)
        analysisManager->SetH1XAxisTitle(histoId, "Energy deposition (eV)");
        analysisManager->SetH1YAxisTitle(histoId, "Number of events");

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
        G4int residualId = analysisManager->CreateH1(
            "PrimaryResidualEnergy",
            "Primary residual kinetic energy at end of event",
            primaryBins, // 1 eV per bin across scan max energy
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

        // Create a 2D histogram: event edep vs number of steps in Al2O3
        // ID 0 for H2: EdepPrimaryVsSteps
        G4int edepVsStepsId = analysisManager->CreateH2(
            "EdepPrimaryVsSteps",
            "Primary energy deposition vs steps in Al_{2}O_{3}",
            primaryBins,
            0.,
            maxEnergy / eV,
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

    if (fNPrimaryElectrons > 0) {
        const G4double sey = static_cast<G4double>(fNSecondaryElectrons) /
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
        analysisManager->FillNtupleDColumn(2, fMaxPrimaryEnergy / MeV);
        analysisManager->FillNtupleIColumn(3, fPaiEnabled ? 1 : 0);
        analysisManager->FillNtupleSColumn(4, fPrimaryParticleName);
        analysisManager->FillNtupleSColumn(5, fEmModel);
        analysisManager->FillNtupleIColumn(6, fLivermoreAtomicDeexcitation);
        analysisManager->AddNtupleRow();

        analysisManager->Write();
        analysisManager->CloseFile();
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

