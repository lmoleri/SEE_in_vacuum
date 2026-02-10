#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "RunAction.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::string ReadFile(const std::string& path)
{
    std::ifstream file(path);
    std::ostringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

std::string StripWhitespace(std::string input)
{
    input.erase(std::remove_if(input.begin(), input.end(),
                               [](unsigned char c) { return std::isspace(c); }),
                input.end());
    return input;
}

std::vector<double> ParseArray(const std::string& content, const std::string& key)
{
    std::vector<double> values;
    const std::string searchKey = "\"" + key + "\"";
    auto keyPos = content.find(searchKey);
    if (keyPos == std::string::npos) {
        return values;
    }
    auto start = content.find('[', keyPos);
    auto end = content.find(']', start);
    if (start == std::string::npos || end == std::string::npos || end <= start) {
        return values;
    }
    std::string list = content.substr(start + 1, end - start - 1);
    std::stringstream ss(list);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) {
            values.push_back(std::stod(item));
        }
    }
    return values;
}

bool ParseNumber(const std::string& content, const std::string& key, double& out)
{
    const std::string searchKey = "\"" + key + "\"";
    auto keyPos = content.find(searchKey);
    if (keyPos == std::string::npos) {
        return false;
    }
    auto colon = content.find(':', keyPos);
    if (colon == std::string::npos) {
        return false;
    }
    auto end = content.find_first_of(",}", colon + 1);
    if (end == std::string::npos) {
        end = content.size();
    }
    std::string number = content.substr(colon + 1, end - colon - 1);
    if (number.empty()) {
        return false;
    }
    out = std::stod(number);
    return true;
}

std::string FormatParam(double value)
{
    std::ostringstream oss;
    const double absValue = std::abs(value);
    if (std::abs(value - std::round(value)) < 1e-6) {
        oss << std::fixed << std::setprecision(0) << value;
        return oss.str();
    }
    if (absValue > 0. && absValue < 0.001) {
        oss << std::fixed << std::setprecision(6) << value;
    } else {
        oss << std::fixed << std::setprecision(3) << value;
    }
    std::string out = oss.str();
    if (out.find('.') != std::string::npos) {
        out.erase(out.find_last_not_of('0') + 1, std::string::npos);
        if (!out.empty() && out.back() == '.') {
            out.pop_back();
        }
    }
    std::replace(out.begin(), out.end(), '.', 'p');
    return out;
}

bool EndsWith(const std::string& value, const std::string& suffix)
{
    if (value.size() < suffix.size()) {
        return false;
    }
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string JoinParams(const std::vector<double>& values, const std::string& suffix)
{
    std::ostringstream oss;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            oss << "-";
        }
        oss << FormatParam(values[i]) << suffix;
    }
    return oss.str();
}

std::string ParseString(const std::string& content, const std::string& key)
{
    const std::string searchKey = "\"" + key + "\"";
    auto keyPos = content.find(searchKey);
    if (keyPos == std::string::npos) {
        return "";
    }
    auto colon = content.find(':', keyPos);
    if (colon == std::string::npos) {
        return "";
    }
    auto firstQuote = content.find('"', colon + 1);
    if (firstQuote == std::string::npos) {
        return "";
    }
    auto secondQuote = content.find('"', firstQuote + 1);
    if (secondQuote == std::string::npos) {
        return "";
    }
    return content.substr(firstQuote + 1, secondQuote - firstQuote - 1);
}

std::string ToLower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return value;
}

bool ParseBool(const std::string& content, const std::string& key, bool& out)
{
    const std::string searchKey = "\"" + key + "\"";
    auto keyPos = content.find(searchKey);
    if (keyPos == std::string::npos) {
        return false;
    }
    auto colon = content.find(':', keyPos);
    if (colon == std::string::npos) {
        return false;
    }
    auto end = content.find_first_of(",}", colon + 1);
    if (end == std::string::npos) {
        end = content.size();
    }
    std::string value = content.substr(colon + 1, end - colon - 1);
    if (value.empty()) {
        return false;
    }
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (value == "true" || value == "1") {
        out = true;
        return true;
    }
    if (value == "false" || value == "0") {
        out = false;
        return true;
    }
    return false;
}
}  // namespace

int main(int argc, char** argv)
{
    const bool hasArg = (argc != 1);
    const std::string arg1 = hasArg ? argv[1] : "";
    const bool isJsonScan = hasArg && EndsWith(arg1, ".json");
    std::string jsonContent;
    bool paiOverride = false;
    bool paiEnabled = true;
    bool livermoreDeexcitationOverride = false;
    bool livermoreDeexcitationEnabled = true;
    bool atomicDeexcitationOverride = false;
    bool atomicDeexcitationEnabled = true;
    bool deexcitationIgnoreCutOverride = false;
    bool deexcitationIgnoreCutEnabled = true;
    std::string emModel = "PAI";
    std::string emModelLower = "pai";
    if (isJsonScan) {
        jsonContent = StripWhitespace(ReadFile(arg1));
        std::string emModelJson = ParseString(jsonContent, "em_model");
        if (!emModelJson.empty()) {
            emModel = emModelJson;
            emModelLower = ToLower(emModelJson);
        }
        if (ParseBool(jsonContent, "pai_enabled", paiEnabled)) {
            paiOverride = true;
        }
        bool disableDeexcitation = false;
        if (ParseBool(jsonContent, "disable_deexcitation", disableDeexcitation)) {
            atomicDeexcitationOverride = true;
            atomicDeexcitationEnabled = !disableDeexcitation;
        } else if (ParseBool(jsonContent, "atomic_deexcitation", atomicDeexcitationEnabled)) {
            atomicDeexcitationOverride = true;
        }
        if (ParseBool(jsonContent, "deexcitation_ignore_cut", deexcitationIgnoreCutEnabled)) {
            deexcitationIgnoreCutOverride = true;
        }
        if (ParseBool(jsonContent, "livermore_atomic_deexcitation", livermoreDeexcitationEnabled)) {
            livermoreDeexcitationOverride = true;
        }
    }

    // Construct the default run manager
    G4RunManager* runManager = new G4RunManager;

    // Set mandatory initialization classes
    auto* detector = new DetectorConstruction();
    runManager->SetUserInitialization(detector);
    auto* physicsList = new PhysicsList(emModel.c_str());
    if (paiOverride) {
        physicsList->SetPaiEnabledOverride(paiEnabled);
    }
    if (atomicDeexcitationOverride) {
        physicsList->SetAtomicDeexcitationOverride(atomicDeexcitationEnabled);
    }
    if (deexcitationIgnoreCutOverride) {
        physicsList->SetDeexcitationIgnoreCutOverride(deexcitationIgnoreCutEnabled);
    }
    if (livermoreDeexcitationOverride) {
        physicsList->SetLivermoreAtomicDeexcitationOverride(livermoreDeexcitationEnabled);
    }
    runManager->SetUserInitialization(physicsList);
    auto* actions = new ActionInitialization();
    runManager->SetUserInitialization(actions);

    // Initialize G4 kernel
    runManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Initialize visualization (only for interactive or macro runs)
    G4VisManager* visManager = nullptr;
    if (!isJsonScan) {
        visManager = new G4VisExecutive;
        visManager->Initialize();
    }

    if (isJsonScan) {
        std::string content = jsonContent;
        auto thicknessNm = ParseArray(content, "sample_thickness_nm");
        auto energiesMeV = ParseArray(content, "primary_energy_MeV");
        auto substrateNm = ParseArray(content, "substrate_thickness_nm");
        auto radiusNm = ParseArray(content, "sample_radius_nm");
        double maxStepNm = 0.0;
        bool maxStepOverride = ParseNumber(content, "max_step_nm", maxStepNm);
        if (!maxStepOverride) {
            maxStepOverride = ParseNumber(content, "al2o3_max_step_nm", maxStepNm);
        }
        if (substrateNm.empty()) {
            double substrateSingle = 0.0;
            if (ParseNumber(content, "substrate_thickness_nm", substrateSingle)) {
                substrateNm.push_back(substrateSingle);
            } else {
                substrateNm.push_back(0.0);
            }
        }
        if (radiusNm.empty()) {
            double radiusSingle = 100.0;
            if (ParseNumber(content, "sample_radius_nm", radiusSingle)) {
                radiusNm.push_back(radiusSingle);
            } else {
                radiusNm.push_back(100.0);
            }
        }
        double events = 100000;
        ParseNumber(content, "events", events);
        std::string outputDir = ParseString(content, "output_dir");
        std::string primaryParticle = ParseString(content, "primary_particle");
        double seyAlphaInvNm = 0.0;
        ParseNumber(content, "sey_alpha_inv_nm", seyAlphaInvNm);
        double verboseStepFrac = 0.0;
        ParseNumber(content, "verbose_step_fraction", verboseStepFrac);
        double verboseStepMax = 0.0;
        ParseNumber(content, "verbose_step_max", verboseStepMax);
        bool verboseSteps = false;
        bool verboseStepsSet = ParseBool(content, "verbose_step_diagnostics", verboseSteps);
        bool verboseStepsEnable = verboseStepsSet && verboseSteps;
        if (primaryParticle.empty()) {
            primaryParticle = "e-";
        }

        if (thicknessNm.empty() || energiesMeV.empty()) {
            G4cerr << "Error: JSON must include arrays for sample_thickness_nm and primary_energy_MeV."
                   << G4endl;
            delete runManager;
            return 1;
        }
        if (!(substrateNm.size() == 1 || substrateNm.size() == thicknessNm.size())) {
            G4cerr << "Error: substrate_thickness_nm must have size 1 or match sample_thickness_nm."
                   << G4endl;
            delete runManager;
            return 1;
        }
        if (!(radiusNm.size() == 1 || radiusNm.size() == thicknessNm.size())) {
            G4cerr << "Error: sample_radius_nm must have size 1 or match sample_thickness_nm."
                   << G4endl;
            delete runManager;
            return 1;
        }

        auto* runAction = actions->GetRunAction();
        if (!runAction) {
            G4cerr << "Error: RunAction not available for scan." << G4endl;
            delete runManager;
            return 1;
        }
        auto* primaryGenerator = actions->GetPrimaryGenerator();
        if (!primaryGenerator) {
            G4cerr << "Error: PrimaryGeneratorAction not available for scan." << G4endl;
            delete runManager;
            return 1;
        }

        primaryGenerator->SetParticleName(primaryParticle);
        runAction->SetPrimaryParticleName(primaryParticle);
        runAction->SetEmModel(emModel);
        if (maxStepOverride && maxStepNm > 0.) {
            detector->SetMaxStep(maxStepNm * nm);
            runAction->SetMaxStep(maxStepNm * nm);
        }
        if (verboseStepsEnable) {
            runAction->SetVerboseStepDiagnostics(true);
            if (verboseStepFrac > 0.0) {
                runAction->SetVerboseStepThresholdFrac(verboseStepFrac);
            }
            if (verboseStepMax > 0.0) {
                runAction->SetVerboseStepMaxCount(static_cast<G4int>(verboseStepMax));
            }
        }
        if (seyAlphaInvNm > 0.0) {
            runAction->SetSeyAlphaInvNm(seyAlphaInvNm);
        }
        const bool emModelIsPai = (emModelLower == "pai");
        const bool emModelIsLivermore =
            (emModelLower == "g4emlivermorephysics" || emModelLower == "livermore" ||
             emModelLower == "livermorephysics");
        const bool effectivePaiEnabled =
            emModelIsPai ? (paiOverride ? paiEnabled : true) : false;
        runAction->SetPaiEnabled(effectivePaiEnabled);
        const int livermoreDeexcitationMeta =
            emModelIsLivermore
                ? (livermoreDeexcitationOverride ? (livermoreDeexcitationEnabled ? 1 : 0) : -1)
                : -1;
        runAction->SetLivermoreAtomicDeexcitation(livermoreDeexcitationMeta);
        if (!energiesMeV.empty()) {
            const auto maxEnergy = *std::max_element(energiesMeV.begin(), energiesMeV.end());
            runAction->SetMaxPrimaryEnergy(maxEnergy * MeV);
        }

        std::string autoDir;
        const std::string substrateList = JoinParams(substrateNm, "nm");
        const std::string substrateSuffix = "_sub" + substrateList;
        const std::string radiusList = JoinParams(radiusNm, "nm");
        const std::string radiusSuffix = "_r" + radiusList;
        std::string stepSuffix;
        if (maxStepOverride && maxStepNm > 0.) {
            stepSuffix = "_step" + FormatParam(maxStepNm) + "nm";
        }
        if (outputDir.empty()) {
            std::string thickList = JoinParams(thicknessNm, "nm");
            std::string energyList = JoinParams(energiesMeV, "MeV");
            autoDir = "scan_thick" + thickList + "_particle" + primaryParticle +
                      "_energy" + energyList +
                      substrateSuffix +
                      radiusSuffix +
                      stepSuffix +
                      "_events" + FormatParam(events);
        } else if (outputDir.find("_sub") == std::string::npos) {
            outputDir += substrateSuffix;
        }
        if (outputDir.find("_r") == std::string::npos) {
            outputDir += radiusSuffix;
        }
        if (!stepSuffix.empty() && outputDir.find(stepSuffix) == std::string::npos) {
            outputDir += stepSuffix;
        }
        std::filesystem::path baseDir = std::filesystem::current_path();
        if (baseDir.filename() == "build") {
            baseDir = baseDir.parent_path();
        }
        if (outputDir.empty()) {
            outputDir = (baseDir / "results" / autoDir).string();
        } else {
            std::filesystem::path outPath(outputDir);
            if (outPath.is_relative()) {
                outPath = baseDir / outPath;
            }
            outputDir = outPath.string();
        }
        std::filesystem::create_directories(outputDir);

        const int eventsInt = static_cast<int>(events);

        for (size_t idx = 0; idx < thicknessNm.size(); ++idx) {
            const double thickness = thicknessNm[idx];
            const double substrate = (substrateNm.size() == 1) ? substrateNm[0] : substrateNm[idx];
            const double radius = (radiusNm.size() == 1) ? radiusNm[0] : radiusNm[idx];
            detector->SetSampleThickness(thickness * nm);
            detector->SetSubstrateThickness(substrate * nm);
            detector->SetSampleRadius(radius * nm);
            runAction->SetSampleThickness(thickness * nm);
            runAction->SetSubstrateThickness(substrate * nm);
            runAction->SetSampleRadius(radius * nm);
            runManager->ReinitializeGeometry(true);
            UImanager->ApplyCommand("/run/initialize");

            for (double energy : energiesMeV) {
                runAction->SetPrimaryEnergy(energy * MeV);
                G4String energyCmd = "/gun/energy " + std::to_string(energy) + " MeV";
                UImanager->ApplyCommand(energyCmd);

                std::string tag = outputDir + "/SEE_in_vacuum_thick" +
                                  FormatParam(thickness) + "nm_sub" +
                                  FormatParam(substrate) + "nm_r" +
                                  FormatParam(radius) + "nm_particle" +
                                  primaryParticle + "_energy" +
                                  FormatParam(energy) + "MeV_events" +
                                  FormatParam(events);
                runAction->SetOutputTag(tag.c_str());

                G4String runCmd = "/run/beamOn " + std::to_string(eventsInt);
                UImanager->ApplyCommand(runCmd);
            }
        }
    } else if (hasArg) {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    } else {
        // interactive mode : define UI session
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        // Allow running from build/ while keeping macros in the project root
        UImanager->ApplyCommand("/control/macroPath ..");
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    delete visManager;
    delete runManager;

    return 0;
}
