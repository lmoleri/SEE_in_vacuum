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
    if (std::abs(value - std::round(value)) < 1e-6) {
        oss << std::fixed << std::setprecision(0) << value;
        return oss.str();
    }
    oss << std::fixed << std::setprecision(3) << value;
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
}  // namespace

int main(int argc, char** argv)
{
    // Construct the default run manager
    G4RunManager* runManager = new G4RunManager;

    // Set mandatory initialization classes
    auto* detector = new DetectorConstruction();
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList());
    auto* actions = new ActionInitialization();
    runManager->SetUserInitialization(actions);

    // Initialize G4 kernel
    runManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    const bool hasArg = (argc != 1);
    const std::string arg1 = hasArg ? argv[1] : "";
    const bool isJsonScan = hasArg && EndsWith(arg1, ".json");

    // Initialize visualization (only for interactive or macro runs)
    G4VisManager* visManager = nullptr;
    if (!isJsonScan) {
        visManager = new G4VisExecutive;
        visManager->Initialize();
    }

    if (isJsonScan) {
        std::string content = StripWhitespace(ReadFile(arg1));
        auto thicknessNm = ParseArray(content, "sample_thickness_nm");
        auto energiesMeV = ParseArray(content, "primary_energy_MeV");
        double events = 100000;
        ParseNumber(content, "events", events);
        std::string outputDir = ParseString(content, "output_dir");

        if (thicknessNm.empty() || energiesMeV.empty()) {
            G4cerr << "Error: JSON must include arrays for sample_thickness_nm and primary_energy_MeV."
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

        std::string autoDir;
        if (outputDir.empty()) {
            std::string thickList = JoinParams(thicknessNm, "nm");
            std::string energyList = JoinParams(energiesMeV, "MeV");
            autoDir = "scan_thick" + thickList + "_energy" + energyList +
                      "_events" + FormatParam(events);
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

        for (double thickness : thicknessNm) {
            detector->SetSampleThickness(thickness * nm);
            runAction->SetSampleThickness(thickness * nm);
            runManager->ReinitializeGeometry(true);
            UImanager->ApplyCommand("/run/initialize");

            for (double energy : energiesMeV) {
                runAction->SetPrimaryEnergy(energy * MeV);
                G4String energyCmd = "/gun/energy " + std::to_string(energy) + " MeV";
                UImanager->ApplyCommand(energyCmd);

                std::string tag = outputDir + "/SEE_in_vacuum_thick" +
                                  FormatParam(thickness) + "nm_energy" +
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
