#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RegionStore.hh"
#include "G4UserLimits.hh"

DetectorConstruction::DetectorConstruction()
    : fWorldMaterial(nullptr),
      fAl2O3Material(nullptr),
      fSiMaterial(nullptr),
      fWorldLogical(nullptr),
      fAl2O3Logical(nullptr),
      fSiLogical(nullptr),
      fWorldPhysical(nullptr),
      fAl2O3Physical(nullptr),
      fSiPhysical(nullptr),
      fAl2O3Region(nullptr)
{
    fSampleThickness = 20.0 * nm;
    fSubstrateThickness = 0.0;
    fSampleRadius = 100.0 * nm;
    fMaxStep = 0.5 * nm;
}

DetectorConstruction::~DetectorConstruction()
{
}

void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();

    // World material - vacuum
    fWorldMaterial = nist->FindOrBuildMaterial("G4_Galactic");

    // Al2O3 (Aluminum Oxide) material
    // Density: 3.95 g/cm³
    fAl2O3Material = G4Material::GetMaterial("Al2O3");
    if (!fAl2O3Material) {
        G4double density = 3.95 * g/cm3;
        G4double a;  // atomic mass
        G4double z;  // atomic number

        // Create Al2O3 material
        fAl2O3Material = new G4Material("Al2O3", density, 2);

        // Aluminum
        a = 26.98 * g/mole;
        G4Element* elAl = new G4Element("Aluminum", "Al", z=13, a);

        // Oxygen
        a = 16.00 * g/mole;
        G4Element* elO = new G4Element("Oxygen", "O", z=8, a);

        // Add elements to Al2O3 (2 Al atoms, 3 O atoms)
        fAl2O3Material->AddElement(elAl, 2);
        fAl2O3Material->AddElement(elO, 3);
    }

    // Silicon substrate material
    fSiMaterial = nist->FindOrBuildMaterial("G4_Si");
    
    G4cout << "\n--- Material properties ---" << G4endl;
    G4cout << "World material: " << fWorldMaterial->GetName() << G4endl;
    G4cout << "Al2O3 density: " << fAl2O3Material->GetDensity() / (g/cm3) << " g/cm³" << G4endl;
    if (fSiMaterial) {
        G4cout << "Si density: " << fSiMaterial->GetDensity() / (g/cm3) << " g/cm³" << G4endl;
    }
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    DefineMaterials();
    return DefineVolumes();
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    // World volume - large box
    G4double worldSize = 1.0 * m;
    G4Box* worldSolid = new G4Box("World", 
                                   0.5 * worldSize, 
                                   0.5 * worldSize, 
                                   0.5 * worldSize);
    
    fWorldLogical = new G4LogicalVolume(worldSolid, 
                                         fWorldMaterial, 
                                         "World");
    
    // World physical volume (root volume)
    fWorldPhysical = new G4PVPlacement(
        nullptr,                    // no rotation
        G4ThreeVector(0, 0, 0),     // position
        fWorldLogical,              // logical volume
        "World",                    // name
        nullptr,                    // no mother volume (this is the world)
        false,                      // no boolean operations
        0,                          // copy number
        true);                      // check overlaps

    // Al2O3 layer - cylindrical disk
    // Thickness and diameter are configurable (defaults: 20 nm, 200 nm)
    G4double thickness = fSampleThickness > 0. ? fSampleThickness : 20.0 * nm;
    G4double radius = (fSampleRadius > 0.) ? fSampleRadius : 100.0 * nm;
    G4double substrateThickness = fSubstrateThickness;

    fSampleThickness = thickness;
    fSampleRadius = radius;
    
    G4Tubs* al2o3Solid = new G4Tubs("Al2O3",
                                    0.0,           // inner radius
                                    radius,        // outer radius
                                    0.5 * thickness,  // half-length in z
                                    0.0,           // start angle
                                    2.0 * M_PI);   // span angle
    
    fAl2O3Logical = new G4LogicalVolume(al2o3Solid,
                                        fAl2O3Material,
                                        "Al2O3");

    // Constrain step length inside Al2O3 for depth-resolved energy deposition.
    // Default is 0.5 nm unless overridden by configuration.
    if (fMaxStep > 0.) {
        fAl2O3Logical->SetUserLimits(new G4UserLimits(fMaxStep));
    }
    
    // Position the Al2O3 layer at the origin, inside the world
    fAl2O3Physical = new G4PVPlacement(
        nullptr,                    // no rotation
        G4ThreeVector(0, 0, 0),     // position
        fAl2O3Logical,              // logical volume
        "Al2O3",                    // name
        fWorldLogical,              // mother volume
        false,                      // no boolean operations
        0,                          // copy number
        true);                      // check overlaps

    // Silicon substrate below the Al2O3 coating (positive z direction)
    if (fSiMaterial && substrateThickness > 0.) {
        G4Tubs* siSolid = new G4Tubs("Si",
                                     0.0,
                                     radius,
                                     0.5 * substrateThickness,
                                     0.0,
                                     2.0 * M_PI);
        fSiLogical = new G4LogicalVolume(siSolid,
                                         fSiMaterial,
                                         "Si");

        // Align top of Si substrate with bottom of Al2O3 layer
        const G4double siCenterZ = 0.5 * thickness + 0.5 * substrateThickness;
        fSiPhysical = new G4PVPlacement(
            nullptr,
            G4ThreeVector(0, 0, siCenterZ),
            fSiLogical,
            "Si",
            fWorldLogical,
            false,
            0,
            true);
    }

    // Define a dedicated region for the Al2O3 target so that we can
    // attach the PAI model only in this thin layer.
    auto* regionStore = G4RegionStore::GetInstance();
    if (fAl2O3Region) {
        regionStore->DeRegister(fAl2O3Region);
        delete fAl2O3Region;
        fAl2O3Region = nullptr;
    }
    if (auto* existing = regionStore->GetRegion("Al2O3Region", false)) {
        regionStore->DeRegister(existing);
        delete existing;
    }
    fAl2O3Region = new G4Region("Al2O3Region");
    fAl2O3Region->AddRootLogicalVolume(fAl2O3Logical);

    G4cout << "\n--- Geometry ---" << G4endl;
    G4cout << "Al2O3 layer:" << G4endl;
    G4cout << "  Thickness: " << fSampleThickness / nm << " nm" << G4endl;
    G4cout << "  Diameter: " << 2.0 * radius / nm << " nm" << G4endl;
    G4cout << "  Volume: " << M_PI * radius * radius * thickness / (nm*nm*nm) << " nm³" << G4endl;
    if (fMaxStep > 0.) {
        G4cout << "  Max step: " << fMaxStep / nm << " nm" << G4endl;
    }
    if (fSiLogical && substrateThickness > 0.) {
        G4cout << "Si substrate:" << G4endl;
        G4cout << "  Thickness: " << substrateThickness / nm << " nm" << G4endl;
        G4cout << "  Diameter: " << 2.0 * radius / nm << " nm" << G4endl;
        G4cout << "  Volume: " << M_PI * radius * radius * substrateThickness / (nm*nm*nm) << " nm³" << G4endl;
    }

    return fWorldPhysical;
}

G4double DetectorConstruction::GetSampleThickness() const
{
    return fSampleThickness;
}

void DetectorConstruction::SetSampleThickness(G4double thickness)
{
    if (thickness > 0.) {
        fSampleThickness = thickness;
    }
}

G4double DetectorConstruction::GetSubstrateThickness() const
{
    return fSubstrateThickness;
}

void DetectorConstruction::SetSubstrateThickness(G4double thickness)
{
    if (thickness < 0.) {
        return;
    }
    fSubstrateThickness = thickness;
}

G4double DetectorConstruction::GetSampleRadius() const
{
    return fSampleRadius;
}

void DetectorConstruction::SetSampleRadius(G4double radius)
{
    if (radius > 0.) {
        fSampleRadius = radius;
    }
}

G4double DetectorConstruction::GetMaxStep() const
{
    return fMaxStep;
}

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
    if (maxStep > 0.) {
        fMaxStep = maxStep;
        if (fAl2O3Logical) {
            fAl2O3Logical->SetUserLimits(new G4UserLimits(fMaxStep));
        }
    }
}
