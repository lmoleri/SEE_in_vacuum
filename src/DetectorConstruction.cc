#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RegionStore.hh"

DetectorConstruction::DetectorConstruction()
    : fWorldMaterial(nullptr),
      fAl2O3Material(nullptr),
      fWorldLogical(nullptr),
      fAl2O3Logical(nullptr),
      fWorldPhysical(nullptr),
      fAl2O3Physical(nullptr),
      fAl2O3Region(nullptr)
{
    fSampleThickness = 20.0 * nm;
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
    
    G4cout << "\n--- Material properties ---" << G4endl;
    G4cout << "World material: " << fWorldMaterial->GetName() << G4endl;
    G4cout << "Al2O3 density: " << fAl2O3Material->GetDensity() / (g/cm3) << " g/cm³" << G4endl;
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
    // Thickness and diameter are configurable (defaults: 20 nm, 100 nm)
    G4double thickness = fSampleThickness > 0. ? fSampleThickness : 20.0 * nm;
    G4double radius = 50.0 * nm;  // radius = diameter/2

    fSampleThickness = thickness;
    
    G4Tubs* al2o3Solid = new G4Tubs("Al2O3",
                                    0.0,           // inner radius
                                    radius,        // outer radius
                                    0.5 * thickness,  // half-length in z
                                    0.0,           // start angle
                                    2.0 * M_PI);   // span angle
    
    fAl2O3Logical = new G4LogicalVolume(al2o3Solid,
                                        fAl2O3Material,
                                        "Al2O3");
    
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
