#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Region.hh"
#include "globals.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    G4double GetSampleThickness() const;
    void SetSampleThickness(G4double thickness);
    G4double GetSubstrateThickness() const;
    void SetSubstrateThickness(G4double thickness);

private:
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    G4Material* fWorldMaterial;
    G4Material* fAl2O3Material;
    G4Material* fSiMaterial;
    
    G4LogicalVolume* fWorldLogical;
    G4LogicalVolume* fAl2O3Logical;
    G4LogicalVolume* fSiLogical;
    
    G4VPhysicalVolume* fWorldPhysical;
    G4VPhysicalVolume* fAl2O3Physical;
    G4VPhysicalVolume* fSiPhysical;

    // Region used for applying the PAI model
    G4Region* fAl2O3Region = nullptr;

    G4double fSampleThickness = 0.;
    G4double fSubstrateThickness = 0.;
};

#endif
