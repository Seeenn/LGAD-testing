#include "construction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VSensitiveDetector.hh"
#include "G4PVPlacement.hh"

MyDetectorConstruction::MyDetectorConstruction() {
    
}

MyDetectorConstruction::~MyDetectorConstruction() {
 
}

G4VPhysicalVolume* MyDetectorConstruction::Construct() {

    
    // ------------------------------------------------------------------------
  // 1. Define Elements
  // ------------------------------------------------------------------------
  G4double a;  // atomic mass
  G4double z;  // atomic number

  // NIST manager for built-in materials (Al, Si, etc.)
  G4NistManager* nist = G4NistManager::Instance();

  //1. Elements 
  G4Element* elSi = new G4Element("Silicon",    "Si", z=14., a=28.0855*g/mole);
  G4Element* elAl = new G4Element("Aluminum",   "Al", z=13., a=26.9815*g/mole);
  G4Element* elB  = new G4Element("Boron",      "B",  z=5.,  a=10.81*g/mole);
  G4Element* elP  = new G4Element("Phosphorus", "P",  z=15., a=30.974*g/mole);
  G4Element* elC  = new G4Element("Carbon",     "C",  z=6.,  a=12.011*g/mole);

  // ------------------------------------------------------------------------
  // 2. Define Materials
  // ------------------------------------------------------------------------

 
  G4double densityAl = 2.70*g/cm3;
  G4Material* Al = new G4Material("Aluminum", densityAl, 1, kStateSolid);
  Al->AddElement(elAl, 1.0);

  
  G4double densitySi = 2.33*g/cm3;
  G4Material* Si = new G4Material("Silicon", densitySi, 1, kStateSolid);
  Si->AddElement(elSi, 1.0);

  G4double fractionP_npp = 2.0e-4;
  G4Material* nppSi = new G4Material("nPlusPlusSi", densitySi, 2, kStateSolid);
  nppSi->AddElement(elSi, 1.0 - fractionP_npp);
  nppSi->AddElement(elP,  fractionP_npp);

  G4double fractionB_pplus = 2.0e-6;
  G4Material* pPlusSi = new G4Material("pPlusSi", densitySi, 2, kStateSolid);
  pPlusSi->AddElement(elSi, 1.0 - fractionB_pplus);
  pPlusSi->AddElement(elB,  fractionB_pplus);

  
  G4double fractionB_sub  = 2.0e-9;
  G4double fractionC_sub  = 2.0e-9;
  G4double fractionSi_sub = 1.0 - fractionB_sub - fractionC_sub;  // ~0.999999996
  G4Material* pSubSi = new G4Material("pSubSi", densitySi, 3, kStateSolid);
  pSubSi->AddElement(elSi, fractionSi_sub);
  pSubSi->AddElement(elB,  fractionB_sub);
  pSubSi->AddElement(elC,  fractionC_sub);

  // ------------------------------------------------------------------------
  // 3. Build Geometry
  // ------------------------------------------------------------------------
  
  G4double topAlThickness       = 0.001*mm;  // ~1 µm
  G4double nppThickness         = 0.001*mm;  // ~1 µm
  G4double pPlusThickness       = 0.002*mm;  // ~2 µm
  G4double pSubstrateThickness  = 0.300*mm;  // ~300 µm
  G4double backAlThickness      = 0.001*mm;  // ~1 µm

  G4double totalThickness = topAlThickness
                          + nppThickness
                          + pPlusThickness
                          + pSubstrateThickness
                          + backAlThickness;  // ~0.305 mm

  // Define world size
  G4double worldSizeXY = 10.0*mm;
  G4double worldSizeZ  = 10.0*mm;

  // Create solid for the World
  G4Box* solidWorld = new G4Box("World",
                                0.5*worldSizeXY,
                                0.5*worldSizeXY,
                                0.5*worldSizeZ);

  G4Material* airMat = nist->FindOrBuildMaterial("G4_AIR"); // or vacuum

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,
                        airMat,
                        "WorldLV");


  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0.,0.,0.),
                      logicWorld,
                      "WorldPV",
                      nullptr,
                      false,
                      0,
                      true);


  G4double zPosition = -0.5*totalThickness;

  // ------------------------------------------------------------------------
  // 4. Place each LGAD layer (with colors)
  // ------------------------------------------------------------------------

  {
    G4String layerName = "TopAl";
    G4Box* solidTopAl =
      new G4Box(layerName, 0.5*worldSizeXY, 0.5*worldSizeXY, 0.5*topAlThickness);

    G4LogicalVolume* logicTopAl =
      new G4LogicalVolume(solidTopAl, Al, layerName+"LV");

   
    G4VisAttributes* topAlVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.5));
    logicTopAl->SetVisAttributes(topAlVis);

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0., zPosition + 0.5*topAlThickness),
                      logicTopAl,
                      layerName+"PV",
                      logicWorld,
                      false,
                      0,
                      true);

    zPosition += topAlThickness;
  }

  
  {
    G4String layerName = "nPlusPlus";
    G4Box* solidNpp =
      new G4Box(layerName, 0.5*worldSizeXY, 0.5*worldSizeXY, 0.5*nppThickness);

    G4LogicalVolume* logicNpp =
      new G4LogicalVolume(solidNpp, nppSi, layerName+"LV");

   
    G4VisAttributes* nppVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.4));
    logicNpp->SetVisAttributes(nppVis);

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0., zPosition + 0.5*nppThickness),
                      logicNpp,
                      layerName+"PV",
                      logicWorld,
                      false,
                      0,
                      true);

    zPosition += nppThickness;
  }

  {
    G4String layerName = "pPlus";
    G4Box* solidPplus =
      new G4Box(layerName, 0.5*worldSizeXY, 0.5*worldSizeXY, 0.5*pPlusThickness);

    G4LogicalVolume* logicPplus =
      new G4LogicalVolume(solidPplus, pPlusSi, layerName+"LV");

    G4VisAttributes* pplusVis = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.4));
    logicPplus->SetVisAttributes(pplusVis);

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0., zPosition + 0.5*pPlusThickness),
                      logicPplus,
                      layerName+"PV",
                      logicWorld,
                      false,
                      0,
                      true);

    zPosition += pPlusThickness;
  }


  {
    G4String layerName = "pSubstrate";
    G4Box* solidSub =
      new G4Box(layerName, 0.5*worldSizeXY, 0.5*worldSizeXY, 0.5*pSubstrateThickness);

    G4LogicalVolume* logicSub =
      new G4LogicalVolume(solidSub, pSubSi, layerName+"LV");

    G4VisAttributes* substrateVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.4));
    logicSub->SetVisAttributes(substrateVis);

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0., zPosition + 0.5*pSubstrateThickness),
                      logicSub,
                      layerName+"PV",
                      logicWorld,
                      false,
                      0,
                      true);

    zPosition += pSubstrateThickness;
  }

  
  {
    G4String layerName = "BackAl";
    G4Box* solidBackAl =
      new G4Box(layerName, 0.5*worldSizeXY, 0.5*worldSizeXY, 0.5*backAlThickness);

    G4LogicalVolume* logicBackAl =
      new G4LogicalVolume(solidBackAl, Al, layerName+"LV");

  
    G4VisAttributes* backAlVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.5));
    logicBackAl->SetVisAttributes(backAlVis);

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0., zPosition + 0.5*backAlThickness),
                      logicBackAl,
                      layerName+"PV",
                      logicWorld,
                      false,
                      0,
                      true);

    zPosition += backAlThickness;
  }

  ConstructSDandField();
 
  return physWorld;


    
}




void MyDetectorConstruction::ConstructSDandField(){

    MySensitiveDetector* sensDet = new MySensitiveDetector("SensitiveDetector");
    if (!logicDetector) {
    G4cerr << "Error: logicDetector is null!" << G4endl;
    return;
}
    // Set the sensitive detector
    logicDetector->SetSensitiveDetector(sensDet);
    

    // (B) Add the electric field
    // 1. Define bounding box for the electric field (matching the detector size)
    G4double xmin = -5.0 * mm, xmax = 5.0 * mm;
    G4double ymin = -5.0 * mm, ymax = 5.0 * mm;
    G4double zmin = -0.5 * mm, zmax = 0.5 * mm;

    // 2. Create an instance of MyVoxelizedElectricField
    MyVoxelizedElectricField* electricField = new MyVoxelizedElectricField();
    electricField->SetBoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);

    // 3. Define voxel grid and electric field array
    const G4int Nx = 10, Ny = 10, Nz = 10;
    std::vector<G4ThreeVector> fieldArray(Nx * Ny * Nz);

    // Initialize the electric field array (e.g., a uniform field in the z-direction)
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                int idx = k + Nz * (j + Ny * i);
                fieldArray[idx] = G4ThreeVector(0.0, 0.0, 1.0 * kilovolt / cm); // Field in z-direction
            }
        }
    }

    // 4. Set the field array and voxel grid dimensions
    electricField->SetFieldArray(fieldArray.data(), Nx, Ny, Nz);

    // 5. Connect the electric field to the Geant4 field manager
    G4FieldManager* fieldManager = new G4FieldManager(electricField);

    // 6. Create an equation of motion for the electric field
    G4EqMagElectricField* equation = new G4EqMagElectricField(electricField);

    // 7. Create a Runge-Kutta stepper for field integration
    G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation);

    // 8. Create a chord finder for trajectory calculation
    G4double minStep = 0.01 * mm; // Minimum step size for field integration
    G4ChordFinder* chordFinder = new G4ChordFinder(electricField, minStep, stepper);

    // 9. Assign the field manager and chord finder
    fieldManager->SetChordFinder(chordFinder);
    logicDetector->SetFieldManager(fieldManager, true);

    G4cout << "Electric field initialized and connected to the detector!" << G4endl;





}
   

