#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "detector.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction 
{
public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();
    
    virtual G4VPhysicalVolume* Construct();

private:
    G4LogicalVolume *logicDetector;
    virtual void ConstructSDandField();
    

};

#endif // CONSTRUCTION_HH
