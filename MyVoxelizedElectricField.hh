#ifndef MyVoxelizedElectricField_hh
#define MyVoxelizedElectricField_hh 

#include "G4ElectroMagneticField.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

class MyVoxelizedElectricField : public G4ElectroMagneticField
{
public:
  MyVoxelizedElectricField();
  virtual ~MyVoxelizedElectricField();

  // The key method: returns E and B at a given (x,y,z,t)
  virtual void GetFieldValue(const G4double point[4], G4double* field ) const override;

  // Set bounding box, dimension
  void SetBoundingBox(G4double xmin, G4double xmax,
                      G4double ymin, G4double ymax,
                      G4double zmin, G4double zmax);

  // Set the field array (or pass references, or load from file, etc.)
  // For simplicity, let's expose pointers to arrays:
  void SetFieldArray(const G4ThreeVector* fieldArray, 
                     G4int Nx, G4int Ny, G4int Nz);

private:
  // bounding box
  G4double fXmin, fXmax;
  G4double fYmin, fYmax;
  G4double fZmin, fZmax;

  // voxel counts
  G4int fNx, fNy, fNz;

  // pointer to array of G4ThreeVector
  // We'll keep it flat, i.e. fieldArray[k + Nz*(j + Ny*i)]
  const G4ThreeVector* fFieldArray; 

  // helper function to convert (x,y,z) -> (i,j,k)
  void GetIndices(G4double x, G4double y, G4double z,
                  G4int& i, G4int& j, G4int& k) const;
};

#endif
