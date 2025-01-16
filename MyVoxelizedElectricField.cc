#include "MyVoxelizedElectricField.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include <cmath>

MyVoxelizedElectricField::MyVoxelizedElectricField()
: G4ElectroMagneticField(),
  fXmin(0.), fXmax(0.),
  fYmin(0.), fYmax(0.),
  fZmin(0.), fZmax(0.),
  fNx(0), fNy(0), fNz(0),
  fFieldArray(nullptr)
{}

MyVoxelizedElectricField::~MyVoxelizedElectricField() {}

void MyVoxelizedElectricField::SetBoundingBox(G4double xmin, G4double xmax,
                                              G4double ymin, G4double ymax,
                                              G4double zmin, G4double zmax)
{
  fXmin = xmin; fXmax = xmax;
  fYmin = ymin; fYmax = ymax;
  fZmin = zmin; fZmax = zmax;
}

void MyVoxelizedElectricField::SetFieldArray(const G4ThreeVector* fieldArray,
                                             G4int Nx, G4int Ny, G4int Nz)
{
  fFieldArray = fieldArray;
  fNx = Nx; 
  fNy = Ny; 
  fNz = Nz;
}

// Convert (x,y,z) to (i,j,k)
void MyVoxelizedElectricField::GetIndices(G4double x, G4double y, G4double z,
                                          G4int& i, G4int& j, G4int& k) const
{
  // Compute voxel size:
  G4double dx = (fXmax - fXmin) / fNx;
  G4double dy = (fYmax - fYmin) / fNy;
  G4double dz = (fZmax - fZmin) / fNz;

  // For safety, clamp x,y,z inside bounding box
  if(x < fXmin) x = fXmin;
  if(x > fXmax) x = fXmax - 1e-9*dx; // small offset
  if(y < fYmin) y = fYmin;
  if(y > fYmax) y = fYmax - 1e-9*dy;
  if(z < fZmin) z = fZmin;
  if(z > fZmax) z = fZmax - 1e-9*dz;

  i = (G4int) ((x - fXmin) / dx);
  j = (G4int) ((y - fYmin) / dy);
  k = (G4int) ((z - fZmin) / dz);

  // Make sure indices are in [0, Nx-1], etc.
  if(i < 0) i=0; if(i >= fNx) i=fNx-1;
  if(j < 0) j=0; if(j >= fNy) j=fNy-1;
  if(k < 0) k=0; if(k >= fNz) k=fNz-1;
}

// The main method: return E and B at (point[0], point[1], point[2]), ignoring point[3] = time
void MyVoxelizedElectricField::GetFieldValue(const G4double point[4], G4double* field ) const
{
  G4double x = point[0];
  G4double y = point[1];
  G4double z = point[2];

  // 1. find voxel indices
  G4int i, j, k;
  GetIndices(x, y, z, i, j, k);

  // 2. Flatten index = k + Nz*(j + Ny*i)
  G4int idx = k + fNz*( j + fNy*i );
  G4ThreeVector E = fFieldArray[idx]; // E vector in that voxel

  // Set the output array
  field[0] = E.x();  // E_x
  field[1] = E.y();  // E_y
  field[2] = E.z();  // E_z

  // no magnetic field
  field[3] = 0.0;
  field[4] = 0.0;
  field[5] = 0.0;
}
