#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
    fParticleGun = new G4ParticleGun(1); 

    

    
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
    delete fParticleGun;
}
/*
void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("proton");

    G4ThreeVector pos(0.,0.,0);
    G4ThreeVector mom(0.,0.,1.);
    
    
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);
    fParticleGun->SetParticleMomentum(100.*GeV);
    fParticleGun->SetParticleDefinition(particle);
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
		
   }
*/
/*

void MyPrimaryGenerator::GeneratePrimaries(G4Event* anEvent) {
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("proton");

    
    if (!particle) {
        G4cerr << "Error: Particle 'electron' not found in the particle table." << G4endl;
        return;  // Exit the function if the particle is not found to avoid further errors.
    }

    G4ThreeVector pos(0., 0., -2.*cm);
    G4ThreeVector momDirection(0., 0., 1.);

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(momDirection);
    fParticleGun->SetParticleEnergy(10.*MeV);
    fParticleGun->SetParticleDefinition(particle);

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
*/
void MyPrimaryGenerator::GeneratePrimaries(G4Event* anEvent) {
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("e-");

    
    if (!particle) {
        G4cerr << "Error: Particle 'electron' not found in the particle table." << G4endl;
        return;  // Exit the function if the particle is not found to avoid further errors.
    }

    G4ThreeVector pos(0., 0., -2.*cm);
    G4ThreeVector momDirection(0., 0., 1.);

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(momDirection);
    fParticleGun->SetParticleEnergy(2.*MeV);
    fParticleGun->SetParticleDefinition(particle);

    fParticleGun->GeneratePrimaryVertex(anEvent);
}

