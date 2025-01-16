#include "detector.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
    // Constructor body, empty in this case
}

MySensitiveDetector::~MySensitiveDetector()
{
    // Destructor body, empty in this case
}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)


{


    G4cout << "Entering ProcessHits function." << G4endl;
    G4Track *track = aStep->GetTrack();

    // Immediately stops the track
    track->SetTrackStatus(fStopAndKill);

    // Gets pre-step and post-step points
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint(); // Not used in the snippet

    // Retrieves and prints the photon's position at the pre-step point
    G4ThreeVector posPhoton = preStepPoint->GetPosition();
    G4cout << "Photon position: " << posPhoton << G4endl;

    // Accesses touchable and its copy number, then prints it
    const G4VTouchable *touchable = preStepPoint->GetTouchable();
    G4int copyNo = touchable->GetCopyNumber();
    G4cout << "Copy number: " << copyNo << G4endl;

    // Retrieves and prints the detector's position
    G4VPhysicalVolume *physVol = touchable->GetVolume();
    G4ThreeVector posDetector = physVol->GetTranslation();
    G4cout << "Detector position: " << posDetector << G4endl;

    return true;
}

