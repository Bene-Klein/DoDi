#include "detector.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{}
MySensitiveDetector::~MySensitiveDetector()
{}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
    G4Track *track = aStep->GetTrack();

    track->SetTrackStatus(fStopAndKill);

    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

    G4ThreeVector posPhoton = preStepPoint->GetPosition();

    //G4cout << "positition: "<< posPhoton << G4endl;

    const G4VTouchable *touchable =aStep->GetPreStepPoint()->GetTouchable();

    G4int copyNo = touchable->GetCopyNumber();

    G4int parentId = aStep->GetTrack()->GetParentID();
    
    G4int EventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();

    G4double Time = aStep->GetPreStepPoint()->GetGlobalTime();
    
    //G4cout << "Time: "<< Time <<G4endl;
    //G4cout << "ID: "<< EventID <<G4endl;

    G4VPhysicalVolume *physVol = touchable->GetVolume();
    G4ThreeVector posDetector = physVol->GetTranslation();

    //G4cout << "Detector position: "<< posDetector <<G4endl;

    G4AnalysisManager *man = G4AnalysisManager::Instance();

    man->FillNtupleDColumn(0,0,Time);
    man->FillNtupleIColumn(0,1,copyNo);
    man->FillNtupleIColumn(0,2,EventID); //FillNtupleDColumn(Ntuple Number, Entry Number, Entry)

    man->AddNtupleRow();

}