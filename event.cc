#include "event.hh"


MyEventAction::MyEventAction(MyRunAction*)
{
    fEdep = 0.;
    /*
    Time =0.;
    CopyNO =0;
    */
    EventID=0;
    
}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
    fEdep = 0.;
    //EventID = GetEventID();   
}


void MyEventAction::EndOfEventAction(const G4Event*)
{
    //G4cout << "Event ID: " << EventID << G4endl;

    
}
