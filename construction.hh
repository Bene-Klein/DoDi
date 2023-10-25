#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Tet.hh"
#include "G4Polyhedra.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4MultiUnion.hh"
#include "CADMesh.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh" 
#include "G4GenericMessenger.hh"
#include <cmath>

#include "detector.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public: 
    MyDetectorConstruction();
    ~MyDetectorConstruction();

    G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }

    virtual G4VPhysicalVolume *Construct();
private: 
    virtual void ConstructSDandField();

    
    G4Transform3D Rotation(G4double, G4double, G4double, G4double, G4double, G4double, G4double);
    G4Transform3D doubleRotation(G4double, G4double, G4double, G4double, G4double, G4double, G4double);


    G4int nCols, nRows;
    G4double gerscale,thetax,thetay,thetaz;

    G4Polyhedra *solidDoDi, *solidDoDi1;
    G4Tubs *solidDetector,*solidTube;
    G4Box *solidGermanium, *solidWorld;
    G4LogicalVolume *logicWorld,*logicDoDi,*logicDoDi1, *logicGermanium, *logicDetector, *logicObject, *logicmesh, *logicTube;
    G4VPhysicalVolume *physWorld,*physDoDi,*physDoDi1,*physDoDi0,*physDoDi01,*physDetector,*physDetector1,*physDetector2,*physDetector3,*physTube;

    G4Material *worldMat, *germanium, *aluminum;
    G4GenericMessenger *fMessenger;

    void DefineMaterials();

    G4LogicalVolume *fScoringVolume;

    //G4int nRows, nCols;
};
#endif
