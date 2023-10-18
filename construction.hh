#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Tet.hh"
#include "G4Polyhedra.hh"
#include "G4RotationMatrix.hh"
#include "G4MultiUnion.hh"
#include "CADMesh.hh"
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
    
    G4RotationMatrix* Rotation(G4double, G4double, G4double);

    G4int nCols, nRows;
    G4double gerscale,thetax,thetay,thetaz;

    G4Polyhedra *solidDoDi;
    G4Box *solidGermanium, *solidDetector,*solidWorld;
    G4LogicalVolume *logicWorld,*logicDoDi, *logicGermanium, *logicDetector, *logicObject, *logicmesh;
    G4VPhysicalVolume *physWorld,*physDoDi,*physDoDi1,*physDoDi2,*physDoDi3, *physGermanium, *physDetector, *physObject, *physmesh;

    G4Material *worldMat, *germanium;
    G4GenericMessenger *fMessenger;

    void DefineMaterials();

    G4LogicalVolume *fScoringVolume;

    //G4int nRows, nCols;
};
#endif
