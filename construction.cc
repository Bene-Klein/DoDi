#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
    fMessenger = new G4GenericMessenger(this,"/detector/", "Detector Construction");
    

    fMessenger->DeclareProperty("nCols", nCols, "Number of columns");
    fMessenger->DeclareProperty("nRows", nRows, "Number of rows");
    
    fMessenger = new G4GenericMessenger(this,"/germanium/", "Germanium Size");

    fMessenger->DeclareProperty("gerscale", gerscale, "Germanium scale");

    nCols=100;
    nRows=100;
    gerscale=1;
    DefineMaterials();
}

MyDetectorConstruction::~MyDetectorConstruction()
{}

void MyDetectorConstruction::DefineMaterials()
{    
    G4NistManager *nist = G4NistManager::Instance();
    
    worldMat = nist->FindOrBuildMaterial("G4_AIR");
    germanium = nist->FindOrBuildMaterial("G4_Ge");
}

G4RotationMatrix* MyDetectorConstruction::Rotation(G4double thetax,G4double thetay,G4double thetaz)
{
    G4RotationMatrix* Rot =new G4RotationMatrix;
    Rot->rotateX(thetax*deg);
    Rot->rotateY(thetay*deg);
    Rot->rotateZ(thetaz*deg);

    return Rot;
}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
    //test//
    G4double xGermanium = 0.05*m;
    G4double yGermanium = 0.05*m;
    G4double zGermanium = 0.05*m;
    G4double alpha = 72*M_PI/180;
    G4double beta = 54*M_PI/180;
    G4double dWorld = 0.5*m;
    G4double lWorld = dWorld/(2*cos(beta));

    G4double r_i = (dWorld/20)*sqrt(250+110*sqrt(5));
    G4double r_k = (dWorld/4)*(3+sqrt(5));
    G4double Beta = cos(r_i/r_k)*180/M_PI;

    G4double phiStart = 0;
    G4double phiStart1 =36*M_PI/180;
    G4double phiTotal = 2*M_PI;
    G4int numSide =5;
    G4int numZPlanes =2;
    G4double rInner[] = {0,0};
    G4double rOuter[] = {0,lWorld};
    G4double zPlane[] = {0,1.114*dWorld}; 

    G4cout << "a: " << alpha << G4endl;
    G4cout << "b: " << beta << G4endl;
    G4cout << "d: " << dWorld << G4endl;
    G4cout << "l: " << lWorld << G4endl;
    G4cout << "beta: " << Beta << G4endl;
    G4cout << "pRot: " << Rotation(0,0,36) << G4endl;
    //G4Box(*name,*size)
    solidWorld =new G4Box("solidWorld", 1*m,1*m,1*m);

    solidDoDi =new G4Polyhedra("solidDoDi",phiStart,phiTotal,numSide,numZPlanes,zPlane,rInner,rOuter);
    /*Multi union structre
    G4RotationMatrix rotm = G4RotationMatrix();
    G4ThreeVector pos = G4ThreeVector(0,0,0);
    G4Transform3D tr1 = G4Transform3D(rotm,pos);
    
    G4MultiUnion* munion_solid = new G4MultiUnion("DoDi");

    munion_solid->AddNode(*solidDoDi,tr1);
    munion_solid->AddNode(*solidDoDi,tr2);

    munion_solid->Voxelize();*/
    //auto mesh = CADMesh::TessellatedMesh::FromSTL("dodecahedron.stl");
    //auto solidMesh = mesh->GetSolid();
    //G4LogicalVolume(*solidVolume,*material,*name)
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");

    logicDoDi = new G4LogicalVolume(solidDoDi, worldMat, "logicDoDi");
    //logicmesh = new G4LogicalVolume(solidMesh, worldMat, "logicmesh");
  
    //logicObject = new G4LogicalVolume(munion_solid,worldMat,"logicObject");

    //G4PVPlacement(*Rotation,*Offset in Threevector,*logic Volume,*name,*Mothervolume,*boolean operation, *copynumber,*check for overlaps)
    physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),logicWorld,"physWorld",0,false,0,true);
    //physmesh = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),logicmesh,"physmesh",logicWorld,false,0,true);

    //physObject = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),logicObject,"physObject",logicWorld,false,0,true);

    physDoDi = new G4PVPlacement(Rotation(0,0,0), G4ThreeVector(0.,0.,0.),logicDoDi,"physDoDi",logicWorld,false,0,true);
    physDoDi1 = new G4PVPlacement(Rotation(0,2*Beta,36), G4ThreeVector(0.,0.,0.),logicDoDi,"physDoDi1",logicWorld,false,0,true);
    physDoDi2 = new G4PVPlacement(Rotation(Beta,Beta,18), G4ThreeVector(0.,0.,0.),logicDoDi,"physDoDi2",logicWorld,false,0,true);
    //physDoDi3 = new G4PVPlacement(Rotation(0,4*cos(1.114/1.309)*180/M_PI,0), G4ThreeVector(0.,0.,0.),logicDoDi,"physDoDi3",logicWorld,false,0,true);

    /*
    for(G4int i=0; i<=6;i++)
    {
        physDoDi = new G4PVPlacement(Rotation(,,36), G4ThreeVector(0.,0.,0.),logicDoDi,"physDoDi",logicWorld,false,0,true);
    }
   */ 
    return physWorld;
    //return physWorld1;

}
void MyDetectorConstruction::ConstructSDandField()
{
    //MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

    //logicDetector->SetSensitiveDetector(sensDet);
}
