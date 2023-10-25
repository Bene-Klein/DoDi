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
    aluminum = nist->FindOrBuildMaterial("G4_WATER");

}
G4Transform3D MyDetectorConstruction::Rotation(G4double theta,G4double x_1,G4double y_1,G4double z_1,G4double x_2,G4double y_2,G4double z_2)
{
    G4Rotate3D Rot = G4Rotate3D(theta,G4ThreeVector(x_1,y_1,z_1));
    G4Translate3D pos = G4Translate3D(x_2,y_2,z_2);
    G4Transform3D Trans = G4Rotate3D(theta,G4ThreeVector(x_1,y_1,z_1))*G4Translate3D(x_2,y_2,z_2);
    /*
    Rot->rotateX(thetax);
    Rot->rotateY(thetay);
    Rot->rotateZ(thetaz);
    */
    return Trans;
}
G4Transform3D MyDetectorConstruction::doubleRotation(G4double theta,G4double x_1,G4double y_1,G4double z_1,G4double x_2,G4double y_2,G4double z_2)
{
    G4Rotate3D Rot = G4Rotate3D(theta,G4ThreeVector(x_1,y_1,z_1));
    G4Translate3D pos = G4Translate3D(x_2,y_2,z_2);
    G4Transform3D Trans = G4Rotate3D(theta,G4ThreeVector(x_1,y_1,z_1))*G4Rotate3D(180*degree,G4ThreeVector(0,1,0))*G4Translate3D(x_2,y_2,z_2);
    /*
    Rot->rotateX(thetax);
    Rot->rotateY(thetay);
    Rot->rotateZ(thetaz);
    */
    return Trans;
}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
    //test//
    G4double alpha = 72*M_PI/180;
    G4double beta = 54*M_PI/180;
    G4double aWorld = 0.5*m;//mm
    G4double lWorld = aWorld/(2*sin(36*degree));//mm

    G4double r_i = (aWorld/20)*sqrt(250+110*sqrt(5));//mm
    G4double r_k = (aWorld/4)*(3+sqrt(5));//mm
    G4double r_e = (aWorld/4)*sqrt(3)*(1+sqrt(5)); //mm   
    G4double d_l = aWorld/(2*tan(36*degree));//mm

    G4double R_i = (1.113516364)*aWorld;//mm
    G4double R_k = 1.309016994*aWorld;//mm
    
    G4double Beta = acos(r_i/r_k);//Rad
    G4double Alpha = acos(r_i/r_e);//Rad
    G4double Gamma_div2 = asin(aWorld/(2*r_e));//Rad
    G4double Theta = M_PI-2*atan((sqrt(5)-1)/2);//Rad
    G4double Phi =(M_PI/2)+atan(((sqrt(5)-1)/2));//Rad
    G4double Kappa =2*Beta-acos(d_l*cos(36*degree/r_i));//Rad

 

    G4double world = 1*m;//mm

    G4double phiStart = 0;
    G4double phiStart1 =36*M_PI/180;
    G4double phiTotal = 2*M_PI;
    G4int numSide =5;
    G4int numZPlanes =2;
    G4double rInner[] = {0,0};
    G4double rOuter[] = {0,d_l};
    G4double zPlane[] = {0,r_i}; 

    G4double pRMin =0;
    G4double pRMax =10*mm;
    G4double pRMax1 =100*mm;
    G4double pDz = 1*mm;
    G4double pDz1 = 5*mm;
    G4double pSPhi= 0;
    G4double pDPhi= 2*M_PI;

    G4cout << "a: " << aWorld << G4endl;
    G4cout << "l: " << lWorld << G4endl;
    G4cout << "d_l: " << d_l << G4endl;
    G4cout << "r_i/aWorld: " << r_i/aWorld << G4endl;
    G4cout << "r_i: " << r_i << G4endl;
    G4cout << "r_k/aWorld: " << r_k/aWorld << G4endl;
    G4cout << "beta: " << Beta*180/M_PI << G4endl;
    G4cout << "theta: " << Theta*180/M_PI << G4endl;
    G4cout << "phi: " << Phi*180/M_PI << G4endl;
    G4cout << " pi?:" << 2*(Alpha+Beta+Gamma_div2) << G4endl;

    //G4Box(*name,*size)
    solidWorld =new G4Box("solidWorld", world,world,world);

    solidDoDi =new G4Polyhedra("solidDoDi",phiStart,phiTotal,numSide,numZPlanes,zPlane,rInner,rOuter);
    
    solidTube =new G4Tubs("solidtube",pRMin,pRMax,pDz,pSPhi,pDPhi);

    solidDetector = new G4Tubs("solidDetector",pRMin,pRMax1,pDz1,pSPhi,pDPhi);

    //G4LogicalVolume(*solidVolume,*material,*name)
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");

    logicDoDi = new G4LogicalVolume(solidDoDi, aluminum, "logicDoDi");

    logicTube = new G4LogicalVolume(solidTube, worldMat, "logicTube");

    logicDetector = new G4LogicalVolume(solidDetector, worldMat, "logicDetector");

    //G4PVPlacement(*Rotation,*Offset in Threevector,*logic Volume,*name,*Mothervolume,*boolean operation, *copynumber,*check for overlaps)
    physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),logicWorld,"physWorld",0,false,0,true);
       
    for(G4int i=0; i<=1;i++)
        {
            physDoDi0 = new G4PVPlacement(
                Rotation(i*180*degree,-i*r_i*tan(Beta),0,r_i,0,0,0),
                logicDoDi,"physDoDi0",logicWorld,false,i,true);
        }

    for(G4int j=1; j<=2;j++)
    {
            for(G4int i=1; i<=2;i++)
        {
            physDoDi = new G4PVPlacement(
                Rotation(180*degree,-r_i*tan(Beta)*cos(i*72*degree),pow(-1,j)*r_i*tan(Beta)*sin(i*72*degree),r_i,0,0,0),
                logicDoDi,"physDoDi",logicWorld,false,i+j*10,true);
        }
       
    }

    for(G4int i=0; i<=1;i++)
    {
        physDoDi01 = new G4PVPlacement(
            doubleRotation(i*180*degree,-i*r_i*tan(Beta),0,r_i,0,0,0),
            logicDoDi,"physDoDi01",logicWorld,false,i+2,true);
    }

    for(G4int j=0; j<=1;j++)
    {
            for(G4int i=1; i<=2;i++)
        {
            physDoDi1 = new G4PVPlacement(
                doubleRotation(180*degree,-r_i*tan(Beta)*cos(i*72*degree),pow(-1,j)*r_i*tan(Beta)*sin(i*72*degree),r_i,0,0,0),
                logicDoDi,"physDoDi1",logicWorld,false,i+j*10+100,true);
        }
       
    }
    for(G4int i=0; i<=1;i++)
        {
            physDetector = new G4PVPlacement(
                doubleRotation(i*180*degree,-i*r_i*tan(Beta),0,r_i,0,0,r_i+5*mm),
                logicDetector,"physDetector",logicWorld,false,i,true);
        }

        physDetector1 = new G4PVPlacement(
                Rotation(180*degree,-r_i*tan(Beta),0,r_i,0,0,r_i+5*mm),
                logicDetector,"physDetector",logicWorld,false,3,true);

        for(G4int j=1; j<=2;j++)
    {
            for(G4int i=1; i<=2;i++)
        {
            physDetector2 = new G4PVPlacement(
                Rotation(180*degree,-r_i*tan(Beta)*cos(i*72*degree),pow(-1,j)*r_i*tan(Beta)*sin(i*72*degree),r_i,0,0,r_i+5*mm),
                logicDetector,"physDetector",logicWorld,false,i+j*10+1000,true);
        }
       
    }

        for(G4int j=1; j<=2;j++)
    {
            for(G4int i=1; i<=2;i++)
        {
            physDetector3 = new G4PVPlacement(
                doubleRotation(180*degree,-r_i*tan(Beta)*cos(i*72*degree),pow(-1,j)*r_i*tan(Beta)*sin(i*72*degree),r_i,0,0,r_i+5*mm),
                logicDetector,"physDetector",logicWorld,false,i+j*10+10000,true);
        }
       
    }
    
    /*
    for(G4int i=0; i<=360;i=i+10)
    {
        physDoDi1 = new G4PVPlacement(Rotation(i*degree,-cos(Beta*M_PI/180)*r_i,0,r_i,0,0,0),logicDoDi,"physDoDi1",logicWorld,false,0,true);
    }
    */
    return physWorld;
    //return physWorld1;

}
void MyDetectorConstruction::ConstructSDandField()
{
    MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

    logicDetector->SetSensitiveDetector(sensDet);
}
