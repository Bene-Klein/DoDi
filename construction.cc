#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
    fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");

    fMessenger->DeclareProperty("nCols", nCols, "Number of columns");
    fMessenger->DeclareProperty("nRows", nRows, "Number of rows");

    fMessenger = new G4GenericMessenger(this, "/germanium/", "Germanium Size");

    fMessenger->DeclareProperty("gerscale", gerscale, "Germanium scale");

    nCols = 100;
    nRows = 100;
    gerscale = 1;
    DefineMaterials();
}

MyDetectorConstruction::~MyDetectorConstruction()
{
}

void MyDetectorConstruction::DefineMaterials()
{
    G4NistManager *nist = G4NistManager::Instance();

    worldMat = nist->FindOrBuildMaterial("G4_AIR");
    germanium = nist->FindOrBuildMaterial("G4_Ge");
    water = nist->FindOrBuildMaterial("G4_WATER");
}
G4Transform3D MyDetectorConstruction::Rotation(G4double theta, G4double x_1, G4double y_1, G4double z_1, G4double x_2, G4double y_2, G4double z_2)
{
    G4Rotate3D Rot = G4Rotate3D(theta, G4ThreeVector(x_1, y_1, z_1));
    G4Translate3D pos = G4Translate3D(x_2, y_2, z_2);
    G4Transform3D Trans = G4Rotate3D(theta, G4ThreeVector(x_1, y_1, z_1)) * G4Translate3D(x_2, y_2, z_2);
    /*
    Rot->rotateX(thetax);
    Rot->rotateY(thetay);
    Rot->rotateZ(thetaz);
    */
    return Trans;
}
G4Transform3D MyDetectorConstruction::doubleRotation(G4double theta, G4double x_1, G4double y_1, G4double z_1, G4double x_2, G4double y_2, G4double z_2)
{
    G4Rotate3D Rot = G4Rotate3D(theta, G4ThreeVector(x_1, y_1, z_1));
    G4Translate3D pos = G4Translate3D(x_2, y_2, z_2);
    G4Transform3D Trans = G4Rotate3D(theta, G4ThreeVector(x_1, y_1, z_1)) * G4Rotate3D(180 * degree, G4ThreeVector(0, 1, 0)) * G4Translate3D(x_2, y_2, z_2);
    /*
    Rot->rotateX(thetax);
    Rot->rotateY(thetay);
    Rot->rotateZ(thetaz);
    */
    return Trans;
}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
    // test//
    G4double alpha = 72 * M_PI / 180;
    G4double beta = 54 * M_PI / 180;
    G4double aWorld = 0.5 * m;                         // mm
    G4double lWorld = aWorld / (2 * sin(36 * degree)); // mm

    G4double r_i = (aWorld / 20) * sqrt(250 + 110 * sqrt(5)); // mm
    G4double r_k = (aWorld / 4) * (3 + sqrt(5));              // mm
    G4double r_e = (aWorld / 4) * sqrt(3) * (1 + sqrt(5));    // mm
    G4double d_l = aWorld / (2 * tan(36 * degree));           // mm

    G4double R_i = (1.113516364) * aWorld; // mm
    G4double R_k = 1.309016994 * aWorld;   // mm

    G4double Beta = acos(r_i / r_k);                                // Rad
    G4double Alpha = acos(r_i / r_e);                               // Rad
    G4double Gamma_div2 = asin(aWorld / (2 * r_e));                 // Rad
    G4double Theta = M_PI - 2 * atan((sqrt(5) - 1) / 2);            // Rad
    G4double Phi = (M_PI / 2) + atan(((sqrt(5) - 1) / 2));          // Rad
    G4double Kappa = 2 * Beta - acos(d_l * cos(36 * degree / r_i)); // Rad

    G4double energy[2] = {1.239841939 * eV / 0.9, 1.239841939 * eV / 0.2};
    G4double rindexWorld[2] = {1.0, 1.0};
    G4double rindexWater[2] = {1.33, 1.33};
    G4double transmittance[2] = {0, 0};
    G4double reflectivity[2] = {0.9, 0.9};
    G4double absH2O[2] = {20 * m, 20 * m};

    G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
    mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);
    worldMat->SetMaterialPropertiesTable(mptWorld);

    G4MaterialPropertiesTable *mptwater = new G4MaterialPropertiesTable();
    mptwater->AddProperty("RINDEX", energy, rindexWater, 2);
    mptwater->AddProperty("ABSLENGTH", energy, absH2O, 2);
    water->SetMaterialPropertiesTable(mptwater);

    G4double world = 1 * m; // mm

    G4double phiStart = 0;
    G4double phiStart1 = 36 * M_PI / 180;
    G4double phiTotal = 2 * M_PI;
    G4int numSide = 5;
    G4int numZPlanes = 2;
    G4double rInner[] = {0, 0};
    G4double rOuter[] = {0, d_l};
    G4double zPlane[] = {0, r_i};

    G4double pRMin = 0;
    G4double pRMax = 10 * mm;
    G4double pRMax1 = 100 * mm;
    G4double pDz = 1 * mm;
    G4double pDz1 = 5 * mm;
    G4double pSPhi = 0;
    G4double pDPhi = 2 * M_PI;

    G4cout << "a: " << aWorld << G4endl;
    G4cout << "l: " << lWorld << G4endl;
    G4cout << "d_l: " << d_l << G4endl;
    G4cout << "r_i/aWorld: " << r_i / aWorld << G4endl;
    G4cout << "r_i: " << r_i << G4endl;
    G4cout << "r_k/aWorld: " << r_k / aWorld << G4endl;
    G4cout << "beta: " << Beta * 180 / M_PI << G4endl;
    G4cout << "theta: " << Theta * 180 / M_PI << G4endl;
    G4cout << "phi: " << Phi * 180 / M_PI << G4endl;
    G4cout << " pi?:" << 2 * (Alpha + Beta + Gamma_div2) << G4endl;

    // G4Box(*name,*size)
    solidWorld = new G4Box("solidWorld", world, world, world);

    solidDoDi = new G4Polyhedra("solidDoDi", phiStart, phiTotal, numSide, numZPlanes, zPlane, rInner, rOuter);

    solidTube = new G4Tubs("solidtube", pRMin, pRMax, pDz, pSPhi, pDPhi);

    solidDetector = new G4Tubs("solidDetector", pRMin, pRMax1, pDz1, pSPhi, pDPhi);

    // G4LogicalVolume(*solidVolume,*material,*name)
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");

    logicDoDi = new G4LogicalVolume(solidDoDi, water, "logicDoDi");

    logicTube = new G4LogicalVolume(solidTube, worldMat, "logicTube");

    logicDetector = new G4LogicalVolume(solidDetector, water, "logicDetector");

    // G4PVPlacement(*Rotation,*Offset in Threevector,*logic Volume,*name,*Mothervolume,*boolean operation, *copynumber,*check for overlaps)
    physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);

    physDoDi1 = new G4PVPlacement(
        Rotation(0, 0, 0, r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 1, true);
    physDoDi2 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta), 0, r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 2, true);
    physDoDi3 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(72 * degree), r_i * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 3, true);
    physDoDi4 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 4, true);
    physDoDi5 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(72 * degree), -r_i * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 5, true);
    physDoDi6 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), -r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 6, true);
    physDoDi7 = new G4PVPlacement(
        doubleRotation(0, 0, 0, r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 7, true);
    physDoDi8 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta), 0, r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 8, true);
    physDoDi9 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(72 * degree), r_i * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 9, true);
    physDoDi10 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 10, true);
    physDoDi11 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(72 * degree), -r_i * tan(Beta) * sin(72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 11, true);
    physDoDi12 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), -r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, 0),
        logicDoDi, "physDoDi", logicWorld, false, 12, true);

    physDetector1 = new G4PVPlacement(
        doubleRotation(0, -0 * r_i * tan(Beta), 0, r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 01, true);
    physDetector2 = new G4PVPlacement(
        doubleRotation(1 * 180 * degree, -1 * r_i * tan(Beta), 0, r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 02, true);
    physDetector3 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta), 0, r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 03, true);
    physDetector4 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(1 * 72 * degree), pow(-1, 1) * r_i * tan(Beta) * sin(1* 72 * degree), r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 04, true);
    physDetector5 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(1 * 72 * degree), pow(-1, 2) * r_i * tan(Beta) * sin(1 * 72 * degree), r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 05, true);
    physDetector6 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), pow(-1, 1) * r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 06, true);
    physDetector7 = new G4PVPlacement(
        Rotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), pow(-1, 2) * r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 07, true);
    physDetector8 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(1 * 72 * degree), pow(-1, 1) * r_i * tan(Beta) * sin(1* 72 * degree), r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 8, true);
    physDetector9 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(1 * 72 * degree), pow(-1, 2) * r_i * tan(Beta) * sin(1 * 72 * degree), r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 9, true);
    physDetector10 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), pow(-1, 1) * r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 10, true);
    physDetector11 = new G4PVPlacement(
        doubleRotation(180 * degree, -r_i * tan(Beta) * cos(2 * 72 * degree), pow(-1, 2) * r_i * tan(Beta) * sin(2 * 72 * degree), r_i, 0, 0, r_i + 5 * mm),
        logicDetector, "physDetector", logicWorld, false, 11, true);

    G4OpticalSurface *surfDoDi = new G4OpticalSurface("WaterSurface1");
    surfDoDi->SetType(dielectric_dielectric);
    surfDoDi->SetModel(unified);
    surfDoDi->SetFinish(groundfrontpainted);

    G4LogicalBorderSurface *waterSurface1 = new G4LogicalBorderSurface("WaterSurface", physDoDi1, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface2 = new G4LogicalBorderSurface("WaterSurface", physDoDi2, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface3 = new G4LogicalBorderSurface("WaterSurface", physDoDi3, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface4 = new G4LogicalBorderSurface("WaterSurface", physDoDi4, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface5 = new G4LogicalBorderSurface("WaterSurface", physDoDi5, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface6 = new G4LogicalBorderSurface("WaterSurface", physDoDi6, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface7 = new G4LogicalBorderSurface("WaterSurface", physDoDi7, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface8 = new G4LogicalBorderSurface("WaterSurface", physDoDi8, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface9 = new G4LogicalBorderSurface("WaterSurface", physDoDi9, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface10 = new G4LogicalBorderSurface("WaterSurface", physDoDi10, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface11 = new G4LogicalBorderSurface("WaterSurface", physDoDi11, physWorld, surfDoDi);
    G4LogicalBorderSurface *waterSurface12 = new G4LogicalBorderSurface("WaterSurface", physDoDi12, physWorld, surfDoDi);

    G4MaterialPropertiesTable *mptDoDi = new G4MaterialPropertiesTable();
    mptDoDi->AddProperty("TRANSMITTANCE", energy, transmittance, 2);
    mptDoDi->AddProperty("REFLECTIVITY", energy, reflectivity, 2);
    surfDoDi->SetMaterialPropertiesTable(mptDoDi);

    return physWorld;
    // return physWorld1;
}
void MyDetectorConstruction::ConstructSDandField()
{
    MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

    logicDetector->SetSensitiveDetector(sensDet);
}
