//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 
/// \file BDDetectorConstruction.cc
/// \brief Implementation of the BDDetectorConstruction class

#include "BDDetectorConstruction.hh"
#include "BDCalorimeterSD.hh"
#include "BDParameterisation.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//______________________________________________________________________________

G4ThreadLocal 
G4GlobalMagFieldMessenger* BDDetectorConstruction::fMagFieldMessenger = 0; 

//______________________________________________________________________________
BDDetectorConstruction::BDDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true),
   fNofLayers(-1),
   fNDiffLayers(0)
{

}
//______________________________________________________________________________
BDDetectorConstruction::~BDDetectorConstruction()
{ 

}  
//______________________________________________________________________________
G4VPhysicalVolume* BDDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}
//______________________________________________________________________________
void BDDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
 
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_Al");
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}
//______________________________________________________________________________
G4VPhysicalVolume* BDDetectorConstruction::DefineVolumes()
{

  // DFlay modifications 
  // - Change absorber from lead to aluminum, NLayers = 3 
  // - Change gap to vacuum and turn off sensitive detector functionality  

  // Geometry parameters
  fNofLayers = 3; // 10; 
  G4double absoThickness = 10.*mm;
  G4double gapThickness =  5.*mm;
  G4double calorSizeXY  = 10.*cm;

  auto layerThickness = absoThickness + gapThickness;
  auto calorThickness = fNofLayers * layerThickness;
  auto worldSizeXY    = 5*m; // 1.2*calorSizeXY;
  auto worldSizeZ     = 5*m; // 1.2*calorThickness; 
  
  // Get materials
  auto defaultMaterial  = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_Al");
  auto gapMaterial      = G4Material::GetMaterial("Galactic");
  
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  // World
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  // Calorimeter
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
                         
  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
   
  // Layer
  auto layerS 
    = new G4Box("Layer",           // its name
                 calorSizeXY/2, calorSizeXY/2, layerThickness/2); //its size
                         
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 fNofLayers,       // number of replica
                 layerThickness);  // witdth of replica
  
  // Absorber
  auto absorberS 
    = new G4Box("Abso",            // its name
                 calorSizeXY/2, calorSizeXY/2, absoThickness/2); // its size
                         
  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV");        // its name

  // visualization 
  G4VisAttributes *visAbso = new G4VisAttributes(); 
  visAbso->SetColor( G4Color::Blue() );
  visAbso->SetForceWireframe(false); 
  absorberLV->SetVisAttributes(visAbso); 
                                   
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,-gapThickness/2.), // its position
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

   // Gap
   auto gapS 
     = new G4Box("Gap",             // its name
                  calorSizeXY/2, calorSizeXY/2, gapThickness/2); // its size
                          
   auto gapLV
     = new G4LogicalVolume(
                  gapS,             // its solid
                  gapMaterial,      // its material
                  "GapLV");         // its name
                                    
   new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., absoThickness/2), // its position
                  gapLV,            // its logical volume                         
                  "Gap",            // its name
                  layerLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps 
  
 
  // A case for diffuser
  // - made of vacuum 
  // - allows placement of the volume in same mother as the calorimeter
  //   (can't have two replicas or parameterised volumes in same mother...)  
  G4double diffCase_xy = 0.5*m; 
  G4double diffCase_z  = 0.5*m; 

  auto diffCaseS 
    = new G4Box("diffCase",           // its name
                 diffCase_xy/2.,diffCase_xy/2.,diffCase_z/2.); // its size
                         
  auto diffCaseLV
    = new G4LogicalVolume(
                 diffCaseS,        // its solid
                 defaultMaterial,  // its material
                 "diffCase");      // its name
                                   
  new G4PVPlacement(
                 0,                        // no rotation
                 G4ThreeVector(0,0,2.*m),  // location in mother volume 
                 diffCaseLV,               // its logical volume                         
                 "diffCase",               // its name
                 worldLV,                  // its mother  volume
                 false,                    // no boolean operation
                 0,                        // copy number
                 fCheckOverlaps);          // checking overlaps 

  BuildDiffuser(diffCaseLV,'A'); 
  // BuildDiffuser_HallA(worldLV);
  // BuildDiffuser_HallC(worldLV);
 
  // print parameters
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << fNofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  // Visualization attributes
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  diffCaseLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);

  // Always return the physical World
  return worldPV;
}
//______________________________________________________________________________
void BDDetectorConstruction::BuildDiffuser(G4LogicalVolume *logicMother,char Hall){
   // parameterised build of the diffuser
   // build first plate (same for Hall A or C)  
   G4double inch     = 25.4*mm;
   G4double r_min    = 2.*inch; 
   G4double r_max    = 5.*inch;
   G4double thk      = 0.125*inch;
   G4double startPhi = 255.*deg;
   G4double dPhi     = 30.*deg; 
  
   // choose the origin of the device (where the first plate starts, relative to the mother volume)  
   G4ThreeVector P0 = G4ThreeVector(0,0,0);
 
   auto plateMaterial = G4Material::GetMaterial("G4_Al");

   if(Hall=='A') fNDiffLayers = 15; 
   if(Hall=='C') fNDiffLayers = 16; 

   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Blue() ); 
 
   G4VSolid *solidPlate     = new G4Tubs("plate",r_min,r_max,thk/2.,startPhi,dPhi);
   G4LogicalVolume *plateLV = new G4LogicalVolume(solidPlate,plateMaterial,"plateLV");
   plateLV->SetVisAttributes(vis);

   G4VPVParameterisation *plateParam = new BDParameterisation(Hall,P0); 
   // G4VPhysicalVolume *physDiffuser   = new G4PVParameterised("Diffuser",plateLV,logicMother,kZAxis,fNDiffLayers,plateParam); 
   new G4PVParameterised("Diffuser",plateLV,logicMother,kZAxis,fNDiffLayers,plateParam); 

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildDiffuser_HallA(G4LogicalVolume *logicMother){
   // Hall A beam dump rocking-vane diffuser
   // - three layers of aluminum (alloy 6601).  each has a different thickness
   // - plates have 5-mm gaps between them  
   // - based on JLab-TN-16-024
   // FIXME: - Is this the right implementation memory wise? Need to attach a *single* 
   //          sensitive detector to the logical volumes.  I think giving them all the same 
   //          logical volume name should accomplish this...
   //        - What are the right sweep angles, and inner/outer radii?  These are 
   //          not specified in the technote.    
   //        - SOLUTION: G4VPVParameterisation class defining the geometry and shape of each object based on copy number 
 
   G4double z0      = 2.*m;      // arbitrary
   G4double gap     = 5.*mm; 
   G4double dz      = 0.;  

   G4double inch    = 25.4*mm; 
   G4double thk_100 = 0.100*inch; 
   G4double thk_125 = 0.125*inch;
   G4double width   = 2.*inch;    // arbitrary  
 
   G4double r_min  = 5.*inch;     // arbitrary
   G4double r_max  = r_min + width;  
 
   // sweep angles 
   G4double startPhi_30  = 255.*deg;
   G4double startPhi_90  = 225.*deg;
   G4double dPhi_30      = 30.*deg;   // step size
   G4double dPhi_60      = 60.*deg;   
   G4double dPhi_90      = 90.*deg;   
   // offset towards beam left 
   G4double startPhi_60l = 255.*deg;
   
   auto plateMaterial = G4Material::GetMaterial("G4_Al");

   // 1/8" plate
   G4Tubs *plate_t125_a30 = new G4Tubs("plate_t125_a30",r_min,r_max,thk_125/2.,startPhi_30,dPhi_30); 

   // 1/10" plate
   G4Tubs *plate_t100_a30 = new G4Tubs("plate_t100_a30",r_min,r_max,thk_100/2.,startPhi_30,dPhi_30); 

   // 1/8" plate, rotated to beam right    
   G4Tubs *plate_t125_a60l = new G4Tubs("plate_t125_a60l",r_min,r_max,thk_125/2.,startPhi_60l,dPhi_60); 

   // 1/10" plate, rotated to beam right    
   G4Tubs *plate_t100_a60l = new G4Tubs("plate_t100_a60l",r_min,r_max,thk_100/2.,startPhi_60l,dPhi_60); 

   // 1/8" plate, 90 deg sweep    
   G4Tubs *plate_t125_a90 = new G4Tubs("plate_t125_a90",r_min,r_max,thk_125/2.,startPhi_90,dPhi_90);
 
   // 1/10" plate, 90 deg sweep    
   G4Tubs *plate_t100_a90 = new G4Tubs("plate_t100_a90",r_min,r_max,thk_100/2.,startPhi_90,dPhi_90); 

   // make an array of logical volumes to be placed 
   const int NP = 15; 
   G4LogicalVolume **logicPlate = new G4LogicalVolume*[NP];
   logicPlate[0]  = new G4LogicalVolume(plate_t125_a30 ,plateMaterial,"logicPlateHallA"); // 1/8" plate, 30 deg wide   
   logicPlate[1]  = new G4LogicalVolume(plate_t125_a30 ,plateMaterial,"logicPlateHallA"); // 1/8" plate, 30 deg wide  
   logicPlate[2]  = new G4LogicalVolume(plate_t100_a30 ,plateMaterial,"logicPlateHallA"); // 1/10" plate, 30 deg wide 
   logicPlate[3]  = new G4LogicalVolume(plate_t125_a60l,plateMaterial,"logicPlateHallA"); // 1/8" plate, 60 deg wide, rotated to beam right
   logicPlate[4]  = new G4LogicalVolume(plate_t125_a60l,plateMaterial,"logicPlateHallA"); // 1/8" plate, 60 deg wide, rotated to beam right 
   logicPlate[5]  = new G4LogicalVolume(plate_t100_a60l,plateMaterial,"logicPlateHallA"); // 1/10" plate, 60 deg wide, rotated to beam right 
   logicPlate[6]  = new G4LogicalVolume(plate_t125_a90 ,plateMaterial,"logicPlateHallA"); // 1/8" plate, 90 deg wide 
   logicPlate[7]  = new G4LogicalVolume(plate_t100_a90 ,plateMaterial,"logicPlateHallA"); // 1/10" plate, 90 deg wide 
   logicPlate[8]  = new G4LogicalVolume(plate_t125_a90 ,plateMaterial,"logicPlateHallA"); // 1/8" plate, 90 deg wide 
   logicPlate[9]  = new G4LogicalVolume(plate_t100_a60l,plateMaterial,"logicPlateHallA"); // 1/10" plate, 60 deg wide, rotated to beam right 
   logicPlate[10] = new G4LogicalVolume(plate_t125_a60l,plateMaterial,"logicPlateHallA"); // 1/8" plate, 60 deg wide, rotated to beam right
   logicPlate[11] = new G4LogicalVolume(plate_t125_a60l,plateMaterial,"logicPlateHallA"); // 1/8" plate, 60 deg wide, rotated to beam right
   logicPlate[12] = new G4LogicalVolume(plate_t100_a30 ,plateMaterial,"logicPlateHallA"); // 1/10" plate, 30 deg wide 
   logicPlate[13] = new G4LogicalVolume(plate_t125_a30 ,plateMaterial,"logicPlateHallA"); // 1/8" plate, 30 deg wide   
   logicPlate[14] = new G4LogicalVolume(plate_t125_a30 ,plateMaterial,"logicPlateHallA"); // 1/8" plate, 30 deg wide   
  
   G4VisAttributes **vis = new G4VisAttributes*[NP];
   for(int i=0;i<NP;i++) vis[i] = new G4VisAttributes();  
   vis[0]->SetColour( G4Colour::Red()     ); 
   vis[1]->SetColour( G4Colour::Red()     ); 
   vis[2]->SetColour( G4Colour::Yellow()  ); 
   vis[3]->SetColour( G4Colour::Red()     ); 
   vis[4]->SetColour( G4Colour::Red()     ); 
   vis[5]->SetColour( G4Colour::Yellow()  ); 
   vis[6]->SetColour( G4Colour::Red()     ); 
   vis[7]->SetColour( G4Colour::Yellow()  ); 
   vis[8]->SetColour( G4Colour::Red()     ); 
   vis[9]->SetColour( G4Colour::Yellow()  ); 
   vis[10]->SetColour( G4Colour::Red()    ); 
   vis[11]->SetColour( G4Colour::Red()    ); 
   vis[12]->SetColour( G4Colour::Yellow() ); 
   vis[13]->SetColour( G4Colour::Red()    ); 
   vis[14]->SetColour( G4Colour::Red()    ); 
 
   G4double THK[NP] = {thk_125,thk_125,thk_100,thk_125,thk_125,
                       thk_100,thk_125,thk_100,thk_125,thk_100,
                       thk_125,thk_125,thk_100,thk_125,thk_125};

   fNDiffLayers = NP; 

   bool isBoolean = false;
   char physName[200]; 
   // char logicName[200];
   // char msg[200];
   
   // place volumes
   for(int i=0;i<NP;i++){
      sprintf(physName ,"physPlateHallA");
      // sprintf(logicName,"logicPlateHallA_%02d",i);
      // logicPlate[i]->SetName(logicName); 
      if(i==0){
	 dz = 0;
      }else{
	 dz = dz + gap + 0.5*(THK[i-1] + THK[i]); 
      }
      // sprintf(msg,"plate %02d, dz = %.3lf mm",i,dz/mm);
      // std::cout << msg << std::endl; 
      logicPlate[i]->SetVisAttributes(vis[i]);  
      new G4PVPlacement(0,                           // rotation relative to mother
	                G4ThreeVector(0,0,z0+dz),    // position relative to mother 
                        logicPlate[i],               // logical volume
                        physName,                    // physical name 
                        logicMother,                 // logic mother
                        isBoolean,                   // boolean solid?
                        i,                           // copy no 
                        fCheckOverlaps);             // check overlaps
   }

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildDiffuser_HallC(G4LogicalVolume *logicMother){
   // Hall C beam dump rocking-vane diffuser
   // - three layers of aluminum (alloy 6601).  each has a different thickness
   // - plates have 5-mm gaps between them  
   // - based on JLab-TN-16-024
   // FIXME: - Is this the right implementation memory wise? Need to attach a *single* 
   //          sensitive detector to the logical volumes.  I think giving them all the same 
   //          logical volume name should accomplish this...
   //        - What are the right sweep angles, and inner/outer radii?  These are 
   //          not specified in the technote.     
   //        - SOLUTION: G4VPVParameterisation class defining the geometry and shape of each object based on copy number 
   
   fNDiffLayers     = 16; 
 
   G4double inch    = 25.4*mm; 
   G4double gap     = 0.195*inch; // spacing  

   G4double z0      = 2.*m;      // arbitrary
   G4double dz      = 0.;  

   // plate thickness 
   G4double thk_031 = 0.03125*inch; 
   G4double thk_062 = 0.06250*inch; 
   G4double thk_094 = 0.09375*inch; 
   G4double thk_125 = 0.12500*inch;
   // plate radii 
   G4double width   = 2.*inch;    // arbitrary  
   G4double r_min   = 5.*inch;    // arbitrary
   G4double r_max   = r_min + width;  
 
   // sweep angles 
   G4double startPhi_30  = 255.*deg;
   G4double startPhi_90  = 225.*deg;
   G4double startPhi_60l = 255.*deg;  
   // step sizes 
   G4double dPhi_30      = 30.*deg;   
   G4double dPhi_60      = 60.*deg;   
   G4double dPhi_90      = 90.*deg;   
   
   auto plateMaterial = G4Material::GetMaterial("G4_Al");

   // 1/8" plate, 30 deg sweep
   G4Tubs *plate_t125_a30  = new G4Tubs("plate_t125_a30",r_min,r_max,thk_125/2.,startPhi_30,dPhi_30); 

   // 1/8" plate, 60 deg sweep, rotated to beam left
   G4Tubs *plate_t125_a60l = new G4Tubs("plate_t125_a60l",r_min,r_max,thk_125/2.,startPhi_60l,dPhi_60); 

   // 1/32" plate, 30 deg sweep
   G4Tubs *plate_t031_a30  = new G4Tubs("plate_t031_a30",r_min,r_max,thk_031/2.,startPhi_30,dPhi_30); 

   // 1/16" plate, 30 deg sweep
   G4Tubs *plate_t062_a30  = new G4Tubs("plate_t062_a30",r_min,r_max,thk_062/2.,startPhi_30,dPhi_30); 

   // 3/32" plate, rotated to beam left    
   G4Tubs *plate_t094_a60l = new G4Tubs("plate_t094_a60l",r_min,r_max,thk_094/2.,startPhi_60l,dPhi_60); 

   // 3/32" plate, 90 deg sweep    
   G4Tubs *plate_t094_a90  = new G4Tubs("plate_t094_a90",r_min,r_max,thk_094/2.,startPhi_90,dPhi_90);


   // make an array of logical volumes to be placed 
   const int NP = 16; 
   G4LogicalVolume **logicPlate = new G4LogicalVolume*[NP];
   logicPlate[0]  = new G4LogicalVolume(plate_t125_a30 ,plateMaterial,"logicPlateHallC"); // 1/8" plate, 30 deg wide   
   logicPlate[1]  = new G4LogicalVolume(plate_t125_a30 ,plateMaterial,"logicPlateHallC"); // 1/8" plate, 30 deg wide 
   logicPlate[2]  = new G4LogicalVolume(plate_t125_a30 ,plateMaterial,"logicPlateHallC"); // 1/8" plate, 30 deg wide
   logicPlate[3]  = new G4LogicalVolume(plate_t125_a30 ,plateMaterial,"logicPlateHallC"); // 1/8" plate, 30 deg wide
   logicPlate[4]  = new G4LogicalVolume(plate_t125_a60l,plateMaterial,"logicPlateHallC"); // 1/8" plate, 60 deg wide, rotated to beam left 
   logicPlate[5]  = new G4LogicalVolume(plate_t125_a60l,plateMaterial,"logicPlateHallC"); // 1/8" plate, 60 deg wide, rotated to beam left  
   logicPlate[6]  = new G4LogicalVolume(plate_t094_a90 ,plateMaterial,"logicPlateHallC"); // 3/32" plate, 90 deg wide 
   logicPlate[7]  = new G4LogicalVolume(plate_t094_a90 ,plateMaterial,"logicPlateHallC"); // 3/32" plate, 90 deg wide 
   logicPlate[8]  = new G4LogicalVolume(plate_t094_a90 ,plateMaterial,"logicPlateHallC"); // 3/32" plate, 90 deg wide 
   logicPlate[9]  = new G4LogicalVolume(plate_t094_a60l,plateMaterial,"logicPlateHallC"); // 3/32" plate, 60 deg wide, rotated to beam left 
   logicPlate[10] = new G4LogicalVolume(plate_t094_a60l,plateMaterial,"logicPlateHallC"); // 3/32" plate, 60 deg wide, rotated to beam left 
   logicPlate[11] = new G4LogicalVolume(plate_t094_a60l,plateMaterial,"logicPlateHallC"); // 3/32" plate, 60 deg wide, rotated to beam left 
   logicPlate[12] = new G4LogicalVolume(plate_t062_a30 ,plateMaterial,"logicPlateHallC"); // 1/16" plate, 30 deg wide
   logicPlate[13] = new G4LogicalVolume(plate_t062_a30 ,plateMaterial,"logicPlateHallC"); // 1/16" plate, 30 deg wide 
   logicPlate[14] = new G4LogicalVolume(plate_t031_a30 ,plateMaterial,"logicPlateHallC"); // 1/32" plate, 30 deg wide  
   logicPlate[15] = new G4LogicalVolume(plate_t031_a30 ,plateMaterial,"logicPlateHallC"); // 1/32" plate, 30 deg wide   
  
   G4VisAttributes **vis = new G4VisAttributes*[NP];
   for(int i=0;i<NP;i++) vis[i] = new G4VisAttributes();  
   vis[0]->SetColour( G4Colour::Blue()    ); 
   vis[1]->SetColour( G4Colour::Blue()    ); 
   vis[2]->SetColour( G4Colour::Blue()    ); 
   vis[3]->SetColour( G4Colour::Blue()    ); 
   vis[4]->SetColour( G4Colour::Blue()    ); 
   vis[5]->SetColour( G4Colour::Blue()    ); 
   vis[6]->SetColour( G4Colour::Green()   ); 
   vis[7]->SetColour( G4Colour::Green()   ); 
   vis[8]->SetColour( G4Colour::Green()   ); 
   vis[9]->SetColour( G4Colour::Green()   ); 
   vis[10]->SetColour( G4Colour::Green()  ); 
   vis[11]->SetColour( G4Colour::Green()  ); 
   vis[12]->SetColour( G4Colour::Yellow() ); 
   vis[13]->SetColour( G4Colour::Yellow() ); 
   vis[14]->SetColour( G4Colour::Red()    ); 
   vis[15]->SetColour( G4Colour::Red()    ); 
 
   G4double THK[NP] = {thk_125,thk_125,thk_125,thk_125,
                       thk_125,thk_125,thk_094,thk_094,
                       thk_094,thk_094,thk_094,thk_094,
                       thk_062,thk_062,thk_031,thk_031};

   fNDiffLayers = NP;  

   bool isBoolean = false;
   char physName[200];
   // char logicName[200];
   // char msg[200];
   
   // place volumes
   for(int i=0;i<NP;i++){
      sprintf(physName ,"physPlateHallC");
      // sprintf(logicName,"logicPlateHallC_%02d",i);
      // logicPlate[i]->SetName(logicName); 
      if(i==0){
	 dz = 0;
      }else{
	 dz = dz + gap + 0.5*(THK[i-1] + THK[i]); 
      }
      // sprintf(msg,"plate %02d, dz = %.3lf mm",i,dz/mm);
      // std::cout << msg << std::endl; 
      logicPlate[i]->SetVisAttributes(vis[i]);  
      new G4PVPlacement(0,                           // rotation relative to mother
	                G4ThreeVector(0,0,z0+dz),    // position relative to mother 
                        logicPlate[i],               // logical volume
                        physName,                    // physical name 
                        logicMother,                 // logic mother
                        isBoolean,                   // boolean solid?
                        i,                           // copy no 
                        fCheckOverlaps);             // check overlaps
   }
   
}
//______________________________________________________________________________
void BDDetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // Sensitive detectors
  auto absoSD     = new BDCalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("AbsoLV",absoSD);

  auto diffuserSD = new BDCalorimeterSD("DiffuserSD","DiffuserHitsCollection",fNDiffLayers); 
  G4SDManager::GetSDMpointer()->AddNewDetector(diffuserSD);
  SetSensitiveDetector("plateLV",diffuserSD); 

  // auto gapSD      = new BDCalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  // G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  // SetSensitiveDetector("GapLV",gapSD);

  // Magnetic field
  // - Create global magnetic field messenger.
  // - Uniform magnetic field is then created automatically if
  //   the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}
//______________________________________________________________________________
