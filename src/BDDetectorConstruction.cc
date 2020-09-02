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
#include "BeamDiffuserSD.hh"
#include "BDParameterisation.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
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
  // fNofLayers = 3; // 10; 
  // G4double absoThickness = 10.*mm;
  // G4double gapThickness =  5.*mm;
  // G4double calorSizeXY  = 10.*cm;
  // auto layerThickness = absoThickness + gapThickness;
  // auto calorThickness = fNofLayers * layerThickness;

  auto worldSizeXY    = 5*m; // 1.2*calorSizeXY;
  auto worldSizeZ     = 5*m; // 1.2*calorThickness; 
   
  // Get materials
  auto defaultMaterial  = G4Material::GetMaterial("Galactic");
  // auto absorberMaterial = G4Material::GetMaterial("G4_Al");
  // auto gapMaterial      = G4Material::GetMaterial("Galactic");
   
  // if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
  //   G4ExceptionDescription msg;
  //   msg << "Cannot retrieve materials already defined."; 
  //   G4Exception("B4DetectorConstruction::DefineVolumes()",
  //     "MyCode0001", FatalException, msg);
  // }  
    
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

  // // Calorimeter
  // auto calorimeterS
  //   = new G4Box("Calorimeter",     // its name
  //                calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
  //                        
  // auto calorLV
  //   = new G4LogicalVolume(
  //                calorimeterS,     // its solid
  //                defaultMaterial,  // its material
  //                "Calorimeter");   // its name
  //                                  
  // new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(),  // at (0,0,0)
  //                calorLV,          // its logical volume                         
  //                "Calorimeter",    // its name
  //                worldLV,          // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps);  // checking overlaps
  //  
  // // Layer
  // auto layerS 
  //   = new G4Box("Layer",           // its name
  //                calorSizeXY/2, calorSizeXY/2, layerThickness/2); //its size
  //                        
  // auto layerLV
  //   = new G4LogicalVolume(
  //                layerS,           // its solid
  //                defaultMaterial,  // its material
  //                "Layer");         // its name

  // new G4PVReplica(
  //                "Layer",          // its name
  //                layerLV,          // its logical volume
  //                calorLV,          // its mother
  //                kZAxis,           // axis of replication
  //                fNofLayers,       // number of replica
  //                layerThickness);  // witdth of replica
  // 
  // // Absorber
  // auto absorberS 
  //   = new G4Box("Abso",            // its name
  //                calorSizeXY/2, calorSizeXY/2, absoThickness/2); // its size
  //                        
  // auto absorberLV
  //   = new G4LogicalVolume(
  //                absorberS,        // its solid
  //                absorberMaterial, // its material
  //                "AbsoLV");        // its name

  // // visualization 
  // G4VisAttributes *visAbso = new G4VisAttributes(); 
  // visAbso->SetColor( G4Color::Blue() );
  // visAbso->SetForceWireframe(false); 
  // absorberLV->SetVisAttributes(visAbso); 
  //                                  
  //  new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0.,-gapThickness/2.), // its position
  //                absorberLV,       // its logical volume                         
  //                "Abso",           // its name
  //                layerLV,          // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps);  // checking overlaps 

  //  // Gap
  //  auto gapS 
  //    = new G4Box("Gap",             // its name
  //                 calorSizeXY/2, calorSizeXY/2, gapThickness/2); // its size
  //                         
  //  auto gapLV
  //    = new G4LogicalVolume(
  //                 gapS,             // its solid
  //                 gapMaterial,      // its material
  //                 "GapLV");         // its name
  //                                   
  //  new G4PVPlacement(
  //                 0,                // no rotation
  //                 G4ThreeVector(0., 0., absoThickness/2), // its position
  //                 gapLV,            // its logical volume                         
  //                 "Gap",            // its name
  //                 layerLV,          // its mother  volume
  //                 false,            // no boolean operation
  //                 0,                // copy number
  //                 fCheckOverlaps);  // checking overlaps 
  
  BuildDiffuser(worldLV,'A'); 
  // BuildDiffuser_HallA(worldLV);
  // BuildDiffuser_HallC(worldLV);
 
  // print parameters
  // G4cout
  //   << G4endl 
  //   << "------------------------------------------------------------" << G4endl
  //   << "---> The calorimeter is " << fNofLayers << " layers of: [ "
  //   << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
  //   << " + "
  //   << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
  //   << "------------------------------------------------------------" << G4endl;
  // 

  // Visualization attributes
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  BuildCell(worldLV); 
  BuildBeamPipe(worldLV); 

  // auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  // simpleBoxVisAtt->SetVisibility(true);
  // calorLV->SetVisAttributes(simpleBoxVisAtt);

  // Always return the physical World
  return worldPV;
}
//______________________________________________________________________________
void BDDetectorConstruction::BuildBeamPipe(G4LogicalVolume *logicMother){
   // a mock beam pipe
   G4double r1_min   = 5.*cm;  
   G4double r1_max   = 10.*cm;  
   G4double len1     = 50.*cm;
   G4double r2_min   = 10.*cm;  
   G4double r2_max   = 20.*cm;  
   G4double len2     = 100.*cm;
   G4double startPhi = 0*deg;  
   G4double dPhi     = 360*deg;  
 
   // define materials 
   G4Material *Al  = G4Material::GetMaterial("G4_Al");

   G4Tubs *tube1 = new G4Tubs("tube1",r1_min,r1_max,len1/2.,startPhi,dPhi);   
   G4Tubs *tube2 = new G4Tubs("tube2",r2_min,r2_max,len2/2.,startPhi,dPhi);    

   G4ThreeVector P = G4ThreeVector(0,0,len2/2.); 
   G4UnionSolid *tubeAssembly = new G4UnionSolid("beamPipe",tube1,tube2,0,P);

   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Green() ); 
 
   G4LogicalVolume *pipeLV = new G4LogicalVolume(tubeAssembly,Al,"pipeLV");
   pipeLV->SetVisAttributes(vis);  

   // placement 
   new G4PVPlacement(0,
	             G4ThreeVector(0,0,1*m),
                     pipeLV,
                     "physPipe",
                     logicMother,
                     false,
                     0,
                     fCheckOverlaps);


}
//______________________________________________________________________________
void BDDetectorConstruction::BuildCell(G4LogicalVolume *logicMother){
   // a simple tube of aluminum, Cu, or 3He 
   G4double r_min    = 0.*cm;  
   G4double r_max    = 2.*cm;  
   G4double length   = 50.*cm;
   G4double startPhi = 0*deg;  
   G4double dPhi     = 360*deg;  
  
   // define materials 
   G4Material *Al  = G4Material::GetMaterial("G4_Al");
   G4Material *Cu  = G4Material::GetMaterial("G4_Cu");

   G4int ncomponents=0;
   G4double Z=0,N=0,abundance=0; 
   G4Isotope *iso_3He = new G4Isotope( "He3", Z=2, N=3 );
   G4Element *el3He   = new G4Element("Helium3", "3He", ncomponents=1 ); //Define isotopically pure Helium-3 
   el3He->AddIsotope( iso_3He, abundance=100.0*perCent );

   G4double rho       = 10.77*atmosphere*(3.016*g/Avogadro)/(300*kelvin*k_Boltzmann);
   G4Material *pol3He = new G4Material("pol3He",rho,1);
   pol3He->AddElement(el3He, 1);

   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::White() ); 

   // define solid 
   G4Tubs *solidTube = new G4Tubs("tube",r_min,r_max,length/2.,startPhi,dPhi);

   // define logical volume 
   G4LogicalVolume *tubeLV = new G4LogicalVolume(solidTube,pol3He,"tubeLV");
   tubeLV->SetVisAttributes(vis); 

   // placement 
   new G4PVPlacement(0,
	             G4ThreeVector(0,0,0.5*m),
                     tubeLV,
                     "physTube",
                     logicMother,
                     false,
                     0,
                     fCheckOverlaps);

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildDiffuser(G4LogicalVolume *logicMother,char Hall){
 
   // A case for diffuser
   // - made of vacuum 
   // - allows placement of the volume in same mother as the calorimeter
   //   (can't have two replicas or parameterised volumes in same mother...)  
 
   auto defaultMaterial  = G4Material::GetMaterial("Galactic");

   G4double inch        = 25.4*mm;
   G4double diffCase_x  = 12.*inch;  
   G4double diffCase_y  = 6.*inch;  // 0.5*m; 
   G4double diffCase_z  = 15.*cm;   // 0.5*m; 

   auto diffCaseS 
      = new G4Box("diffCase",           // its name
	    diffCase_x/2.,diffCase_y/2.,diffCase_z/2.); // its size

   auto diffCaseLV
      = new G4LogicalVolume(
	    diffCaseS,        // its solid
	    defaultMaterial,  // its material
	    "diffCase");      // its name
  
   // where to place the diffuser 
   // note: the (x,y) center of the diffuser plates is centered on this logical volume 
   G4double xd = 0.;
   G4double yd = 0.; 
   G4double zd = 2.5*m; 
   G4ThreeVector P_case = G4ThreeVector(xd,yd,zd); 

   new G4PVPlacement(0,                        // no rotation
	             P_case,                   // location in mother volume 
	             diffCaseLV,               // its logical volume                         
	             "diffCase",               // its name
	             logicMother,              // its mother  volume
	             false,                    // no boolean operation
	             0,                        // copy number
	             fCheckOverlaps);          // checking overlaps 
  
   G4VisAttributes *visCase = new G4VisAttributes();
   visCase->SetForceWireframe(); 
 
   // diffCaseLV->SetVisAttributes(G4VisAttributes::GetInvisible());
   diffCaseLV->SetVisAttributes(visCase); 

   // parameterised build of the diffuser
   // build first plate (same for Hall A or C)  
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
   new G4PVParameterised("Diffuser",plateLV,diffCaseLV,kZAxis,fNDiffLayers,plateParam); 

   // auto diffuserSD = new BeamDiffuserSD("DiffuserSD","DiffuserHitsCollection"); 
   // G4SDManager::GetSDMpointer()->AddNewDetector(diffuserSD);
   // plateLV->SetSensitiveDetector(diffuserSD); 

}
//______________________________________________________________________________
void BDDetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // Sensitive detectors
  // auto absoSD     = new BDCalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  // G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  // SetSensitiveDetector("AbsoLV",absoSD);

  auto diffuserSD = new BeamDiffuserSD("DiffuserSD","DiffuserHitsCollection"); 
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
