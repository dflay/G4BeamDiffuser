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
#include "IonChamberSD.hh"
#include "BeamDiffuserSD.hh"
#include "BDParameterisation.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
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

  auto worldSizeXY    = 50*m; // 1.2*calorSizeXY;
  auto worldSizeZ     = 50*m; // 1.2*calorThickness; 
   
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

  BuildHe3(worldLV); 
  BuildGlassCell(worldLV);

  BuildIonChamber(worldLV); 
 
  // BuildBeamPipe(worldLV); 

  // DrawAxis(worldLV,'x'); 
  // DrawAxis(worldLV,'y'); 
  DrawAxis(worldLV,'z'); 

  // G4double z_bd = 2.5*m; 
  // BuildBeamDump(worldLV,z_bd);

  // BuildCollimators(worldLV); 
 
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

  // BuildCell(worldLV); 
  // BuildBeamPipe(worldLV); 

  // auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  // simpleBoxVisAtt->SetVisibility(true);
  // calorLV->SetVisAttributes(simpleBoxVisAtt);

  // Always return the physical World
  return worldPV;
}
//______________________________________________________________________________
void BDDetectorConstruction::DrawAxis(G4LogicalVolume *logicMother,char axis,G4double z){

   // draw an axis centered on the origin 

   int validAxis=0;
   if(axis=='x') validAxis = 1; 
   if(axis=='y') validAxis = 1; 
   if(axis=='z') validAxis = 1;
   if(validAxis!=1){
      std::cout << "[BDDetectorConstruction::DrawAxis]: WARNING: Invalid axis = " << axis << std::endl;
      return;
   } 
   
   char axisName[200],axisName_LV[200],axisName_PHY[200];
   sprintf(axisName    ,"%c-axis",axis    ); 
   sprintf(axisName_LV ,"%s_LV"  ,axisName); 
   sprintf(axisName_PHY,"%s_PHY" ,axisName); 

   G4double r_min    = 0.*mm; 
   G4double r_max    = 2.*mm;
   G4double len      = 5*m; 
   G4double startPhi = 0*deg;
   G4double dPhi     = 360*deg;
   G4Tubs *solidTube = new G4Tubs(axisName,r_min,r_max,len/2.,startPhi,dPhi);

   G4VisAttributes *vis = new G4VisAttributes(); 
   if(axis=='x') vis->SetColour( G4Colour::Red()   ); 
   if(axis=='y') vis->SetColour( G4Colour::Green() ); 
   if(axis=='z') vis->SetColour( G4Colour::Blue()  );

   // define materials and logical volume 
   G4Material *Vacuum = G4Material::GetMaterial("Galactic");

   G4LogicalVolume *tube_LV = new G4LogicalVolume(solidTube,Vacuum,axisName_LV); 
   tube_LV->SetVisAttributes(vis); 

   // rotation and placement 
   G4double RX=0,RY=0,RZ=0;
   G4RotationMatrix *rm = new G4RotationMatrix(); 
   if(axis=='x'){
      RX = 0.*deg; RY = 90.*deg; RZ = 0.*deg; 
   }else if(axis=='y'){
      RX = 90.*deg; RY = 0.*deg; RZ = 0.*deg; 
   }else if(axis=='z'){
      RX = 0.*deg; RY = 0.*deg; RZ = 0.*deg; 
   }
   rm->rotateX(RX); rm->rotateY(RY); rm->rotateZ(RZ); 

   new G4PVPlacement(rm,                        // rotation
                     G4ThreeVector(0,0,z),      // position 
                     tube_LV,                   // logical volume   
                     axisName_PHY,              // physical name 
                     logicMother,               // logical mother
                     false,                     // boolean? 
                     0,                         // copy no 
                     fCheckOverlaps);           // check overlaps

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildGlassCell(G4LogicalVolume *logicMother){
   // a glass cell 

   // materials
   G4NistManager *man = G4NistManager::Instance();

   G4Element *elAl = man->FindOrBuildElement("Al");
   G4Element *elSi = man->FindOrBuildElement("Si");
   G4Element *elCa = man->FindOrBuildElement("Ca");
   G4Element *elSr = man->FindOrBuildElement("Sr");
   G4Element *elBa = man->FindOrBuildElement("Ba");
   G4Element *elO  = man->FindOrBuildElement("O" );

   // Cell Glass - GE180 Aluminosilicate Glass
   G4double bigden = 1E9*g/cm3; 
   // Changed names of materials in this composition since the G4Materials used here are only used to make up GE180:
   // SiO2 60.3%
   G4Material* SiO2 = new G4Material("GE180_SiO2", 2.2*g/cm3, 2 );
   SiO2->AddElement(elSi, 1);
   SiO2->AddElement(elO, 2);
   // BaO  18.2%
   G4Material* BaO = new G4Material("GE180_BaO", bigden, 2 );
   BaO->AddElement(elBa, 1);
   BaO->AddElement(elO, 1);
   // Al2O3 14.3%
   G4Material* Al2O3 = new G4Material("GE180_Al2O3", bigden, 2 );
   Al2O3->AddElement(elAl, 2);
   Al2O3->AddElement(elO, 3);
   // CaO   6.5%
   G4Material* CaO = new G4Material("GE180_CaO", bigden, 2 );
   CaO->AddElement(elCa, 1);
   CaO->AddElement(elO, 1);
   // SrO   0.25%
   G4Material* SrO = new G4Material("GE180_SrO", bigden, 2 );
   SrO->AddElement(elSr, 1);
   SrO->AddElement(elO, 1);

   // Density 2.76 g/cm^3
   // Index of Refraction 1.536
   G4Material* GE180 = new G4Material("GE180", 2.76*g/cm3, 5);
   GE180->AddMaterial(SiO2 ,0.6039);
   GE180->AddMaterial(BaO  ,0.1829);
   GE180->AddMaterial(Al2O3,0.1439);
   GE180->AddMaterial(CaO  ,0.0659);
   GE180->AddMaterial(SrO  ,0.0034);

   // build the tube 
   G4double wall     = 1.0*mm;
   G4double OD       = 2.1*cm;  
   G4double r_max    = 0.5*OD; 
   G4double r_min    = r_max - wall; 
   G4double length   = 60.*cm;  
   G4double startPhi = 0;
   G4double dPhi     = 360.*deg;
   G4Tubs *solidTube = new G4Tubs("solidTube_GE180",r_min,r_max,length/2.,startPhi,dPhi); 

   // visualization 
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::White() ); 
   vis->SetForceWireframe(); 

   G4LogicalVolume *tube_LV = new G4LogicalVolume(solidTube,GE180,"GE180_LV");
   tube_LV->SetVisAttributes(vis);  
  
   // placement 
   new G4PVPlacement(0,
	             G4ThreeVector(0,0,0),
                     tube_LV,
                     "GE180_PHY",
                     logicMother,
                     false,
                     0,
                     fCheckOverlaps);

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
void BDDetectorConstruction::BuildHe3(G4LogicalVolume *logicMother){
   // a simple tube of 3He
   G4double OD       = 19.*mm;  
   G4double r_min    = 0.*cm;  
   G4double r_max    = OD/2.;  
   G4double length   = 60.*cm;
   G4double startPhi = 0*deg;  
   G4double dPhi     = 360*deg;  
   G4Tubs *solidTube = new G4Tubs("tube",r_min,r_max,length/2.,startPhi,dPhi);
  
   // define materials 
   // G4Material *Al  = G4Material::GetMaterial("G4_Al");
   // G4Material *Cu  = G4Material::GetMaterial("G4_Cu");

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
   vis->SetColour( G4Colour::Yellow() ); 

   // define logical volume 
   G4LogicalVolume *tubeLV = new G4LogicalVolume(solidTube,pol3He,"tubeLV");
   tubeLV->SetVisAttributes(vis); 

   // placement 
   new G4PVPlacement(0,
	             G4ThreeVector(0,0,0),
                     tubeLV,
                     "He3_PHY",
                     logicMother,
                     false,
                     0,
                     fCheckOverlaps);

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildIonChamber(G4LogicalVolume *logicMother){
   // a dummy ion chamber that's a sensitive detector to see 
   // if we can notice events scattering from glass wall, and not endcaps 

   // material
   auto man = G4NistManager::Instance();

   // Argon
   // G4Element* elAr = man->FindOrBuildElement("Ar");
   // G4double density_Ar = 1.7823*mg/cm3;
   // G4Material *Argon = new G4Material("Argon",density_Ar,1);
   // Argon->AddElement(elAr,1);

   // Nitrogen.  Take this from my notes on ion chambers. 
   // This density is from a typical LHC device 
   G4Element* elN  = man->FindOrBuildElement("N");
   G4double gasden = 1.08*atmosphere*(14.0067*2*g/Avogadro)/(300*kelvin*k_Boltzmann);
   G4Material *icN2 = new G4Material("icN2",gasden,1);
   icN2->AddElement(elN,1);

   // solid 
   G4double vol      = 0.11E-3*m3;    // 0.11 litres 
   G4double length   = 10.*cm;         
   G4double r_max    = sqrt( vol/(pi*length) );  
   G4double r_min    = 0;
   G4double startPhi = 0.*deg;  
   G4double dPhi     = 360.*deg;  
   G4Tubs *solidIC = new G4Tubs("solidIC",r_min,r_max,length/2.,startPhi,dPhi);   
 
   // G4double x_len = 10.*cm; 
   // G4double y_len = 10.*cm; 
   // G4double z_len = 10.*cm;
   // G4Box *solidIC = new G4Box("solidIC",x_len/2.,y_len/2.,z_len/2.); 

   // logical volume
   G4LogicalVolume *ic_LV = new G4LogicalVolume(solidIC,icN2,"ic_LV");

   // placement (arbitrary)  
   G4ThreeVector P = G4ThreeVector(-20.*cm,0.*cm,30.*cm);
   G4RotationMatrix *rm = new G4RotationMatrix();
   rm->rotateX(0); rm->rotateY(10*deg); rm->rotateZ(0); 

   new G4PVPlacement(rm,                    // rotation
	             P,                     // position      
                     ic_LV,                 // logical volume     
                     "ionChamber_PHY",      // name      
                     logicMother,           // mother logical volume      
                     false,                 // is it a boolean solid?    
                     0,                     // copy number    
                     fCheckOverlaps);       // check for overlaps       

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildBeamDump(G4LogicalVolume *logicMother,G4double z0){
   // build the Hall A beam dump 
   // z0 = global offset to all parts   
   // location of beam diffuser front face relative to target pivot (from Ron Lassiter Sept 2020)
   // - target pivot to upstream pipe conical flange:                         1052.4797 inches 
   // - upstream pipe conical flange to ISO wall weldment:                    207.1108 inches
   // - upstream ISO wall weldment to upstream face of diffuser:              24.56 inches
   // - upstream face of beam diffuser to upstream face of downstream flange: 17.3943 inches
   G4double inch  = 25.4*mm; 
   G4double z_us  = z0; // + 1052.4797*inch;
   G4double z_iso = z_us  + 207.1108*inch;
   G4double z_bd  = z_iso + 24.56*inch;
   G4double z_ds  = z_bd  + 17.3943*inch;
   BuildBeamDump_UpstreamPipe(logicMother,z_us);  
   BuildBeamDump_ISOWallWeldment(logicMother,z_iso);
   BuildBeamDump_Diffuser(logicMother,'A',z_bd);
   BuildBeamDump_DownstreamPipe(logicMother,z_ds);  
}
//______________________________________________________________________________
void BDDetectorConstruction::BuildBeamDump_Diffuser(G4LogicalVolume *logicMother,char Hall,G4double z0){
 
   // A case for diffuser
   // - made of vacuum 
   // - allows placement of the volume in same mother as the calorimeter
   //   (can't have two replicas or parameterised volumes in same mother...)  
 
   auto defaultMaterial  = G4Material::GetMaterial("Galactic");

   G4double inch = 25.4*mm;

   // G4double diffCase_x = 12.*inch;  
   // G4double diffCase_y = 6.*inch;  // 0.5*m; 
   // G4double diffCase_z = 15.*cm;   // 0.5*m; 
   // auto diffCaseS = new G4Box("diffCase",diffCase_x/2.,diffCase_y/2.,diffCase_z/2.); 

   G4double BDLength=0;  
   if(Hall=='A') BDLength = 11.44*cm; 
   if(Hall=='C') BDLength = 11.0*cm; // FIXME 

   G4double r_min    = 0.;
   G4double r_max    = 25.*inch;
   G4double len      = BDLength; 
   G4double startPhi = 0.*deg; 
   G4double dPhi     = 360.*deg; 
   auto diffCaseS    = new G4Tubs("diffCase",r_min,r_max,len/2.,startPhi,dPhi);
   auto diffCaseLV   = new G4LogicalVolume(diffCaseS,defaultMaterial,"diffCase"); // its name
   
   // place the diffuser 
   // note: the (x,y) center of the diffuser plates is centered on this logical volume 
   // z0 = location of FRONT FACE of the beam diffuser
   // zz = location of CENTER of the beam diffuser case (coincides with the BD)  
   G4double zz = z0 + BDLength/2.;
   G4ThreeVector P_case = G4ThreeVector(0,0,zz); 

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
   r_min        = 17.67*inch;  
   r_max        = 24.*inch; 
   dPhi         = 38.*deg; 
   startPhi     = 270.*deg - dPhi/2.; 
   G4double thk = 0.125*inch;
  
   // choose the origin of the device (where the first plate starts, relative to the mother volume)
   // origin of BD is the origin of the case  
   zz = -BDLength/2.;  
   G4ThreeVector P0 = G4ThreeVector(0,0,zz);
 
   auto plateMaterial = G4Material::GetMaterial("G4_Al");

   if(Hall=='A') fNDiffLayers = 15; 
   if(Hall=='C') fNDiffLayers = 16; 

   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Blue() ); 
 
   G4VSolid *solidPlate     = new G4Tubs("plate",r_min,r_max,thk/2.,startPhi,dPhi);
   G4LogicalVolume *plateLV = new G4LogicalVolume(solidPlate,plateMaterial,"plateLV");
   plateLV->SetVisAttributes(vis);

   G4VPVParameterisation *plateParam = new BDParameterisation(Hall,P0); 
   new G4PVParameterised("Diffuser",plateLV,diffCaseLV,kZAxis,fNDiffLayers,plateParam); 

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildBeamDump_ISOWallWeldment(G4LogicalVolume *logicMother,G4double z0){
   // Hall A Beam Dump: ISO Wall Weldment 
   // Drawing: JL0015694, JL0015725 

   G4double inch     = 25.4*mm;
   G4double startPhi = 0.*deg;
   G4double dPhi     = 360.*deg;

   // vacuum tube hub [drawing JL0016212] 
   // - outer component  
   G4double r_min_vth   = 0.5*16.*inch;
   G4double r_max_vth   = 0.5*20.*inch;
   G4double len_vth     = 2.75*inch;
   G4Tubs *solidVTH_cyl = new G4Tubs("solidVTH_cyl",r_min_vth,r_max_vth,len_vth/2.,startPhi,dPhi);
   // - cut component   
   G4double r_min_vth_cc = 0.5*16.5*inch;
   G4double r_max_vth_cc = 0.5*25.0*inch;
   G4double len_vth_cc   = 2.50*inch;
   G4Tubs *solidVTH_cc = new G4Tubs("solidVTH_cc",r_min_vth_cc,r_max_vth_cc,len_vth_cc/2.,startPhi,dPhi);
   // subtraction 
   G4ThreeVector Pcc = G4ThreeVector(0,0,len_vth_cc-len_vth);
   G4SubtractionSolid *solidVTH = new G4SubtractionSolid("solidVTH",solidVTH_cyl,solidVTH_cc,0,Pcc);

   // wall [drawings JL0015694, JL0015725]  
   G4double x_len = 77.0*inch;   // from JL0015694  
   G4double y_len = 117.0*inch;  // from JL0015694
   G4double z_len = 0.25*inch;   // from JL0015725

   auto material  = G4Material::GetMaterial("G4_Al"); // might not be aluminum... 

   // solid box 
   G4Box *solidBox = new G4Box("isoBox",x_len/2.,y_len/2.,z_len/2.);  

   // cut a circular hole 
   G4double r_min    = 0.*mm;
   G4double r_max    = 0.5*16.*inch;
   G4double len      = z_len + 5.*inch; // make sure it cuts through  
   G4Tubs *solidTube = new G4Tubs("isoTube",r_min,r_max,len/2.,startPhi,dPhi); 

   // cut the hole in the box
   // - relative to bottom of object, vertical distance is 74.48 inches in drawing JL0015725, 
   //   where the plate is 114.25 inches tall.  Note that this drawing is the inner portion of JL0015694.
   //   However, the height of the whole wall is actually 117 inches (JL0015694), 
   //   so we add 1.375 inches to get 75.86 inches from the bottom of the wall. 
   // - this means we have a distance from the center of the wall: 
   G4double yc     = y_len/2. - (y_len - 75.86*inch);
   G4ThreeVector P = G4ThreeVector(0.*inch,yc,z_len/2.);
   G4SubtractionSolid *solidWall = new G4SubtractionSolid("solidWall",solidBox,solidTube,0,P);

   // now union the solidVTH and solidWall
   G4double zw = len_vth/2. + z_len/2.;
   G4ThreeVector Pw = G4ThreeVector(0,-yc,zw);
   G4UnionSolid *solidISO = new G4UnionSolid("solidISO",solidVTH,solidWall,0,Pw);
 
   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Red() );
   vis->SetForceWireframe(); 

   // logical volume
   G4LogicalVolume *isoWallLV = new G4LogicalVolume(solidISO,material,"isoWallLV");
   isoWallLV->SetVisAttributes(vis); 

   // placement
   G4double z = z0 - 0.5*len_vth; // move center of weldment forward so front face of wall is at z0
   G4ThreeVector P_wall = G4ThreeVector(0,0,z);
   new G4PVPlacement(0,                        // no rotation
	             P_wall,                   // location in mother volume 
	             isoWallLV,                // its logical volume                         
	             "isoWallPHY",             // its name
	             logicMother,              // its mother  volume
	             true,                     // boolean operation? 
	             0,                        // copy number
	             fCheckOverlaps);          // checking overlaps    

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildBeamDump_UpstreamPipe(G4LogicalVolume *logicMother,G4double z0){
   // Hall A Beam Dump: Pipe upstream of ISO Weldment
   // Drawing: JL0009934-C-VAC SPOOL REGION UPPER LEVEL
 
   // materials 
   auto man = G4NistManager::Instance();
  
   G4Element *elAl = man->FindOrBuildElement("Al");
   G4Element *elFe = man->FindOrBuildElement("Fe");
   G4Element *elSi = man->FindOrBuildElement("Si");
   G4Element *elCu = man->FindOrBuildElement("Cu");
   G4Element *elMn = man->FindOrBuildElement("Mn");
   G4Element *elMg = man->FindOrBuildElement("Mg");
   G4Element *elCr = man->FindOrBuildElement("Cr");
   G4Element *elZn = man->FindOrBuildElement("Zn");
   G4Element *elTi = man->FindOrBuildElement("Ti");

   // Aluminum alloy 5052 [https://unitedaluminum.com/united-aluminum-alloy-5052/] 
   G4Material *Aluminum_5052 = new G4Material("Aluminum_5052",2.68*g/cm3,8);
   Aluminum_5052->AddElement(elAl,0.9635);
   Aluminum_5052->AddElement(elSi,0.0025);
   Aluminum_5052->AddElement(elFe,0.0040);
   Aluminum_5052->AddElement(elCu,0.0010);
   Aluminum_5052->AddElement(elMn,0.0010);
   Aluminum_5052->AddElement(elMg,0.0250);
   Aluminum_5052->AddElement(elCr,0.0020);
   Aluminum_5052->AddElement(elZn,0.0010);
   
   // Aluminum alloy 6061 [https://unitedaluminum.com/united-aluminum-alloy-6061/] 
   G4Material *Aluminum_6061 = new G4Material("Aluminum_6061",2.70*g/cm3,9);
   Aluminum_6061->AddElement(elAl,0.9668);
   Aluminum_6061->AddElement(elSi,0.0060);
   Aluminum_6061->AddElement(elFe,0.0070);
   Aluminum_6061->AddElement(elCu,0.0028);
   Aluminum_6061->AddElement(elMn,0.0015);
   Aluminum_6061->AddElement(elMg,0.0100);
   Aluminum_6061->AddElement(elCr,0.0019);
   Aluminum_6061->AddElement(elZn,0.0025);
   Aluminum_6061->AddElement(elTi,0.0015);
 
   G4double inch         = 25.4*mm; 
   G4double TOTAL_LENGTH = 196.91*inch;
   G4double startPhi     = 0.*deg; 
   G4double dPhi         = 360.*deg; 

   // large conical tube [item 1]
   G4double delta      = 0.005*inch;  // FIXME: arbitrary!   
   G4double r_min1_lgc = 0.5*23.54*inch; 
   G4double r_max1_lgc = r_min1_lgc + delta;   
   G4double r_min2_lgc = 0.5*37.54*inch; 
   G4double r_max2_lgc = r_min2_lgc + delta;  
   G4double len_lgc    =  9.95*inch; 
   G4Cons *solidConeLG = new G4Cons("solidConeLG",r_min1_lgc,r_max1_lgc,r_min2_lgc,r_max2_lgc,len_lgc/2.,startPhi,dPhi); 

   // large conical tube [item 1, vacuum] 
   G4Cons *solidConeLG_vac = new G4Cons("solidConeLG_vac",0,r_min1_lgc,0,r_min2_lgc,len_lgc/2.,startPhi,dPhi);  
 
   // small conical tube [item 3] 
   G4double r_min1_smc = 0.5*11.99*inch; 
   G4double r_max1_smc = r_min1_smc + delta;   
   G4double r_min2_smc = 0.5*23.49*inch; 
   G4double r_max2_smc = r_min2_smc + delta;  
   G4double len_smc    =  6.13*inch; 
   G4Cons *solidConeSM = new G4Cons("solidConeSM",r_min1_smc,r_max1_smc,r_min2_smc,r_max2_smc,len_smc/2.,startPhi,dPhi); 
   
   // small conical tube [item 3, vacuum] 
   G4Cons *solidConeSM_vac = new G4Cons("solidConeSM_vac",0,r_min1_smc,0,r_min2_smc,len_smc/2.,startPhi,dPhi); 

   // vacuum window tube [item 4] 
   G4double r_min_vwt = 0.5*12.12*inch; 
   G4double r_max_vwt = 0.5*12.50*inch; 
   G4double len_vwt   = 1.88*inch;
   G4Tubs *solidVWTube = new G4Tubs("solidVWTube",r_min_vwt,r_max_vwt,len_vwt/2.,startPhi,dPhi);
   
   // vacuum window tube [item 4, vacuum] 
   G4Tubs *solidVWTube_vac = new G4Tubs("solidVWTube_vac",0,r_min_vwt,len_vwt/2.,startPhi,dPhi);

   // main tube [item 6] 
   G4double delta2   = 0.4775*inch;  // fudge factor to get the length to match TOTAL_LENGTH 
   G4double r_min_m6 = 0.5*23.90*inch;   // FIXME: This is arbitrary! 
   G4double r_max_m6 = 0.5*24.00*inch;   // FIXME: This is arbitrary! old value = 12.12*inch;   
   G4double len_m6   = 39.85*inch - delta2;    // ESTIMATE: total length of item 6 is 85.36*inch; split into two parts 
   G4Tubs *solidTubeM6 = new G4Tubs("solidTubeM6",r_min_m6,r_max_m6,len_m6/2.,startPhi,dPhi);

   // main tube [item 6, vacuum] 
   G4Tubs *solidTubeM6_vac = new G4Tubs("solidTubeM6_vac",0,r_min_m6,len_m6/2.,startPhi,dPhi);

   // main tube [item 6b, ESTIMATE]  
   G4double r_min_m6b = 0.5*23.90*inch;   // FIXME: This is arbitrary! 
   G4double r_max_m6b = 0.5*24.00*inch;   // FIXME: This is arbitrary! old value = 12.12*inch;   
   G4double len_m6b   = 45.51*inch - delta2;    // derived number.   
   G4Tubs *solidTubeM6b = new G4Tubs("solidTubeM6b",r_min_m6b,r_max_m6b,len_m6b/2.,startPhi,dPhi);
   
   // main tube [item 6b, vacuum] 
   G4Tubs *solidTubeM6b_vac = new G4Tubs("solidTubeM6b_vac",0,r_min_m6b,len_m6b/2.,startPhi,dPhi);

   // main tube [item 2]  
   G4double r_min_m2 = 0.5*23.90*inch;   // FIXME: This is arbitrary! 
   G4double r_max_m2 = 0.5*24.00*inch;   // FIXME: This is arbitrary! old value = 12.12*inch;   
   G4double len_m2   = 52.69*inch - delta2;    // ESTIMATE: total length of item 2 is 91.38*inch; split into two parts 
   G4Tubs *solidTubeM2 = new G4Tubs("solidTubeM2",r_min_m2,r_max_m2,len_m2/2.,startPhi,dPhi);
   
   // main tube [item 2, vacuum] 
   G4Tubs *solidTubeM2_vac = new G4Tubs("solidTubeM2_vac",0,r_min_m2,len_m2/2.,startPhi,dPhi);

   // main tube [item 2b, ESTIMATE]  
   G4double r_min_m2b = 0.5*23.90*inch;   // FIXME: This is arbitrary! 
   G4double r_max_m2b = 0.5*24.00*inch;   // FIXME: This is arbitrary! old value = 12.12*inch;   
   G4double len_m2b   = 38.69*inch - delta2;       // derived number.   
   G4Tubs *solidTubeM2b = new G4Tubs("solidTubeM2b",r_min_m2b,r_max_m2b,len_m2b/2.,startPhi,dPhi);

   // main tube [item 2b, vacuum] 
   G4Tubs *solidTubeM2b_vac = new G4Tubs("solidTubeM2b_vac",0,r_min_m2b,len_m2b/2.,startPhi,dPhi);

   // large flange [item 13, drawing JL0012786] 
   G4double r_min_lgf  = 0.5*37.50*inch; 
   G4double r_max_lgf  = 0.5*46.00*inch;
   G4double len_lgf    = 0.500*inch; 
   G4Tubs *solidLGF = new G4Tubs("solidLGF",r_min_lgf,r_max_lgf,len_lgf/2.,startPhi,dPhi);
  
   // large flange [item 13, vacuum]  
   G4Tubs *solidLGF_vac = new G4Tubs("solidLGF_vac",0,r_min_lgf,len_lgf/2.,startPhi,dPhi);

   // flange with o-ring [item 9, drawing JL0029536] 
   G4double r_min_for = 0.5*23.25*inch; 
   G4double r_max_for = 0.5*29.53*inch; 
   G4double len_for   = 1.00*inch; 
   G4Tubs *solidFORing = new G4Tubs("solidFORing",r_min_for,r_max_for,len_for/2.,startPhi,dPhi);
     
   // flange with o-ring [item 9, vacuum] 
   G4Tubs *solidFORing_vac = new G4Tubs("solidFORing_vac",0,r_min_for,len_for/2.,startPhi,dPhi);

   // aperture plate flange [item 14, drawing JL0058855]
   // simplified approach: single material of aluminum 
   G4double r_min_apf = 0.5*3.50*inch; 
   G4double r_max_apf = 0.5*30.43*inch; 
   G4double len_apf   = 1.50*inch; 
   G4Tubs *solidAPF = new G4Tubs("solidAPF",r_min_apf,r_max_apf,len_apf/2.,startPhi,dPhi);
   // // - outer plate, Aluminum 6061 
   // r_min = 0.5*20.00*inch; 
   // r_max = 0.5*30.43*inch; 
   // len   = 1.50*inch; 
   // G4Tubs *solidAPFo = new G4Tubs("solidAPFo",r_min,r_max,len/2.,startPhi,dPhi);
   // // - inner plate, Aluminum 5052 
   // r_min = 0.5*3.50*inch; 
   // r_max = 0.5*20.00*inch; 
   // len   = 1.50*inch; 
   // G4Tubs *solidAPFi = new G4Tubs("solidAPFi",r_min,r_max,len/2.,startPhi,dPhi);

   // vacuum insert into aperture plate flange 
   G4Tubs *solidAPF_vac = new G4Tubs("solidAPF_vac",0,r_min_apf,len_apf/2.,startPhi,dPhi); 
   
   // vacuum tube support ring [item 12] 
   G4double r_min_vtsr = 0.5*24.01*inch; 
   G4double r_max_vtsr = 0.5*27.01*inch;
   G4double len_vtsr   = 0.25*inch;  
   G4Tubs *solidVTSRing = new G4Tubs("solidVTSRing",r_min_vtsr,r_max_vtsr,len_vtsr/2.,startPhi,dPhi);
   
   // vacuum tube support ring [item 12, vacuum] 
   G4Tubs *solidVTSRing_vac = new G4Tubs("solidVTSRing_vac",0,r_min_vtsr,len_vtsr/2.,startPhi,dPhi);

   // vacuum step-down flange [item 5, drawing JL0009940] 
   G4double r_min_vsdf = 0.5*12.25*inch;
   G4double r_max_vsdf = 0.5*16.00*inch; 
   G4double len_vsdf   = 0.620*inch;  
   G4Tubs *solidVSDF = new G4Tubs("solidVSDF",r_min_vsdf,r_max_vsdf,len_vsdf/2.,startPhi,dPhi);

   // vacuum step-down flange [item 5, vacuum] 
   G4Tubs *solidVSDF_vac = new G4Tubs("solidVSDF_vac",0,r_min_vsdf,len_vsdf/2.,startPhi,dPhi);

   // union solid
   // - start with large flange + cone [origin is center of LGF]  
   G4double zz = 0.5*len_lgf + 0.5*len_lgc; 
   G4ThreeVector P_0 = G4ThreeVector(0,0,-zz);
   G4UnionSolid *upstrPipe = new G4UnionSolid("lgf_lgc",solidLGF,solidConeLG,0,P_0); 
   // - attach m2b 
   zz  = 0.5*len_lgf + len_lgc + 0.5*len_m2b; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b",upstrPipe,solidTubeM2b,0,P_0); 
   // - vacuum tube support ring [item 12]  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + 0.5*len_vtsr; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr",upstrPipe,solidVTSRing,0,P_0); 
   // - attach m2 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + 0.5*len_m2; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2",upstrPipe,solidTubeM2,0,P_0);
   // - attach aperture plate flange [item 14] 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + 0.5*len_apf; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf",upstrPipe,solidAPF,0,P_0);
   // - attach flange with o-ring [item 9] 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + 0.5*len_for; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for",upstrPipe,solidFORing,0,P_0);
   // - attach m6b 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + 0.5*len_m6b; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b",upstrPipe,solidTubeM6b,0,P_0);
   // - vacuum tube support ring [item 12]  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + 0.5*len_vtsr; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr",upstrPipe,solidVTSRing,0,P_0);
   // - attach m6 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + 0.5*len_m6; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6",upstrPipe,solidTubeM6,0,P_0);
   // - attach small cone 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6 
       + 0.5*len_smc; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6_smc",upstrPipe,solidConeSM,0,P_0);
   // - attach vacuum window tube  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6 
       + len_smc + 0.5*len_vwt; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6_smc_vwt",upstrPipe,solidVWTube,0,P_0);
   // - attach vacuum window stepdown flange  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6 
       + len_smc + len_vwt + 0.5*len_vsdf; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("upstreamPipe",upstrPipe,solidVSDF,0,P_0);
 
   // this is a check  
   G4double TOT_LEN_SUM = len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b 
                        + len_vtsr + len_m6 + len_smc + len_vwt + len_vsdf;

   if(TOT_LEN_SUM!=TOTAL_LENGTH){ 
      std::cout << "TOTAL_LENGTH = " << TOTAL_LENGTH/inch 
                << " inches, sum of parts = " << TOT_LEN_SUM/inch << " inches" << std::endl;
   }
 
   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Green() ); 
   // vis->SetForceWireframe(true); 
   
   // logical volume
   G4LogicalVolume *tubeLV = new G4LogicalVolume(upstrPipe,Aluminum_5052,"beamDump_usPipe_LV");
   tubeLV->SetVisAttributes(vis);

   // placement
   G4double z = z0 + 0.5*len_lgf;  // upstream face is at z0 
   G4ThreeVector P = G4ThreeVector(0.,0.,z);
   G4RotationMatrix *rm = new G4RotationMatrix();
   rm->rotateY(180*deg); 
   new G4PVPlacement(rm,                       // no rotation
	             P,                        // location in mother volume 
	             tubeLV,                   // its logical volume                         
	             "beamDump_usPipe_PHY",    // its name
	             logicMother,              // its mother  volume
	             false,                    // boolean operation? 
	             0,                        // copy number
	             fCheckOverlaps);          // checking overlaps    

   // union solid [vacuum] 
   // - start with large flange + cone [origin is center of LGF]  
   zz = 0.5*len_lgf + 0.5*len_lgc; 
   P_0 = G4ThreeVector(0,0,-zz);
   G4UnionSolid *upstrPipe_vac = new G4UnionSolid("lgf_lgc",solidLGF_vac,solidConeLG_vac,0,P_0); 
   // - attach m2b 
   zz  = 0.5*len_lgf + len_lgc + 0.5*len_m2b; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b",upstrPipe_vac,solidTubeM2b_vac,0,P_0); 
   // - vacuum tube support ring [item 12]  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + 0.5*len_vtsr; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr",upstrPipe_vac,solidVTSRing_vac,0,P_0); 
   // - attach m2 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + 0.5*len_m2; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2",upstrPipe_vac,solidTubeM2_vac,0,P_0);
   // - attach aperture plate flange [item 14] 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + 0.5*len_apf; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf",upstrPipe_vac,solidAPF_vac,0,P_0);
   // - attach flange with o-ring [item 9] 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + 0.5*len_for; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for",upstrPipe_vac,solidFORing_vac,0,P_0);
   // - attach m6b 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + 0.5*len_m6b; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b",upstrPipe_vac,solidTubeM6b_vac,0,P_0);
   // - vacuum tube support ring [item 12]  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + 0.5*len_vtsr; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr",upstrPipe_vac,solidVTSRing_vac,0,P_0);
   // - attach m6 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + 0.5*len_m6; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6",upstrPipe_vac,solidTubeM6_vac,0,P_0);
   // - attach small cone 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6 
       + 0.5*len_smc; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6_smc",upstrPipe_vac,solidConeSM_vac,0,P_0);
   // - attach vacuum window tube  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6 
       + len_smc + 0.5*len_vwt; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6_smc_vwt",upstrPipe_vac,solidVWTube_vac,0,P_0);
   // - attach vacuum window stepdown flange  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6 
       + len_smc + len_vwt + 0.5*len_vsdf; 
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("upstreamPipe_vac",upstrPipe_vac,solidVSDF_vac,0,P_0);

   // visualization
   G4VisAttributes *vis_vac = new G4VisAttributes(); 
   vis_vac->SetColour( G4Colour::White() ); 
   vis_vac->SetForceWireframe(true); 
   
   // logical volume
   std::string name; 
   G4double a=0,density=0;
   G4Material *Vacuum          = new G4Material(name="Vacuum", z=1., a=1.0*g/mole, density=1e-9*g/cm3);
   G4LogicalVolume *tubeLV_vac = new G4LogicalVolume(upstrPipe_vac,Vacuum,"beamDump_usPipe_LV");
   tubeLV_vac->SetVisAttributes(vis_vac);

   new G4PVPlacement(rm,                       // no rotation
	             P,                        // location in mother volume 
	             tubeLV_vac,               // its logical volume                         
	             "beamDump_usPipeVac_PHY", // its name
	             logicMother,              // its mother  volume
	             false,                    // boolean operation? 
	             0,                        // copy number
	             fCheckOverlaps);          // checking overlaps    

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildBeamDump_DownstreamPipe(G4LogicalVolume *logicMother,G4double z0){
   // Hall A Beam Dump: Pipe downstream of ISO Weldment
   // Drawing: JL0011756_27020E0145 [modified]

   G4double inch     = 25.4*mm; 
   G4double startPhi = 0.*deg; 
   G4double dPhi     = 360.*deg; 
   G4double TOTAL_LENGTH = 328.04*inch; 

   // main tube
   G4double r_min_m  = 0.5*12.11*inch;   // FIXME: This is arbitrary! 
   G4double r_max_m  = 0.5*12.12*inch;   
   G4double len_m    = 327.66*inch; 
   G4Tubs *solidTubeM = new G4Tubs("solidTube",r_min_m,r_max_m,len_m/2.,startPhi,dPhi);

   // vacuum insert 
   G4Tubs *solidVacuumInsert = new G4Tubs("solidVacuumInsert",0.,r_min_m,TOTAL_LENGTH/2.,startPhi,dPhi); 

   // bookend flange [upstream]  
   G4double r_min_us = 0.5*12.11*inch; 
   G4double r_max_us = 0.5*14.00*inch;   // FIXME: This is arbitrary!
   G4double len_us   = 0.188*inch;  
   G4Tubs *solidTubeUS = new G4Tubs("solidTubeUS",r_min_us,r_max_us,len_us/2.,startPhi,dPhi);

   // bookend flange [downstream]  
   G4double r_min_ds = 0.5*12.11*inch; 
   G4double r_max_ds = 0.5*14.00*inch;   // FIXME: This is arbitrary!
   G4double len_ds   = 0.188*inch;  
   G4Tubs *solidTubeDS = new G4Tubs("solidTubeDS",r_min_ds,r_max_ds,len_ds/2.,startPhi,dPhi);

   // union solid 
   // - start with upstream flange and main tube 
   G4double zz = 0.5*len_us + 0.5*len_m;
   G4ThreeVector P_0 = G4ThreeVector(0,0,-zz);  
   G4UnionSolid *dwnstrPipe = new G4UnionSolid("us_m",solidTubeUS,solidTubeM,0,P_0);
   // - downstream flange 
   zz = 0.5*len_us + len_m + 0.5*len_ds;
   P_0 = G4ThreeVector(0,0,-zz);  
   dwnstrPipe = new G4UnionSolid("dwnstrPipe",dwnstrPipe,solidTubeDS,0,P_0);
 
   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Yellow() ); 
   // vis->SetForceWireframe(); 
   
   // logical volume
   auto material  = G4Material::GetMaterial("G4_Al"); // might not be aluminum... 
   G4LogicalVolume *tubeLV = new G4LogicalVolume(dwnstrPipe,material,"beamDump_dsPipe_LV");
   tubeLV->SetVisAttributes(vis); 

   // placement
   G4double delta = 5.0*cm; // FIXME: arbitrary!  
   G4double Z = z0 + TOTAL_LENGTH + 0.5*len_us; // upstream face is at z0
   G4ThreeVector P = G4ThreeVector(0,0,Z);
   new G4PVPlacement(0,                        // no rotation
	             P,                        // location in mother volume 
	             tubeLV,                   // its logical volume                         
	             "beamDump_dsPipe_PHY",    // its name
	             logicMother,              // its mother  volume
	             false,                    // boolean operation? 
	             0,                        // copy number
	             fCheckOverlaps);          // checking overlaps   

   // vacuum insert 
   G4VisAttributes *visV = new G4VisAttributes(); 
   visV->SetForceWireframe();

   std::string name;
   G4double z=0,a=0,density=0;
 
   // logical volume
   G4Material *Vacuum     = new G4Material(name="Vacuum", z=1., a=1.0*g/mole, density=1e-9*g/cm3);
   G4LogicalVolume *vacLV = new G4LogicalVolume(solidVacuumInsert,Vacuum,"vacuum_dsPipe_LV");
   vacLV->SetVisAttributes(visV);  

   Z = z0 + TOTAL_LENGTH/2.;
   P = G4ThreeVector(0,0,Z);
   new G4PVPlacement(0,                        // no rotation
	             P,                        // location in mother volume 
	             vacLV,                    // its logical volume                         
	             "vacuum_dsPipe_PHY",      // its name
	             logicMother,              // its mother  volume
	             false,                    // boolean operation? 
	             0,                        // copy number
	             fCheckOverlaps);          // checking overlaps   

}
//______________________________________________________________________________
void BDDetectorConstruction::MakeBeamExit_TargetToMidPipe(G4LogicalVolume *logicMother,G4double z0){
   // SBS exit beam pipe.  This is immediately upstream of the mid pipe
   // z0 = position of upstream face of this part
   // Added by D. Flay (JLab) in Sept 2020
   // Drawings: 
   // - Target to mid pipe: ARC10540  

   G4double inch     = 2.54*cm;
   G4double startPhi = 0.*deg;
   G4double dPhi     = 360.*deg;
   G4double len_fl   = 0.5*inch; // FIXME: arbitrary
   G4double TOTAL_LENGTH = 0;

   // visualization 
   G4VisAttributes *AlColor = new G4VisAttributes( G4Colour(0.3,0.3,1.0) );
   // G4VisAttributes *vis_vac = new G4VisAttributes( G4Colour(0.1,0.5,0.9) );
   G4VisAttributes *vis_vac = new G4VisAttributes();

   // pipe [part 11] 
   G4double r_min_11          = 0.5*36.00*inch;
   G4double r_max_11          = 0.5*42.00*inch;  // based on 0.060 x 3 corrugation 
   G4double len_11            = 216.375*inch;
   G4Tubs *solidTube11        = new G4Tubs("solidTube11",r_min_11,r_max_11,len_11/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_11;
   // vacuum insert  
   G4Tubs *solidTube11_vac    = new G4Tubs("solidTube11_vac",0,r_min_11,len_11/2.,startPhi,dPhi);

   // flange [part 9] 
   G4double r_min_09          = 0.5*24.00*inch;
   G4double r_max_09          = 0.5*40.00*inch;
   G4double len_09            = len_fl;
   G4Tubs *solidTube09        = new G4Tubs("solidTube09",r_min_09,r_max_09,len_09/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_09;
   // vacuum insert  
   G4Tubs *solidTube09_vac    = new G4Tubs("solidTube09_vac",0,r_min_09,len_09/2.,startPhi,dPhi);

   // flange [part 10] 
   G4double r_min_10          = 0.5*36.00*inch;
   G4double r_max_10          = 0.5*43.00*inch;
   G4double len_10            = len_fl;
   G4Tubs *solidTube10        = new G4Tubs("solidTube10",r_min_10,r_max_10,len_10/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_10;
   // vacuum insert 
   G4Tubs *solidTube10_vac    = new G4Tubs("solidTube10_vac",0,r_min_10,len_10/2.,startPhi,dPhi);

   // pipe [part 8] 
   G4double r_min_08          = 0.5*24.00*inch;
   G4double r_max_08          = 0.5*29.334*inch;  // based on 0.5 x 2-2/3 corrugation  
   G4double len_08            = 217.00*inch;
   G4Tubs *solidTube08        = new G4Tubs("solidTube08",r_min_08,r_max_08,len_08/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_08;
   // vacuum insert 
   G4Tubs *solidTube08_vac    = new G4Tubs("solidTube08_vac",0,r_min_08,len_08/2.,startPhi,dPhi);

   // flange [part 6] 
   G4double r_min_06          = 0.5*12.00*inch;
   G4double r_max_06          = 0.5*26.00*inch;
   G4double len_06            = len_fl;
   G4Tubs *solidTube06        = new G4Tubs("solidTube06",r_min_06,r_max_06,len_06/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_06;
   // vacuum insert  
   G4Tubs *solidTube06_vac    = new G4Tubs("solidTube06_vac",0,r_min_06,len_06/2.,startPhi,dPhi);

   // pipe [part 5] 
   G4double r_min_05          = 0.5*12.00*inch;
   G4double r_max_05          = 0.5*17.334*inch; // based on 0.5 x 2-2/3 corrugation  
   G4double len_05            = 41.8125*inch;
   G4Tubs *solidTube05        = new G4Tubs("solidTube05",r_min_05,r_max_05,len_05/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_05;
   // vacuum insert
   G4Tubs *solidTube05_vac    = new G4Tubs("solidTube05_vac",0,r_min_05,len_05/2.,startPhi,dPhi);

   // flange [part 7] 
   G4double r_min_07          = 0.5*10.00*inch;
   G4double r_max_07          = 0.5*13.00*inch;
   G4double len_07            = 0.125*inch;
   G4Tubs *solidTube07        = new G4Tubs("solidTube07",r_min_07,r_max_07,len_07/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_07;
   // vacuum insert  
   G4Tubs *solidTube07_vac    = new G4Tubs("solidTube07_vac",0,r_min_07,len_07/2.,startPhi,dPhi);

   // pipe [part 3]  
   G4double r_min_03          = 0.5*10.00*inch;
   G4double r_max_03          = 0.5*10.25*inch;
   G4double len_03            = 12.00*inch;
   G4Tubs *solidTube03        = new G4Tubs("solidTube03",r_min_03,r_max_03,len_03/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_03;
   // vacuum insert
   G4Tubs *solidTube03_vac    = new G4Tubs("solidTube03_vac",0,r_min_03,len_03/2.,startPhi,dPhi);

   // flange [part 4] 
   G4double r_min_04          = 0.5*8.25*inch;
   G4double r_max_04          = 0.5*10.00*inch;
   G4double len_04            = 0.125*inch;
   G4Tubs *solidTube04        = new G4Tubs("solidTube04",r_min_04,r_max_04,len_04/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_04;
   // vacuum insert  
   G4Tubs *solidTube04_vac    = new G4Tubs("solidTube04_vac",0,r_min_04,len_04/2.,startPhi,dPhi);

   // pipe [part 2]  
   G4double r_min_02          = 0.5*8.00*inch;
   G4double r_max_02          = 0.5*8.25*inch;
   G4double len_02            = 36.875*inch;
   G4Tubs *solidTube02        = new G4Tubs("solidTube02",r_min_02,r_max_02,len_02/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_02;
   // vacuum insert
   G4Tubs *solidTube02_vac    = new G4Tubs("solidTube02_vac",0,r_min_02,len_02/2.,startPhi,dPhi);

   // flange [part 1] 
   G4double r_min_01          = 0.5*8.25*inch;
   G4double r_max_01          = 0.5*10.00*inch;
   G4double len_01            = 0.125*inch;
   G4Tubs *solidTube01        = new G4Tubs("solidTube01",r_min_01,r_max_01,len_01/2.,startPhi,dPhi);
   TOTAL_LENGTH += len_01;
   // vacuum insert  
   G4Tubs *solidTube01_vac    = new G4Tubs("solidTube01_vac",0,r_min_01,len_01/2.,startPhi,dPhi);

   // union: put it all together
   // - start with 01 and 02  
   G4double zp = 0.5*len_01 + 0.5*len_02;
   G4ThreeVector P = G4ThreeVector(0,0,zp);
   G4UnionSolid *tgtToMidPipe = new G4UnionSolid("t01",solidTube01,solidTube02,0,P);
   // - add 04 
   zp = 0.5*len_01 + len_02 + 0.5*len_04;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe = new G4UnionSolid("t04",tgtToMidPipe,solidTube04,0,P);
   // - add 03 
   zp = 0.5*len_01 + len_02 + len_04 + 0.5*len_03;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe = new G4UnionSolid("t03",tgtToMidPipe,solidTube03,0,P);
   // - add 07 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + 0.5*len_07;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe = new G4UnionSolid("t07",tgtToMidPipe,solidTube07,0,P);
   // - add 05 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + 0.5*len_05;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe = new G4UnionSolid("t05",tgtToMidPipe,solidTube05,0,P);
   // - add 06 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + 0.5*len_06;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe = new G4UnionSolid("t06",tgtToMidPipe,solidTube06,0,P);
   // - add 08 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + 0.5*len_08;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe = new G4UnionSolid("t08",tgtToMidPipe,solidTube08,0,P);
   // - add 09 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + 0.5*len_09;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe = new G4UnionSolid("t09",tgtToMidPipe,solidTube09,0,P);
   // - add 11 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + len_09 + 0.5*len_11;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe = new G4UnionSolid("t11",tgtToMidPipe,solidTube11,0,P);
   // - add 10 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + len_09 + len_11 + 0.5*len_10;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe = new G4UnionSolid("tgtToMidPipe",tgtToMidPipe,solidTube10,0,P);

   auto material    = G4Material::GetMaterial("G4_Al");
   auto materialVac = G4Material::GetMaterial("Galactic");

   G4LogicalVolume *tgtMP_LV = new G4LogicalVolume(tgtToMidPipe,material,"tgtMP_LV");
   tgtMP_LV->SetVisAttributes(AlColor);

   G4double z = z0 + 0.5*len_01;  // upstream face at z0
   G4ThreeVector Pz = G4ThreeVector(0,0,z);
   new G4PVPlacement(0,                // no rotation
                     Pz,               // location in mother volume 
                     tgtMP_LV,         // its logical volume                         
                     "tgtMidPipe_PHY", // its name
                     logicMother,      // its mother  volume
                     true,             // boolean operation? 
                     0,                // copy number
                     true);            // checking overlaps 

   // union: put it all together [vacuum] 
   // - start with 01 and 02  
   zp = 0.5*len_01 + 0.5*len_02;
   P = G4ThreeVector(0,0,zp);
   G4UnionSolid *tgtToMidPipe_vac = new G4UnionSolid("t01v",solidTube01_vac,solidTube02_vac,0,P);
   // - add 04 
   zp = 0.5*len_01 + len_02 + 0.5*len_04;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe_vac = new G4UnionSolid("t04v",tgtToMidPipe_vac,solidTube04_vac,0,P);
   // - add 03 
   zp = 0.5*len_01 + len_02 + len_04 + 0.5*len_03;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe_vac = new G4UnionSolid("t03v",tgtToMidPipe_vac,solidTube03_vac,0,P);
   // - add 07 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + 0.5*len_07;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe_vac = new G4UnionSolid("t07v",tgtToMidPipe_vac,solidTube07_vac,0,P);
   // - add 05 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + 0.5*len_05;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe_vac = new G4UnionSolid("t05v",tgtToMidPipe_vac,solidTube05_vac,0,P);
   // - add 06 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + 0.5*len_06;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe_vac = new G4UnionSolid("t06v",tgtToMidPipe_vac,solidTube06_vac,0,P);
   // - add 08 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + 0.5*len_08;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe_vac = new G4UnionSolid("t08v",tgtToMidPipe_vac,solidTube08_vac,0,P);
   // - add 09 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + 0.5*len_09;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe_vac = new G4UnionSolid("t09v",tgtToMidPipe_vac,solidTube09_vac,0,P);
   // - add 11 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + len_09 + 0.5*len_11;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe_vac = new G4UnionSolid("t11v",tgtToMidPipe_vac,solidTube11_vac,0,P);
   // - add 10 
   zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + len_09 + len_11 + 0.5*len_10;
   P  = G4ThreeVector(0,0,zp);
   tgtToMidPipe_vac = new G4UnionSolid("tgtToMidPipe",tgtToMidPipe_vac,solidTube10_vac,0,P);

   G4LogicalVolume *tgtMP_vac_LV = new G4LogicalVolume(tgtToMidPipe_vac,materialVac,"tgtMP_vac_LV");
   tgtMP_vac_LV->SetVisAttributes(vis_vac);

   new G4PVPlacement(0,                    // no rotation
                     Pz,                   // location in mother volume 
                     tgtMP_vac_LV,         // its logical volume                         
                     "tgtMidPipe_vac_PHY", // its name
                     logicMother,          // its mother  volume
                     true,                 // boolean operation? 
                     0,                    // copy number
                     true);                // checking overlaps  

}
//______________________________________________________________________________
void BDDetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // Sensitive detectors
  // auto absoSD     = new BDCalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  // G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  // SetSensitiveDetector("AbsoLV",absoSD);

  // auto diffuserSD = new BeamDiffuserSD("DiffuserSD","DiffuserHitsCollection"); 
  // G4SDManager::GetSDMpointer()->AddNewDetector(diffuserSD);
  // SetSensitiveDetector("plateLV",diffuserSD); 

  // auto gapSD      = new BDCalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  // G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  // SetSensitiveDetector("GapLV",gapSD);

  // ion chamber 
  auto icSD = new IonChamberSD("IonChamberSD","IonChamberHitsCollection"); 
  G4SDManager::GetSDMpointer()->AddNewDetector(icSD); 
  SetSensitiveDetector("ic_LV",icSD);  

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
void BDDetectorConstruction::GetRotatedCoordinates(std::vector<G4double> RA,
                                                   std::vector<G4double> P,
                                                   std::vector<G4double> &PP){
   // compute rotated coordinates based on arbitrary rotations in all three dimensions
   // inputs 
   // - RA = rotation angles [x,y,z]
   // - P  = initial coordinates [x,y,z] 
   // - PP = rotated coordinates [x,y,z] 

   // rotation angles
   G4double COS_G = cos(RA[0]); G4double COS_B = cos(RA[1]); G4double COS_A = cos(RA[2]);
   G4double SIN_G = sin(RA[0]); G4double SIN_B = sin(RA[1]); G4double SIN_A = sin(RA[2]);
   // compute rotated coordinates 
   G4double xp = COS_A*COS_B*P[0] + (COS_A*COS_B*SIN_G - SIN_A*COS_G)*P[1] + (COS_A*SIN_B*COS_G + SIN_A*SIN_G)*P[2];
   G4double yp = SIN_A*COS_B*P[0] + (SIN_A*SIN_B*SIN_G + COS_A*COS_G)*P[1] + (SIN_A*SIN_B*COS_G - COS_A*SIN_G)*P[2];
   G4double zp =      -SIN_B*P[0] +                       COS_B*SIN_G*P[1] +                       COS_B*COS_G*P[2];
   // fill output vector 
   PP.push_back(xp); 
   PP.push_back(yp); 
   PP.push_back(zp);
 
   // print to screen
   char msg[200]; 
   const int N = P.size(); 
   for(int i=0;i<N;i++){
      sprintf(msg,"[BDDetectorConstruction::GetRotatedCoordinates]: P[%d] = %.3lf mm, P'[%d] = %.3lf mm",i,P[i]/mm,i,PP[i]/mm); 
      std::cout << msg << std::endl;
   }
 
}
//______________________________________________________________________________
void BDDetectorConstruction::BuildCollimators(G4LogicalVolume *logicMother,G4double z0){
   // Collimators near target (beam left) 
   // Based on drawings from Sebastian Seeds (UConn), derived from Bert Metzger's JT file 
   BuildCollimator_A(logicMother,z0); 
   BuildCollimator_B(logicMother,z0); 
   BuildCollimator_C(logicMother,z0);

   BuildGEnTarget_CollimatorTable(logicMother);  
}
//______________________________________________________________________________
void BDDetectorConstruction::BuildCollimator_A(G4LogicalVolume *logicMother,G4double z0){
   // From drawings made by Sebastian Seeds (UConn) based on JT file
   // - Collimator A1: CollimatorA_1_drawing.JPG
   // - Collimator A2: CollimatorA_2_drawing.JPG
   // - Offsets and rotations: CollimatorA_xzoffset.JPG
   // Note: Collimators are on the LEFT side of the beam, next to the target  

   double inch   = 2.54*cm;
 
   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Red() ); 

   // define materials and logical volume 
   G4Material *Aluminum = G4Material::GetMaterial("G4_Al");

   // collimator A1: right-angle trapezoid  
   // double zl_a1     = 2.247*inch;  // length along z
   // double yl_a1     = 6.000*inch;  // length along y 
   // double xs_a1     = 1.614*inch;  // length along x (short side)
   // double xl_a1     = 3.391*inch;  // length along x (long side)  
   // G4Trap *raSolid_A1 = new G4Trap("raSolid_A1",zl_a1,yl_a1,xl_a1,xs_a1);

   G4double z       = 6.000*inch;  // z length
   G4double y_mz    = 1.614*inch;  // y length at -z
   G4double y_pz    = 3.391*inch + (3.391*inch-y_mz);  // y length at +z should be 5.456 inch 
   G4double x_my_mz = 2.247*inch;  // x length at -y, -z  
   G4double x_py_mz = 2.247*inch;  // x length at +y, -z
   G4double x_my_pz = 3.189*inch;  // x length at -y, +z  
   G4double x_py_pz = 3.189*inch;  // x length at +y, +z
   G4double theta   = 0.*deg;  
   G4double phi     = 0.*deg; 
   G4double alpha1  = 0.5*(x_py_mz-x_my_mz)/y_pz; 
   G4double alpha2  = alpha1;
   G4Trap *trapA1   = new G4Trap("trapA1",z/2.,theta,phi,y_mz/2.,x_my_mz/2.,x_py_mz/2.,alpha1,y_pz/2.,x_my_pz/2.,x_py_pz/2.,alpha2); 

   // G4LogicalVolume *trapA1_LV = new G4LogicalVolume(trapA1,Aluminum,"trapA1_LV");
   // trapA1_LV->SetVisAttributes(vis); 

   // new G4PVPlacement(0,                         // rotation
   //                   G4ThreeVector(0,0,0),      // position 
   //                   trapA1_LV,                   // logical volume   
   //                   "trapA1_PHY",                // physical name 
   //                   logicMother,               // logical mother
   //                   false,                     // boolean? 
   //                   0,                         // copy no 
   //                   fCheckOverlaps);           // check overlaps

   // cut the bottom off 
   G4Box *a1Cut = new G4Box("a1Cut",2.*inch,2.*inch,7.*inch);
   G4ThreeVector P_a1Cut = G4ThreeVector(0,-y_mz-1.2*inch,0); 
   // G4LogicalVolume *a1Cut_LV = new G4LogicalVolume(a1Cut,Aluminum,"a1Cut_LV");

   // new G4PVPlacement(0,                         // rotation
   //                   P_a1Cut,                   // position 
   //                   a1Cut_LV,                  // logical volume   
   //                   "a1Cut_PHY",               // physical name 
   //                   logicMother,               // logical mother
   //                   false,                     // boolean? 
   //                   0,                         // copy no 
   //                   fCheckOverlaps);           // check overlaps

   G4SubtractionSolid *raSolid_A1 = new G4SubtractionSolid("raSolid_A1",trapA1,a1Cut,0,P_a1Cut);

   // cut in A1 
   G4Box *cutA1 = new G4Box("cutA1",0.250*inch/2.,0.250*inch/2.,4.000*inch/2.);
   G4ThreeVector Pca1 = G4ThreeVector(0,-0.682*inch,-2.*inch);
   G4SubtractionSolid *colSolid_A1 = new G4SubtractionSolid("colSolid_A1",raSolid_A1,cutA1,0,Pca1);

   // collimator A2: right-angle trapezoid
   // double z_a2      = 6.000*inch;
   // double x1_a2     = 2.247*inch;  // upstream x size 
   // double x2_a2     = 2.247*inch; // 3.189*inch;  // downstream x size 
   // double y_a2      = 2.500*inch;  // thickness  
   // G4Trd *colSolid_A2 = new G4Trd("colSolid_A2",x1_a2/2.,x2_a2/2.,y_a2/2.,y_a2/2.,z_a2/2.); 

   z       = 6.000*inch;  // z length
   y_mz    = 2.500*inch;  // y length at -z
   y_pz    = 2.500*inch;  // y length at +z should be 5.456 inch 
   x_my_mz = 2.247*inch;  // x length at -y, -z  
   x_py_mz = 2.247*inch;  // x length at +y, -z
   x_my_pz = 3.189*inch;  // x length at -y, +z  
   x_py_pz = 3.189*inch;  // x length at +y, +z
   theta   = 0.*deg;  
   phi     = 0.*deg; 
   alpha1  = 0.5*(x_py_mz-x_my_mz)/y_pz; 
   alpha2  = alpha1;
   G4Trap *colSolid_A2 = new G4Trap("colSolid_A2",z/2.,theta,phi,y_mz/2.,x_my_mz/2.,x_py_mz/2.,alpha1,y_pz/2.,x_my_pz/2.,x_py_pz/2.,alpha2); 

   // G4LogicalVolume *trapA1b_LV = new G4LogicalVolume(colSolid_A2,Aluminum,"trapA1b_LV");
   // trapA1b_LV->SetVisAttributes(vis); 

   // new G4PVPlacement(0,                         // rotation
   //                   G4ThreeVector(0,0,0),      // position 
   //                   trapA1b_LV,                // logical volume   
   //                   "trapA1b_PHY",             // physical name 
   //                   logicMother,               // logical mother
   //                   false,                     // boolean? 
   //                   0,                         // copy no 
   //                   fCheckOverlaps);           // check overlaps

   // // cut the bottom off 
   // G4Box *a1Cut = new G4Box("a1Cut",2.*inch,2.*inch,7.*inch);
   // G4ThreeVector P_a1Cut = G4ThreeVector(0,-y_mz-1.2*inch,0); 
   // G4LogicalVolume *a1Cut_LV = new G4LogicalVolume(a1Cut,Aluminum,"a1Cut_LV");

   // new G4PVPlacement(0,                         // rotation
   //                   P_a1Cut,                   // position 
   //                   a1Cut_LV,                  // logical volume   
   //                   "a1Cut_PHY",               // physical name 
   //                   logicMother,               // logical mother
   //                   false,                     // boolean? 
   //                   0,                         // copy no 
   //                   fCheckOverlaps);           // check overlaps
 

   // union of these objects.  use A2 as the reference point since it's easier
   // - all positions and rotations are relative to A2 center  
   G4ThreeVector P21     = G4ThreeVector(0,y_mz-0.45*inch,0); 
   G4RotationMatrix *r21 = new G4RotationMatrix();
   // r21->rotateX(90*deg); r21->rotateY(90*deg); r21->rotateZ(0);  
   G4UnionSolid *col_A   = new G4UnionSolid("col_A",colSolid_A2,colSolid_A1,r21,P21); 

   // define materials and logical volume 
   G4LogicalVolume *col_A_LV = new G4LogicalVolume(col_A,Aluminum,"col_A_LV"); 
   col_A_LV->SetVisAttributes(vis); 

   // FIXME: Get placement right!
   // placement of the union object in the Hall coordinate system 
   // position 
   double X = 3.75*inch; double Y = -1.243*inch; double Z = z0 - 7.495*inch; 
   G4ThreeVector P = G4ThreeVector(X,Y,Z); 
   // rotation 
   double RX = 0.; double RY = -34.62*deg; double RZ = 0.;
   G4RotationMatrix *rm = new G4RotationMatrix();
   rm->rotateX(RX); rm->rotateY(RY); rm->rotateZ(RZ); 

   new G4PVPlacement(rm,                         // rotation
	             P,                          // position 
                     col_A_LV,                   // logical volume   
                     "col_A_PHY",                // physical name 
                     logicMother,                // logical mother
                     true,                       // boolean? 
                     0,                          // copy no 
                     fCheckOverlaps);            // check overlaps

   // now build the collimator base 
   G4double xb_len = 5.114*inch; 
   G4double yb_len = 0.750*inch; 
   G4double zb_len = 6.250*inch;
   G4Box *solidBase = new G4Box("solidBase",xb_len/2.,yb_len/2.,zb_len/2.);
   G4LogicalVolume *solidBase_LV = new G4LogicalVolume(solidBase,Aluminum,"solidBase_LV");
 
   // new G4PVPlacement(0,                          // rotation
   //                   G4ThreeVector(0,0,0),       // position 
   //                   solidBase_LV,               // logical volume   
   //                   "solidBase_PHY",            // physical name 
   //                   logicMother,                // logical mother
   //                   true,                       // boolean? 
   //                   0,                          // copy no 
   //                   fCheckOverlaps);            // check overlaps

   // base cut 1
   G4double xcc_len = 3.*inch;  
   G4double ycc_len = 2.*inch;  
   G4double zcc_len = 4.749*inch; 
   G4Box *solidBase_cut1 = new G4Box("solidBase_cut1",xcc_len/2.,ycc_len/2.,zcc_len/2.); 
   G4LogicalVolume *solidBase_cut1_LV = new G4LogicalVolume(solidBase_cut1,Aluminum,"solidBase_cut1_LV");
   solidBase_cut1_LV->SetVisAttributes(vis);  

   G4RotationMatrix *rmc1 = new G4RotationMatrix();
   rmc1->rotateY(39.11*deg); 

   G4double XX1 = 2.207*inch; // xb_len/2. - 0.35*inch; 
   G4double YY1 = 0; 
   G4double ZZ1 = 2.25*inch; 
   G4ThreeVector P_bc1 = G4ThreeVector(XX1,YY1,ZZ1);

   // new G4PVPlacement(rmc1,                       // rotation
   //                   P_bc1,                      // position 
   //                   solidBase_cut1_LV,          // logical volume   
   //                   "solidBase_cut1_PHY",       // physical name 
   //                   logicMother,                // logical mother
   //                   true,                       // boolean? 
   //                   0,                          // copy no 
   //                   fCheckOverlaps);            // check overlaps

   // base cut2 
   G4double xcc2_len = 3.*inch;  
   G4double ycc2_len = 2.*inch;  
   G4double zcc2_len = 4.968*inch; 
   G4Box *solidBase_cut2 = new G4Box("solidBase_cut2",xcc2_len/2.,ycc2_len/2.,zcc2_len/2.); 
   G4LogicalVolume *solidBase_cut2_LV = new G4LogicalVolume(solidBase_cut2,Aluminum,"solidBase_cut2_LV");
   solidBase_cut2_LV->SetVisAttributes(vis);  

   G4RotationMatrix *rmc2 = new G4RotationMatrix();
   rmc2->rotateY(25.53*deg); 

   G4double XX2 = -2.85*inch; // xb_len/2. - 0.35*inch; 
   G4double YY2 = 0; 
   G4double ZZ2 = -1.5*inch; 
   // std::cout << XX/inch << " " << YY/inch << " " << ZZ/inch << std::endl;
   G4ThreeVector P_bc2 = G4ThreeVector(XX2,YY2,ZZ2);

   // new G4PVPlacement(rmc2,                       // rotation
   //                   P_bc2,                      // position 
   //                   solidBase_cut2_LV,          // logical volume   
   //                   "solidBase_cut2_PHY",       // physical name 
   //                   logicMother,                // logical mother
   //                   false,                      // boolean? 
   //                   0,                          // copy no 
   //                   fCheckOverlaps);            // check overlaps

   // base cut3 
   G4double xcc3_len = 1.600*inch;  
   G4double ycc3_len = 0.188*inch;  
   G4double zcc3_len = 7.000*inch; 
   G4Box *solidBase_cut3 = new G4Box("solidBase_cut3",xcc3_len/2.,ycc3_len/2.,zcc3_len/2.); 
   G4LogicalVolume *solidBase_cut3_LV = new G4LogicalVolume(solidBase_cut3,Aluminum,"solidBase_cut3_LV");
   solidBase_cut3_LV->SetVisAttributes(vis);  

   G4double XX3 = 1.817*inch;  
   G4double YY3 = yb_len/2. - ycc3_len/2.; 
   G4double ZZ3 = 0.*inch; 
   G4ThreeVector P_bc3 = G4ThreeVector(XX3,YY3,ZZ3);

   // subtraction solid 
   // base - cut 1
   G4SubtractionSolid *colBase   = new G4SubtractionSolid("colBase_1",solidBase,solidBase_cut1,rmc1,P_bc1);  
   // cut 2 
   colBase   = new G4SubtractionSolid("colBase_12",colBase,solidBase_cut2,rmc2,P_bc2); 
   // cut 3 
   colBase   = new G4SubtractionSolid("colBase"   ,colBase,solidBase_cut3,0,P_bc3); 

   G4VisAttributes *visBase = new G4VisAttributes(); 
   visBase->SetColour( G4Colour::Blue() ); 
 
   G4LogicalVolume *colBase_A_LV = new G4LogicalVolume(colBase,Aluminum,"colBase_A_LV");
   colBase_A_LV->SetVisAttributes(visBase); 

   // placement 
   G4RotationMatrix *rmb = new G4RotationMatrix(); 
   rmb->rotateX(180.*deg);

   G4ThreeVector Pb = G4ThreeVector(3.693*inch,-2.889*inch,-7.125*inch); 

   new G4PVPlacement(rmb,                        // rotation
                     Pb,                         // position 
                     colBase_A_LV,               // logical volume   
                     "colBase_A_PHY",            // physical name 
                     logicMother,                // logical mother
                     true,                       // boolean? 
                     0,                          // copy no 
                     fCheckOverlaps);            // check overlaps

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildCollimator_B(G4LogicalVolume *logicMother,G4double z0){
   // From drawings made by Sebastian Seeds (UConn) based on JT file
   // - Collimator B: CollimatorB_drawing.JPG
   // - Offsets and rotations: CollimatorB_xzoffset.JPG
   // Note: Collimators are on the LEFT side of the beam, next to the target  

   double inch   = 2.54*cm;
 
   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Red() ); 

   // define materials and logical volume 
   G4Material *Aluminum = G4Material::GetMaterial("G4_Al");

   // collimator B: right-angle trapezoid  
   // double zl_a1     = 0.815*inch;  // length along z
   // double yl_a1     = 3.397*inch;  // length along y 
   // double xs_a1     = 3.620*inch;  // length along x (short side)
   // double xl_a1     = 5.456*inch;  // length along x (long side) 
   // G4Trap *colSolid_B = new G4Trap("colSolid_B",zl_a1,yl_a1,xl_a1,xs_a1);

   // need to use a general trapezoid since there's so many dimensions changing 
   // from Geant4 for producing a right angle trapezoid 
   // 00222 G4Trap::G4Trap( const G4String& pName,
   // 00223                       G4double pZ,
   // 00224                       G4double pY,
   // 00225                       G4double pX, G4double pLTX )
   // 00238   fDz = 0.5*pZ ;
   // 00239   fTthetaCphi = 0 ;
   // 00240   fTthetaSphi = 0 ;
   // 00241 
   // 00242   fDy1 = 0.5*pY;
   // 00243   fDx1 = 0.5*pX ;
   // 00244   fDx2 = 0.5*pLTX;
   // 00245   fTalpha1 =  0.5*(pLTX - pX)/pY;
   // 00246 
   // 00247   fDy2 = fDy1 ;
   // 00248   fDx3 = fDx1;
   // 00249   fDx4 = fDx2 ;
   // 00250   fTalpha2 = fTalpha1;
  
   // general constructor  
   // 00084 G4Trap::G4Trap( const G4String& pName,
   // 00085                       G4double pDz,
   // 00086                       G4double pTheta, G4double pPhi,
   // 00087                       G4double pDy1, G4double pDx1, G4double pDx2,
   // 00088                       G4double pAlp1,
   // 00089                       G4double pDy2, G4double pDx3, G4double pDx4,
   // 00090                       G4double pAlp2)
   // 00106   fDz=pDz;
   // 00107   fTthetaCphi=std::tan(pTheta)*std::cos(pPhi);
   // 00108   fTthetaSphi=std::tan(pTheta)*std::sin(pPhi);
   // 00109       
   // 00110   fDy1=pDy1;
   // 00111   fDx1=pDx1;
   // 00112   fDx2=pDx2;
   // 00113   fTalpha1=std::tan(pAlp1);
   // 00114      
   // 00115   fDy2=pDy2;
   // 00116   fDx3=pDx3;
   // 00117   fDx4=pDx4;
   // 00118   fTalpha2=std::tan(pAlp2);
   G4double z       = 3.397*inch;  // z length
   G4double y_mz    = 3.620*inch;  // y length at -z
   G4double y_pz    = 5.456*inch + (5.456*inch-y_mz);  // y length at +z should be 5.456 inch 
   G4double x_my_mz = 0.815*inch;  // x length at -y, -z  
   G4double x_py_mz = 0.815*inch;  // x length at +y, -z
   G4double x_my_pz = 1.531*inch;  // x length at -y, +z  
   G4double x_py_pz = 1.531*inch;  // x length at +y, +z
   G4double theta   = 0.*deg;  
   G4double phi     = 0.*deg; 
   G4double alpha1  = 0.5*(x_py_mz-x_my_mz)/y_pz; 
   G4double alpha2  = alpha1;
   G4Trap *trapB    = new G4Trap("trapB",z/2.,theta,phi,y_mz/2.,x_my_mz/2.,x_py_mz/2.,alpha1,y_pz/2.,x_my_pz/2.,x_py_pz/2.,alpha2); 

   // cut the bottom off 
   G4Box *bCut = new G4Box("bCut",2.*inch,2.*inch,2.*inch);
   G4ThreeVector P_bCut = G4ThreeVector(0,-y_mz-0.5*cm,0); 
   // G4LogicalVolume *bCut_LV = new G4LogicalVolume(bCut,Aluminum,"bCut_LV");

   // new G4PVPlacement(0,                         // rotation
   //                   P_bCut,                    // position 
   //                   bCut_LV,                   // logical volume   
   //                   "bCut_PHY",                // physical name 
   //                   logicMother,               // logical mother
   //                   false,                     // boolean? 
   //                   0,                         // copy no 
   //                   fCheckOverlaps);           // check overlaps

   G4SubtractionSolid *colSolid_B = new G4SubtractionSolid("colSolid_B",trapB,bCut,0,P_bCut);

   // // collimator B: regular trapezoid
   // double z_a1      = 3.397*inch;
   // double xs_a1     = 3.620*inch;  // upstream x size 
   // double xl_a1     = 5.456*inch;  // 3.189*inch;  // downstream x size 
   // double ys_a1     = 0.815*inch;  // thickness  
   // double yl_a1     = 1.531*inch;  // thickness  
   // G4Trd *colSolid_B = new G4Trd("colSolid_B",xs_a1/2.,xl_a1/2.,ys_a1/2.,yl_a1/2.,z_a1/2.);  

   // define materials and logical volume 
   G4LogicalVolume *col_B_LV = new G4LogicalVolume(colSolid_B,Aluminum,"col_B_LV"); 
   col_B_LV->SetVisAttributes(vis); 

   // placement in the Hall coordinate system 
   // position
   std::vector<G4double> POS;
   POS.push_back( 2.520*inch); POS.push_back(-0.55*inch); POS.push_back(z0 + 14.2*inch); 
   // POS.push_back(0); POS.push_back(0); POS.push_back(0); 
   G4ThreeVector P = G4ThreeVector(POS[0],POS[1],POS[2]); 
   // rotation 
   std::vector<G4double> RA; 
   // RA.push_back(0.*deg); RA.push_back(47.24*deg); RA.push_back(-90.*deg);
   RA.push_back(0.*deg); RA.push_back(-47.24*deg); RA.push_back(0.*deg);
   G4RotationMatrix *rm = new G4RotationMatrix();
   rm->rotateX(RA[0]); rm->rotateY(RA[1]); rm->rotateZ(RA[2]); 

   new G4PVPlacement(rm,                         // rotation
	             P,                          // position 
                     col_B_LV,                   // logical volume   
                     "col_B_PHY",                // physical name 
                     logicMother,                // logical mother
                     false,                      // boolean? 
                     0,                          // copy no 
                     fCheckOverlaps);            // check overlaps

   // now build the collimator base 
   G4double xb_len = 3.750*inch; 
   G4double yb_len = 0.765*inch; 
   G4double zb_len = 3.997*inch;
   G4Box *solidBase = new G4Box("solidBase",xb_len/2.,yb_len/2.,zb_len/2.);
   G4LogicalVolume *solidBase_LV = new G4LogicalVolume(solidBase,Aluminum,"solidBase_LV");
 
   // new G4PVPlacement(0,                          // rotation
   //                   G4ThreeVector(0,0,0),       // position 
   //                   solidBase_LV,               // logical volume   
   //                   "solidBase_PHY",            // physical name 
   //                   logicMother,                // logical mother
   //                   true,                       // boolean? 
   //                   0,                          // copy no 
   //                   fCheckOverlaps);            // check overlaps

   // base cut 1
   G4double xcc_len = 2.0*inch;  
   G4double ycc_len = 2.*inch;  
   G4double zcc_len = 3.700*inch; 
   G4Box *solidBase_cut1 = new G4Box("solidBase_cut1",xcc_len/2.,ycc_len/2.,zcc_len/2.); 
   G4LogicalVolume *solidBase_cut1_LV = new G4LogicalVolume(solidBase_cut1,Aluminum,"solidBase_cut1_LV");
   solidBase_cut1_LV->SetVisAttributes(vis);  

   G4RotationMatrix *rmc1 = new G4RotationMatrix();
   rmc1->rotateY(-36.04*deg); 

   G4double XX1 = 1.6*inch;  
   G4double YY1 = 0; 
   G4double ZZ1 = -1.05*inch;  // diagram says 1.008, but this looks nicer
   G4ThreeVector P_bc1 = G4ThreeVector(XX1,YY1,ZZ1);

   // new G4PVPlacement(rmc1,                       // rotation
   //                   P_bc1,                      // position 
   //                   solidBase_cut1_LV,          // logical volume   
   //                   "solidBase_cut1_PHY",       // physical name 
   //                   logicMother,                // logical mother
   //                   true,                       // boolean? 
   //                   0,                          // copy no 
   //                   fCheckOverlaps);            // check overlaps

   // base cut2  
   double zcc2     = 1.000*inch;  // length along z
   double ycc2     = 0.987*inch;  // length along y 
   double xcc2_s   = 2.439*inch;  // length along x (short side)
   double xcc2_l   = 3.352*inch;  // length along x (long side) 
   G4Trap *solidBase_cut2 = new G4Trap("solidBase_cut2",zcc2,ycc2,xcc2_l,xcc2_s);

   G4LogicalVolume *solidBase_cut2_LV = new G4LogicalVolume(solidBase_cut2,Aluminum,"solidBase_cut2_LV");
   solidBase_cut2_LV->SetVisAttributes(vis);  

   G4RotationMatrix *rmc2 = new G4RotationMatrix();
   rmc2->rotateY(-90*deg); 
   rmc2->rotateX(-30*deg); 

   G4double XX2 = -2.075*inch;  
   G4double YY2 =  0.495*inch; 
   G4double ZZ2 =  0.580*inch;
   // std::cout << XX2/inch << " " << YY2/inch << " " << ZZ2/inch << std::endl;
   G4ThreeVector P_bc2 = G4ThreeVector(XX2,YY2,ZZ2);

   // new G4PVPlacement(rmc2,                       // rotation
   //                   P_bc2,                      // position 
   //                   solidBase_cut2_LV,          // logical volume   
   //                   "solidBase_cut2_PHY",       // physical name 
   //                   logicMother,                // logical mother
   //                   false,                      // boolean? 
   //                   0,                          // copy no 
   //                   fCheckOverlaps);            // check overlaps

   // subtraction solid 
   // base - cut 1
   G4SubtractionSolid *colBase = new G4SubtractionSolid("colBase_1",solidBase,solidBase_cut1,rmc1,P_bc1);  
   // cut 2 
   colBase = new G4SubtractionSolid("colBase_12",colBase,solidBase_cut2,rmc2,P_bc2); 

   G4VisAttributes *visBase = new G4VisAttributes(); 
   visBase->SetColour( G4Colour::Blue() ); 
 
   G4LogicalVolume *colBase_B_LV = new G4LogicalVolume(colBase,Aluminum,"colBase_B_LV");
   colBase_B_LV->SetVisAttributes(visBase); 

   // placement 
   G4RotationMatrix *rmb = new G4RotationMatrix(); 
   rmb->rotateY(180.*deg);

   G4ThreeVector Pb = G4ThreeVector(2.875*inch,-2.815*inch,14.626*inch); 

   new G4PVPlacement(rmb,                        // rotation
                     Pb,                         // position 
                     colBase_B_LV,               // logical volume   
                     "colBase_B_PHY",            // physical name 
                     logicMother,                // logical mother
                     true,                       // boolean? 
                     0,                          // copy no 
                     fCheckOverlaps);            // check overlaps

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildCollimator_C(G4LogicalVolume *logicMother,G4double z0){
   // From drawings made by Sebastian Seeds (UConn) based on JT file
   // - Collimator C: CollimatorC_drawing.JPG
   // - Offsets and rotations: CollimatorC_xzoffset.JPG
   // Note: Collimators are on the LEFT side of the beam, next to the target  

   double inch   = 2.54*cm;
 
   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Red() ); 

   // define materials and logical volume 
   G4Material *Aluminum = G4Material::GetMaterial("G4_Al");

   // // collimator C: right-angle trapezoid  
   // double zl_a1     = 1.498*inch;  // length along z
   // double yl_a1     = 2.252*inch;  // length along y 
   // double xs_a1     = 3.800*inch;  // length along x (short side)
   // double xl_a1     = 4.858*inch;  // length along x (long side) 
   // G4Trap *colSolid_C = new G4Trap("colSolid_C",zl_a1,yl_a1,xl_a1,xs_a1);

   G4double z       = 2.756*inch;  // z length
   G4double y_mz    = 3.800*inch;  // y length at -z
   G4double y_pz    = 4.858*inch + (4.858*inch-y_mz);  // y length at +z should be 5.456 inch 
   G4double x_my_mz = 1.498*inch;  // x length at -y, -z  
   G4double x_py_mz = 1.498*inch;  // x length at +y, -z
   G4double x_my_pz = 2.252*inch;  // x length at -y, +z  
   G4double x_py_pz = 2.252*inch;  // x length at +y, +z
   G4double theta   = 0.*deg;  
   G4double phi     = 0.*deg; 
   G4double alpha1  = 0.5*(x_py_mz-x_my_mz)/y_pz; 
   G4double alpha2  = alpha1;
   G4Trap *trapC    = new G4Trap("trapC",z/2.,theta,phi,y_mz/2.,x_my_mz/2.,x_py_mz/2.,alpha1,y_pz/2.,x_my_pz/2.,x_py_pz/2.,alpha2); 

   // G4LogicalVolume *trapC_LV = new G4LogicalVolume(trapC,Aluminum,"trapC_LV");
   // trapC_LV->SetVisAttributes(vis); 

   // new G4PVPlacement(0,                         // rotation
   //                   G4ThreeVector(0,0,0),      // position 
   //                   trapC_LV,                   // logical volume   
   //                   "trapC_PHY",                // physical name 
   //                   logicMother,               // logical mother
   //                   false,                     // boolean? 
   //                   0,                         // copy no 
   //                   fCheckOverlaps);           // check overlaps

   // cut the bottom off 
   G4Box *cCut = new G4Box("cCut",2.*inch,2.*inch,2.*inch);
   G4ThreeVector P_cCut = G4ThreeVector(0,-y_mz-2*mm,0); 
   // G4LogicalVolume *cCut_LV = new G4LogicalVolume(cCut,Aluminum,"cCut_LV");

   // new G4PVPlacement(0,                         // rotation
   //                   P_cCut,                    // position 
   //                   cCut_LV,                   // logical volume   
   //                   "cCut_PHY",                // physical name 
   //                   logicMother,               // logical mother
   //                   false,                     // boolean? 
   //                   0,                         // copy no 
   //                   fCheckOverlaps);           // check overlaps

   G4SubtractionSolid *colSolid_C = new G4SubtractionSolid("colSolid_C",trapC,cCut,0,P_cCut);

   // define materials and logical volume 
   G4LogicalVolume *col_C_LV = new G4LogicalVolume(colSolid_C,Aluminum,"col_C_LV"); 
   col_C_LV->SetVisAttributes(vis); 

   // placement in the Hall coordinate system 
   // position
   std::vector<G4double> POS;
   POS.push_back(4.070*inch); POS.push_back(-0.35*inch); POS.push_back(z0 + 20.870*inch + 0.4*inch); 
   G4ThreeVector P = G4ThreeVector(POS[0],POS[1],POS[2]); 
   // rotation 
   std::vector<G4double> RA; 
   // RA.push_back(0.*deg); RA.push_back(50.41*deg); RA.push_back(-90.*deg);
   RA.push_back(0.*deg); RA.push_back(-50.41*deg); RA.push_back(0.*deg);
   G4RotationMatrix *rm = new G4RotationMatrix();
   rm->rotateX(RA[0]); rm->rotateY(RA[1]); rm->rotateZ(RA[2]); 

   new G4PVPlacement(rm,                         // rotation
	             P,                          // position 
                     col_C_LV,                   // logical volume   
                     "col_C_PHY",                // physical name 
                     logicMother,                // logical mother
                     false,                      // boolean? 
                     0,                          // copy no 
                     fCheckOverlaps);            // check overlaps

   // now build the collimator base 
   G4double xb_len = 3.625*inch; 
   G4double yb_len = 1.000*inch; 
   G4double zb_len = 3.250*inch;
   G4Box *solidBase = new G4Box("solidBase",xb_len/2.,yb_len/2.,zb_len/2.);
   G4LogicalVolume *solidBase_LV = new G4LogicalVolume(solidBase,Aluminum,"solidBase_LV");
 
   // new G4PVPlacement(0,                          // rotation
   //                   G4ThreeVector(0,0,0),       // position 
   //                   solidBase_LV,               // logical volume   
   //                   "solidBase_PHY",            // physical name 
   //                   logicMother,                // logical mother
   //                   false,                       // boolean? 
   //                   0,                          // copy no 
   //                   fCheckOverlaps);            // check overlaps

   // base cut 1
   G4double xcc_len = 2.420*inch;  
   G4double ycc_len = 4.0*inch;  
   G4double zcc_len = 2.420*inch; 
   G4Box *solidBase_cut1 = new G4Box("solidBase_cut1",xcc_len/2.,ycc_len/2.,zcc_len/2.); 
   G4LogicalVolume *solidBase_cut1_LV = new G4LogicalVolume(solidBase_cut1,Aluminum,"solidBase_cut1_LV");
   solidBase_cut1_LV->SetVisAttributes(vis);  

   G4RotationMatrix *rmc1 = new G4RotationMatrix();
   // rmc1->rotateY(-54.38*deg); 
   rmc1->rotateY(-35.62*deg); 

   G4double XX1 = 1.5*inch;  
   G4double YY1 = 0; 
   G4double ZZ1 = 0.2*inch + 1.5*inch + 0.216*inch;  
   G4ThreeVector P_bc1 = G4ThreeVector(XX1,YY1,ZZ1);

   // new G4PVPlacement(rmc1,                       // rotation
   //                   P_bc1,                      // position 
   //                   solidBase_cut1_LV,          // logical volume   
   //                   "solidBase_cut1_PHY",       // physical name 
   //                   logicMother,                // logical mother
   //                   false,                      // boolean? 
   //                   0,                          // copy no 
   //                   fCheckOverlaps);            // check overlaps

   // base cut2 
   G4double xcc2_len = 4.000*inch;  
   G4double ycc2_len = 0.500*inch;  
   G4double zcc2_len = 0.875*inch; 
   G4Box *solidBase_cut2 = new G4Box("solidBase_cut2",xcc2_len/2.,ycc2_len/2.,zcc2_len/2.); 
   G4LogicalVolume *solidBase_cut2_LV = new G4LogicalVolume(solidBase_cut2,Aluminum,"solidBase_cut2_LV");
   solidBase_cut2_LV->SetVisAttributes(vis);  

   G4RotationMatrix *rmc2 = new G4RotationMatrix();

   G4double XX2 = 0.*inch;  
   G4double YY2 = -yb_len/2. - 0.05*inch; 
   G4double ZZ2 = -zb_len/2. + zcc2_len/2.;
   // std::cout << XX2/inch << " " << YY2/inch << " " << ZZ2/inch << std::endl;
   G4ThreeVector P_bc2 = G4ThreeVector(XX2,YY2,ZZ2);

   // new G4PVPlacement(rmc2,                       // rotation
   //                   P_bc2,                      // position 
   //                   solidBase_cut2_LV,          // logical volume   
   //                   "solidBase_cut2_PHY",       // physical name 
   //                   logicMother,                // logical mother
   //                   false,                      // boolean? 
   //                   0,                          // copy no 
   //                   fCheckOverlaps);            // check overlaps

   // subtraction solid 
   // base - cut 1
   G4SubtractionSolid *colBase = new G4SubtractionSolid("colBase_1",solidBase,solidBase_cut1,rmc1,P_bc1);  
   // cut 2 
   colBase = new G4SubtractionSolid("colBase_12",colBase,solidBase_cut2,rmc2,P_bc2); 

   G4VisAttributes *visBase = new G4VisAttributes(); 
   visBase->SetColour( G4Colour::Blue() ); 
 
   G4LogicalVolume *colBase_C_LV = new G4LogicalVolume(colBase,Aluminum,"colBase_C_LV");
   colBase_C_LV->SetVisAttributes(visBase); 

   // placement 
   G4RotationMatrix *rmb = new G4RotationMatrix(); 
   rmb->rotateY(90.*deg);

   G4ThreeVector Pb = G4ThreeVector(3.812*inch,-2.764*inch,21.50*inch); 

   new G4PVPlacement(rmb,                        // rotation
                     Pb,                         // position 
                     colBase_C_LV,               // logical volume   
                     "colBase_C_PHY",            // physical name 
                     logicMother,                // logical mother
                     true,                       // boolean? 
                     0,                          // copy no 
                     fCheckOverlaps);            // check overlaps

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildGEnTarget_CollimatorTable(G4LogicalVolume *motherLog,G4double z0){
   // Collimator table on which the collimators sit
   // Drawings: C4.JPG--C9.JPG

   G4double inch = 2.54*cm;

   // C4: table platform 
   G4double x_len =  3.875*inch;
   G4double y_len =  0.500*inch;
   G4double z_len = 35.625*inch;
   G4Box *solidC4 = new G4Box("solidC4",x_len/2.,y_len/2.,z_len/2.);

   // C4 cutaway (angled cut 1) 
   G4double xca1_len = 2.000*inch;
   G4double yca1_len = 2.000*inch;
   G4double zca1_len = 3.143*inch;
   G4Box *solidC4_cut_a1 = new G4Box("solidC4_cut_a1",xca1_len/2.,yca1_len/2.,zca1_len/2.);

   G4RotationMatrix *rmca1 = new G4RotationMatrix();
   std::vector<G4double> RA;
   RA.push_back(0*deg); RA.push_back(-25*deg); RA.push_back(0);  
   rmca1->rotateX(RA[0]); rmca1->rotateY(RA[1]); rmca1->rotateZ(RA[2]);

   // position of cut (derived from drawing C4.JPG) 
   G4ThreeVector P_c4c_a1 = G4ThreeVector(-2.165*inch,0*inch,-12.509*inch);

   // C4 cutaway [long stretch] 
   G4double cutLen_x = 1.375*inch;
   G4double xc_len =  2.000*inch;
   G4double yc_len =  2.000*inch;
   G4double zc_len = 11.028*inch;
   G4Box *solidC4_cut = new G4Box("solidC4_cut",xc_len/2.,yc_len/2.,zc_len/2.);
   // position of cut (derived from drawing C4.JPG) 
   G4double xcp = (-1.)*(x_len/2. + xc_len/2. - cutLen_x);
   G4double ycp = 0.;
   G4double zcp = -6.0245*inch;
   G4ThreeVector P_c4c = G4ThreeVector(xcp,ycp,zcp);

   // C4 cutaway (angled cut 2) 
   G4double xca2_len = 2.000*inch;
   G4double yca2_len = 2.000*inch;
   G4double zca2_len = 3.143*inch;
   G4Box *solidC4_cut_a2 = new G4Box("solidC4_cut_a2",xca2_len/2.,yca2_len/2.,zca2_len/2.);

   std::vector<G4double> RA2;
   RA2.push_back(0*deg); RA2.push_back(25*deg); RA2.push_back(0);  
   G4RotationMatrix *rmca2 = new G4RotationMatrix();
   rmca2->rotateX(RA2[0]); rmca2->rotateY(RA2[1]); rmca2->rotateZ(RA2[2]);

   // position of cut (derived from drawing C4.JPG) 
   G4ThreeVector P_c4c_a2 = G4ThreeVector(-2.15*inch,0*inch,0.914*inch-0.42*inch);

   // C4 cutaway (angled cut 3) 
   G4double xca3_len = 2.000*inch;
   G4double yca3_len = 2.000*inch;
   G4double zca3_len = 1.150*inch;
   G4Box *solidC4_cut_a3 = new G4Box("solidC4_cut_a3",xca3_len/2.,yca3_len/2.,zca3_len/2.);

   std::vector<G4double> RA3;
   RA3.push_back(0*deg); RA3.push_back(-35*deg); RA3.push_back(0);  
   G4RotationMatrix *rmca3 = new G4RotationMatrix();
   rmca3->rotateX(RA3[0]); rmca3->rotateY(RA3[1]); rmca3->rotateZ(RA3[2]);

   // position of cut (derived from drawing C4.JPG) 
   G4ThreeVector P_c4c_a3 = G4ThreeVector(-2.40*inch,0*inch,9.0*inch);

   // C4 cutaway [second long stretch] 
   G4double cutLen_x2 = 0.750*inch;
   G4double xc2_len = 2.000*inch;
   G4double yc2_len = 2.000*inch;
   G4double zc2_len = 7.789*inch + 1.25*inch;
   G4Box *solidC4_cut2 = new G4Box("solidC4_cut2",xc2_len/2.,yc2_len/2.,zc2_len/2.);
   // position of cut (derived from drawing C4.JPG) 
   G4double xcp2 = (-1.)*(x_len/2. + xc2_len/2. - cutLen_x2);
   G4double ycp2 = 0.;
   G4double zcp2 = 13.668*inch - 0.25*inch;
   G4ThreeVector P_c4c_2 = G4ThreeVector(xcp2,ycp2,zcp2);

   // subtraction solid
   // base plate - angled cut 
   G4SubtractionSolid *C4c = new G4SubtractionSolid("C4c",solidC4,solidC4_cut_a1,rmca1,P_c4c_a1);
   // subtract first long stretch 
   C4c = new G4SubtractionSolid("C4c",C4c,solidC4_cut,0,P_c4c);
   // subtract second angled cut  
   C4c = new G4SubtractionSolid("C4c",C4c,solidC4_cut_a2,rmca2,P_c4c_a2);
   // subtract third angled cut  
   C4c = new G4SubtractionSolid("C4c",C4c,solidC4_cut_a3,rmca3,P_c4c_a3);
   // subtract third long stretch 
   C4c = new G4SubtractionSolid("C4c",C4c,solidC4_cut2,0,P_c4c_2);

   // C5 
   G4double x5_len = 0.250*inch; 
   G4double y5_len = 1.250*inch; 
   G4double z5_len = 35.625*inch;
   G4Box *solidC5 = new G4Box("solidC5",x5_len/2.,y5_len/2.,z5_len/2.); 
   G4ThreeVector P5 = G4ThreeVector(x_len/2.+x5_len/2.,-0.188*inch,0.);

   // C6--C9 
   // base
   G4double xb_len = 1.250*inch;  
   G4double yb_len = 7.946*inch;  
   G4double zb_len = 0.500*inch;  
   G4Box *solidBase = new G4Box("solidBase",xb_len/2.,yb_len/2.,zb_len/2.);
   // small blocks 
   G4double xbk_len = 2.250*inch;  
   G4double ybk_len = 0.500*inch;  
   G4double zbk_len = 0.500*inch;  
   G4Box *solidBK = new G4Box("solidBK",xbk_len/2.,ybk_len/2.,zbk_len/2.);
   G4ThreeVector P_bk = G4ThreeVector(0,yb_len/2.+ybk_len/2.,0);

   // union solid to build the leg 
   G4UnionSolid *leg = new G4UnionSolid("leg_bbk",solidBase,solidBK,0,P_bk);
   // other leg 
   P_bk = G4ThreeVector(0,-yb_len/2.-ybk_len/2.,0);
   leg = new G4UnionSolid("solidLeg",leg,solidBK,0,P_bk);

   // union solid for all components 
   G4UnionSolid *solidTable = new G4UnionSolid("solidTable_noLegs",C4c,solidC5,0,P5);
   // add C6 leg 
   G4double z0_l = -5.3125*inch;   // same as offset at the very end 
   G4double x6 = 0;
   G4double y6 = (-1.)*(y_len/2. + ybk_len + yb_len/2.);  
   G4double z6 = z0_l + -9.244*inch + ybk_len/2.; 
   G4ThreeVector P6 = G4ThreeVector(x6,y6,z6);
   solidTable = new G4UnionSolid("solidTable_c6",solidTable,leg,0,P6); 
   // add C7 leg  
   G4double x7 = 0;
   G4double y7 = y6;  
   G4double z7 = z0_l + 0*inch; 
   G4ThreeVector P7 = G4ThreeVector(x7,y7,z7);
   solidTable = new G4UnionSolid("solidTable_c7",solidTable,leg,0,P7); 
   // add C8 leg  
   G4double x8 = 0;
   G4double y8 = y6;  
   G4double z8 = z0_l + 13.381*inch + ybk_len/2.; 
   G4ThreeVector P8 = G4ThreeVector(x8,y8,z8);
   solidTable = new G4UnionSolid("solidTable_c8",solidTable,leg,0,P8); 
   // add C9 leg  
   G4double x9 = 0;
   G4double y9 = y6;  
   G4double z9 = z0_l + 20.256*inch + ybk_len/2.; 
   G4ThreeVector P9 = G4ThreeVector(x9,y9,z9);
   solidTable = new G4UnionSolid("solidTable_c9",solidTable,leg,0,P9); 

   G4Material *Aluminum = G4Material::GetMaterial("G4_Al");
   G4LogicalVolume *table_LV = new G4LogicalVolume(solidTable,Aluminum,"logicGEnTarget_colTable");

   // placement 
   G4double x_offset =  2.812*inch; 
   G4double y_offset = -3.514*inch;  
   G4double z_offset =  5.3125*inch; // taking y_len/2 - 12.5*inch
   G4ThreeVector P = G4ThreeVector(x_offset,y_offset,z0+z_offset);

   bool checkOverlaps = true;

   new G4PVPlacement(0,                          // rotation
                     P,                          // position 
                     table_LV,                   // logical volume   
                     "physGEnTarget_colTable",   // physical name 
                     motherLog,                  // logical mother
                     true,                       // boolean? 
                     0,                          // copy no 
                     checkOverlaps);             // check overlaps

    // // test placement 
    // G4LogicalVolume *boxCut_LV = new G4LogicalVolume(solidC4_cut_a1,Aluminum,"logicCut"); 
    // G4VisAttributes *vis = new G4VisAttributes(); 
    // vis->SetColor( G4Colour::Red() ); 
    // boxCut_LV->SetVisAttributes(vis); 

    // // std::cout << (-PP[0]+0.8*inch)/inch << " " << (PP[2]+0.8*inch)/inch << std::endl;

    // new G4PVPlacement(rmca1,                               // rotation
    //                   P_c4c_a1,                            // position 
    //                   boxCut_LV,                           // logical volume   
    //                   "physBoxCut",                        // physical name 
    //                   motherLog,                           // logical mother
    //                   false,                               // boolean? 
    //                   0,                                   // copy no 
    //                   checkOverlaps);                      // check overlaps

    // G4LogicalVolume *boxCut_LV2 = new G4LogicalVolume(solidC4_cut_a2,Aluminum,"logicCut"); 
    // boxCut_LV2->SetVisAttributes(vis); 

    // new G4PVPlacement(rmca2,                               // rotation
    //                   P_c4c_a2,                            // position 
    //                   boxCut_LV2,                           // logical volume   
    //                   "physBoxCut2",                       // physical name 
    //                   motherLog,                           // logical mother
    //                   false,                               // boolean? 
    //                   0,                                   // copy no 
    //                   checkOverlaps);                      // check overlaps

    // G4LogicalVolume *boxCut_LV3 = new G4LogicalVolume(solidC4_cut_a3,Aluminum,"logicCut3"); 
    // boxCut_LV3->SetVisAttributes(vis); 

    // new G4PVPlacement(rmca3,                               // rotation
    //                   P_c4c_a3,                            // position 
    //                   boxCut_LV3,                           // logical volume   
    //                   "physBoxCut3",                       // physical name 
    //                   motherLog,                           // logical mother
    //                   false,                               // boolean? 
    //                   0,                                   // copy no 
    //                   checkOverlaps);                      // check overlaps

    // G4LogicalVolume *boxCut_LV4 = new G4LogicalVolume(solidC4_cut2,Aluminum,"logicCut4"); 
    // boxCut_LV4->SetVisAttributes(vis); 

    // new G4PVPlacement(0,                               // rotation
    //                   P_c4c_2,                            // position 
    //                   boxCut_LV4,                           // logical volume   
    //                   "physBoxCut4",                       // physical name 
    //                   motherLog,                           // logical mother
    //                   false,                               // boolean? 
    //                   0,                                   // copy no 
    //                   checkOverlaps);                      // check overlaps

    // G4double zz = (-1.)*(35.625/2. - 4 - 2.84 - 11.028 - 2.84 - 5.033)*inch; 
    // DrawAxis(motherLog,'y',zz); 

}
