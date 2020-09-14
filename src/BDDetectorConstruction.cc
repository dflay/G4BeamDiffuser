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
#include "G4Cons.hh"
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
  
  BuildDiffuser(worldLV,'A'); 

  BuildBeamDump(worldLV);
 
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
void BDDetectorConstruction::BuildBeamDump(G4LogicalVolume *logicMother){
   // build the Hall A beam dump 

   BuildBeamDump_ISOWallWeldment(logicMother);
   BuildBeamDump_UpstreamPipe(logicMother);  
   BuildBeamDump_DownstreamPipe(logicMother);  

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildBeamDump_ISOWallWeldment(G4LogicalVolume *logicMother){
   // Hall A Beam Dump: ISO Wall Weldment 
   // Drawing: JL0015694, JL0015725 

   G4double inch  = 25.4*mm; 
   G4double x_len = 84.*inch;
   G4double y_len = 117.*inch; 
   G4double z_len = 4.*inch;

   auto material  = G4Material::GetMaterial("G4_Al"); // might not be aluminum... 

   // solid box 
   G4Box *solidBox = new G4Box("isoBox",x_len/2.,y_len/2.,z_len/2.);  

   // cut a circular hole 
   G4double r_min = 0.*mm; 
   G4double r_max = 16.*inch; 
   G4double len   = z_len + 5.*inch; // make sure it cuts through  
   G4double startPhi = 0.*deg; 
   G4double dPhi = 360.*deg; 
   G4Tubs *solidTube = new G4Tubs("isoTube",r_min,r_max,len/2.,startPhi,dPhi); 

   // cut the hole in the box
   G4double y_pos = 13.75*inch; 
   G4ThreeVector P = G4ThreeVector(0.*inch,13.75*inch,z_len/2.); 
   G4SubtractionSolid *solidWall = new G4SubtractionSolid("solidWall",solidBox,solidTube,0,P);
 
   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Red() );
   vis->SetForceWireframe(); 

   // logical volume
   G4LogicalVolume *isoWallLV = new G4LogicalVolume(solidWall,material,"isoWallLV");
   isoWallLV->SetVisAttributes(vis); 

   // placement 
   G4ThreeVector P_wall = G4ThreeVector(0.*m,-y_pos,2.25*m);
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
void BDDetectorConstruction::BuildBeamDump_UpstreamPipe(G4LogicalVolume *logicMother){
   // Hall A Beam Dump: Pipe upstream of ISO Weldment
   // Drawing: JL0009934-C-VAC SPOOL REGION UPPER LEVEL
   
   G4double inch         = 25.4*mm; 
   G4double TOTAL_LENGTH = 196.91*inch;
   G4double startPhi     = 0.*deg; 
   G4double dPhi         = 360.*deg; 

   // large conical tube [item 1]
   G4double delta      = 0.005*inch;  // FIXME: arbitrary!   
   G4double r_min1_lg  = 0.5*23.54*inch; 
   G4double r_max1_lg  = r_min1_lg + delta;   
   G4double r_min2_lg  = 0.5*37.54*inch; 
   G4double r_max2_lg  = r_min2_lg + delta;  
   G4double len_lg     =  9.95*inch; 
   G4Cons *solidConeLG = new G4Cons("solidConeLG",r_min1_lg,r_max1_lg,r_min2_lg,r_max2_lg,len_lg/2.,startPhi,dPhi); 
 
   // small conical tube [item 3] 
   G4double r_min1_sm  = 0.5*11.99*inch; 
   G4double r_max1_sm  = r_min1_sm + delta;   
   G4double r_min2_sm  = 0.5*23.49*inch; 
   G4double r_max2_sm  = r_min2_sm + delta;  
   G4double len_sm     =  6.13*inch; 
   G4Cons *solidConeSM = new G4Cons("solidConeSM",r_min1_sm,r_max1_sm,r_min2_sm,r_max2_sm,len_sm/2.,startPhi,dPhi); 

   // vacuum window tube [item 4] 
   G4double r_min_vwt = 0.5*12.12*inch; 
   G4double r_max_vwt = 0.5*12.50*inch; 
   G4double len_vwt   = 1.88*inch;
   G4Tubs *solidVWTube = new G4Tubs("solidVWTube",r_min_vwt,r_max_vwt,len_vwt/2.,startPhi,dPhi);

   // main tube [item 6, inferred]  
   G4double r_min    = 12.005*inch;   // FIXME: This is arbitrary! 
   G4double r_max    = 12.010*inch;   // FIXME: This is arbitrary! old value = 12.12*inch;   
   G4double len      = TOTAL_LENGTH - len_lg - len_sm - len_vwt; // derived number.   
   G4Tubs *solidTube = new G4Tubs("solidTube",r_min,r_max,len/2.,startPhi,dPhi);

   // union solid
   // - attach the main pipe and the large cone 
   G4double zz = 0.5*len + 0.5*len_lg;  
   G4ThreeVector P_lg = G4ThreeVector(0,0,zz);
   G4UnionSolid *upstrPipe = new G4UnionSolid("pipe_lg",solidTube,solidConeLG,0,P_lg); 
   // - attach the small cone    
   zz = 0.5*len + 0.5*len_sm;  
   G4ThreeVector P_sm = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("pipe_lg_sm",upstrPipe,solidConeSM,0,P_sm); 
   // - attach the vacuum window tube  
   zz = 0.5*len + len_sm + 0.5*len_vwt;   
   G4ThreeVector P_vwt = G4ThreeVector(0,0,-zz);  
   upstrPipe = new G4UnionSolid("upstrPipe",upstrPipe,solidVWTube,0,P_vwt); 
 
   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Green() ); 
   // vis->SetForceWireframe(); 
   
   // logical volume
   auto material  = G4Material::GetMaterial("G4_Al"); // might not be aluminum... 
   G4LogicalVolume *tubeLV = new G4LogicalVolume(upstrPipe,material,"beamDump_usPipe_LV");
   tubeLV->SetVisAttributes(vis);

   // placement
   G4double z_pos = 2.25*m - 0.5*TOTAL_LENGTH;  
   G4ThreeVector P = G4ThreeVector(0.*m,0,z_pos);
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

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildBeamDump_DownstreamPipe(G4LogicalVolume *logicMother){
   // Hall A Beam Dump: Pipe downstream of ISO Weldment
   // Drawing: JL0011756_27020E0145 [modified]

   G4double inch     = 25.4*mm; 
   G4double r_min    = 11.*inch;   // FIXME: This is arbitrary! 
   G4double r_max    = 12.12*inch;   
   G4double len      = 328.04*inch; 
   G4double startPhi = 0.*deg; 
   G4double dPhi     = 360.*deg; 
   G4Tubs *solidTube = new G4Tubs("solidTube",r_min,r_max,len/2.,startPhi,dPhi);
 
   // visualization
   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Yellow() ); 
   // vis->SetForceWireframe(); 
   
   // logical volume
   auto material  = G4Material::GetMaterial("G4_Al"); // might not be aluminum... 
   G4LogicalVolume *tubeLV = new G4LogicalVolume(solidTube,material,"beamDump_usPipe_LV");
   tubeLV->SetVisAttributes(vis); 

   // placement
   G4double z_pos = 2.5*m + 15.*cm + len/2.;  
   G4ThreeVector P = G4ThreeVector(0.*m,0,z_pos);
   new G4PVPlacement(0,                        // no rotation
	             P,                        // location in mother volume 
	             tubeLV,                   // its logical volume                         
	             "beamDump_usPipe_PHY",    // its name
	             logicMother,              // its mother  volume
	             false,                    // boolean operation? 
	             0,                        // copy number
	             fCheckOverlaps);          // checking overlaps    

}
//______________________________________________________________________________
void BDDetectorConstruction::BuildDiffuser(G4LogicalVolume *logicMother,char Hall){
 
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

   G4double r_min    = 0.;
   G4double r_max    = 25.*inch;
   G4double len      = 15.*cm; 
   G4double startPhi = 0.*deg; 
   G4double dPhi     = 360.*deg; 
   auto diffCaseS = new G4Tubs("diffCase",r_min,r_max,len/2.,startPhi,dPhi);

   auto diffCaseLV = new G4LogicalVolume(diffCaseS,defaultMaterial,"diffCase"); // its name
  
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
   r_min    = 17.67*inch;  
   r_max    = 24.*inch; 
   dPhi     = 38.*deg; 
   startPhi = 270.*deg - dPhi/2.; 
   G4double thk = 0.125*inch;
  
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
   new G4PVParameterised("Diffuser",plateLV,diffCaseLV,kZAxis,fNDiffLayers,plateParam); 

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
