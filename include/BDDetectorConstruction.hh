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
/// \file BDDetectorConstruction.hh
/// \brief Definition of the BDDetectorConstruction class

#ifndef BDDetectorConstruction_h
#define BDDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In ConstructSDandField() sensitive detectors of BDCalorimeterSD type
/// are created and associated with the Absorber and Gap volumes.
/// In addition a transverse uniform magnetic field is defined 
/// via G4GlobalMagFieldMessenger class.

class BDDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    BDDetectorConstruction();
    virtual ~BDDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
     
  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    // utilities 
    void DrawAxis(G4LogicalVolume *logicMother,char axis='n',G4double z=0);
    void GetRotatedCoordinates(std::vector<G4double> RA,std::vector<G4double> P,std::vector<G4double> &PP); 

    // test geometries  
    void BuildHe3(G4LogicalVolume *logicMother);  
    void BuildGlassCell(G4LogicalVolume *logicMother);  
    void BuildBeamPipe(G4LogicalVolume *logicMother);
    void BuildIonChamber(G4LogicalVolume *logicMother); 

    // upstream collimators
    void BuildCollimators(G4LogicalVolume *logicMother,G4double z0=0);  
    void BuildCollimator_A(G4LogicalVolume *logicMother,G4double z0=0); 
    void BuildCollimator_B(G4LogicalVolume *logicMother,G4double z0=0); 
    void BuildCollimator_C(G4LogicalVolume *logicMother,G4double z0=0);

    void BuildGEnTarget_CollimatorTable(G4LogicalVolume *motherLog,G4double z0=0); 
   
    // beam dump
    void BuildBeamDump(G4LogicalVolume *logicMother,G4double z0=0);  // default location  
    void BuildBeamDump_Diffuser(G4LogicalVolume *logicMother,char Hall,G4double z0=0);
    void BuildBeamDump_ISOWallWeldment(G4LogicalVolume *logicMother,G4double z0=0); 
    void BuildBeamDump_UpstreamPipe(G4LogicalVolume *logicMother,G4double z0=0); 
    void BuildBeamDump_DownstreamPipe(G4LogicalVolume *logicMother,G4double z0=0);

    // beam exit 
    void MakeBeamExit_TargetToMidPipe(G4LogicalVolume *logicMother,G4double z0=0); 
  
    // data members
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; // magnetic field messenger

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    G4int   fNofLayers;     // number of layers
    G4int   fNDiffLayers;   // number of diffuser layers

    G4double fBDLength;     // length of Beam Diffuser (includes plate-to-plate spacing)    

};

#endif

