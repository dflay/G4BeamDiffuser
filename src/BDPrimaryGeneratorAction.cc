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
/// \file BDPrimaryGeneratorAction.cc
/// \brief Implementation of the BDPrimaryGeneratorAction class

#include "BDPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "BDEventGen.hh"

//______________________________________________________________________________
BDPrimaryGeneratorAction::BDPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  auto particleDef = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particleDef);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(50.*MeV);

  fParticleMass = particleDef->GetPDGMass(); 

  fEventGen = new BDEventGen(); 

}
//______________________________________________________________________________
BDPrimaryGeneratorAction::~BDPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fEventGen;
}
//______________________________________________________________________________
void BDPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("BDPrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 
 
  // use our event generator 
  // - utilizes input for beam rastering and (mis)pointing 
  // - beam energy set from input file
  // - for some reason, sometimes the input macro values are NOT used.  
  //   to fix this, we have an Initialize() function that gives reasonable values 
  //   (that is, identical to the input macro).  this is not good practice, but it's fine for this test.  
  int rc = 1;
  int cntr=0; 
  do {
     rc = fEventGen->GenerateEvent();
  }while(rc!=0);  
 
  // implement beam parameters
  // - raster
  // - mispointing
  G4ThreeVector vert = fEventGen->GetVertex(); 
  vert.setZ(-worldZHalfLength);  // z is fixed 
  // std::cout << "[PrimaryGenerator]: Beam position x = " 
  //           << vert.getX()/mm << " mm, y = " 
  //           << vert.getY()/mm << " mm, z = "  
  //           << vert.getZ()/mm << " mm" << std::endl;
  fParticleGun->SetParticlePosition(vert);

  // beam energy (convert to kinetic energy first!) 
  // T = E - m 
  G4double E = fEventGen->GetBeamEnergy();
  G4double T = E - fParticleMass;
  // std::cout << "[PrimaryGenerator]: Beam E = " << E/GeV << " GeV, T = " << T/GeV << " GeV" << std::endl; 
  fParticleGun->SetParticleEnergy(T);  

  // std::cout << "-------------------------------------------------------" << std::endl; 
  
  // Set gun position
  // fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -worldZHalfLength));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

