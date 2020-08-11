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
/// \file BDCalorHit.cc
/// \brief Implementation of the BDCalorHit class

#include "BDCalorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<BDCalorHit>* BDCalorHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BDCalorHit::BDCalorHit()
 : G4VHit(),
   fEdep(0.),
   fTrackLength(0.),
   fEtot(0.)
{
   fPos.setX(0);    fPos.setY(0);    fPos.setZ(0);
   fLabPos.setX(0); fLabPos.setY(0); fLabPos.setZ(0);
   fMom.setX(0);    fMom.setY(0);    fMom.setZ(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BDCalorHit::~BDCalorHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BDCalorHit::BDCalorHit(const BDCalorHit& right)
  : G4VHit()
{
  fEdep        = right.fEdep;
  fTrackLength = right.fTrackLength;
  fEtot        = right.fEtot;
  fPos         = right.fPos;
  fLabPos      = right.fLabPos;
  fMom         = right.fMom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const BDCalorHit& BDCalorHit::operator=(const BDCalorHit& right)
{
  fEdep        = right.fEdep;
  fTrackLength = right.fTrackLength;
  fEtot        = right.fEtot;
  fPos         = right.fPos;
  fLabPos      = right.fLabPos;
  fMom         = right.fMom;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool BDCalorHit::operator==(const BDCalorHit& right) const
{
  return ( this == &right ) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BDCalorHit::Print()
{
  G4cout
     << "Edep: " 
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << " track length: " 
     << std::setw(7) << G4BestUnit( fTrackLength,"Length")
     << " track x: " 
     << std::setw(7) << G4BestUnit( fPos.getX(),"Length") 
     << " track y: " 
     << std::setw(7) << G4BestUnit( fPos.getY(),"Length") 
     << " track z: " 
     << std::setw(7) << G4BestUnit( fPos.getZ(),"Length") 
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
