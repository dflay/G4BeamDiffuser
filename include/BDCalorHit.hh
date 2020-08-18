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
/// \file BDCalorHit.hh
/// \brief Definition of the BDCalorHit class

#ifndef BDCalorHit_h
#define BDCalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class BDCalorHit : public G4VHit
{
  public:
    BDCalorHit();
    BDCalorHit(const BDCalorHit&);
    virtual ~BDCalorHit();

    // operators
    const BDCalorHit& operator=(const BDCalorHit&);
    G4bool operator==(const BDCalorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de, G4double dl);
    void SetTotalEnergy(G4double E); 
    void SetLayer(G4int i); 

    void SetPos(G4ThreeVector v);
    void SetLabPos(G4ThreeVector v);
    void SetMomentum(G4ThreeVector m);

    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;
    G4double GetTotalEnergy() const; 
    G4int    GetLayer() const; 

    G4ThreeVector GetPos() const; 
    G4ThreeVector GetLabPos() const; 
    G4ThreeVector GetMomentum() const; 
      
  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength; ///< Track length in the sensitive volume
    G4double fEtot;        ///< Total energy (at pre-step) 
    G4int fLayer;          ///< Layer number
  
    G4ThreeVector fPos;    ///< Local hit coordinate 
    G4ThreeVector fLabPos; ///< Global hit coordinate 
    G4ThreeVector fMom;    ///< Momentum 

};

using BDCalorHitsCollection = G4THitsCollection<BDCalorHit>;
extern G4ThreadLocal G4Allocator<BDCalorHit>* BDCalorHitAllocator;

//______________________________________________________________________________
inline void* BDCalorHit::operator new(size_t)
{
  if (!BDCalorHitAllocator) {
    BDCalorHitAllocator = new G4Allocator<BDCalorHit>;
  }
  void *hit;
  hit = (void *) BDCalorHitAllocator->MallocSingle();
  return hit;
}
//______________________________________________________________________________
inline void BDCalorHit::operator delete(void *hit)
{
  if (!BDCalorHitAllocator) {
    BDCalorHitAllocator = new G4Allocator<BDCalorHit>;
  }
  BDCalorHitAllocator->FreeSingle((BDCalorHit*) hit);
}
//______________________________________________________________________________
inline void BDCalorHit::Add(G4double de, G4double dl) {
  fEdep        += de; 
  fTrackLength += dl;
}
//______________________________________________________________________________
inline G4double BDCalorHit::GetEdep() const { 
  return fEdep; 
}
//______________________________________________________________________________
inline G4double BDCalorHit::GetTrackLength() const { 
  return fTrackLength; 
}
//______________________________________________________________________________
inline void BDCalorHit::SetTotalEnergy(G4double E){
   fEtot = E;
}
//______________________________________________________________________________
inline void BDCalorHit::SetPos(G4ThreeVector v){
   fPos = v;
}
//______________________________________________________________________________
inline void BDCalorHit::SetLabPos(G4ThreeVector v){
   fLabPos = v;
}
//______________________________________________________________________________
inline void BDCalorHit::SetMomentum(G4ThreeVector m){
   fMom = m;
}
//______________________________________________________________________________
inline void BDCalorHit::SetLayer(G4int i){
   fLayer = i; 
}
//______________________________________________________________________________
inline G4double BDCalorHit::GetTotalEnergy() const{
   return fEtot;
}
//______________________________________________________________________________
inline G4ThreeVector BDCalorHit::GetPos() const{ 
   return fPos;
}
//______________________________________________________________________________
inline G4ThreeVector BDCalorHit::GetLabPos() const{ 
   return fLabPos;
}
//______________________________________________________________________________
inline G4ThreeVector BDCalorHit::GetMomentum() const{ 
   return fMom;
}
//______________________________________________________________________________
inline G4int BDCalorHit::GetLayer() const{
   return fLayer;
}

#endif
