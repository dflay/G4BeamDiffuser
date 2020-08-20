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
/// \file BDEventAction.hh
/// \brief Definition of the BDEventAction class

#ifndef BDEventAction_h
#define BDEventAction_h 1

#include "G4UserEventAction.hh"

#include "BDHit.hh"
#include "BDoutput.hh"
#include "BDIO.hh"

#include "globals.hh"

#include <cstdlib> 
#include <set> 
#include <map>
#include <iterator>  

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in the hits collections.

class BDEventAction : public G4UserEventAction
{
   public:
      BDEventAction();
      virtual ~BDEventAction();

      virtual void  BeginOfEventAction(const G4Event* event);
      virtual void    EndOfEventAction(const G4Event* event);

      void SetIO(BDIO *io){ fIO = io; }

   private:
      // methods
      BDHitsCollection* GetHitsCollection(G4int hcID,const G4Event* event) const; 

      void PrintEventStatistics(G4double Edep,G4double TrackLength,G4int layer) const;
      void FillBDData(const G4Event *evt,BDHitsCollection *hc,BDoutput &out); 

      // data members                  
      BDIO *fIO;  
 
//      G4int fAbsHCID;
//      G4int fGapHCID;
      G4int fDiffHCID; 
};

#endif

