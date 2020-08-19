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
/// \file BDActionInitialization.cc
/// \brief Implementation of the BDActionInitialization class

#include "BDActionInitialization.hh"
#include "BDPrimaryGeneratorAction.hh"
#include "BDRunAction.hh"
#include "BDEventAction.hh"
#include "BDEventGen.hh"
#include "BDMessenger.hh"
#include "BDIO.hh"

//______________________________________________________________________________
BDActionInitialization::BDActionInitialization()
 : G4VUserActionInitialization()
{
   fMessenger = new BDMessenger();
   fIO        = new BDIO(); 
}
//______________________________________________________________________________
BDActionInitialization::~BDActionInitialization()
{
   delete fMessenger;
   delete fIO;
}
//______________________________________________________________________________
void BDActionInitialization::BuildForMaster() const
{
   // SetUserAction(new BDRunAction);
   G4UserRunAction* run_action = new BDRunAction();
   ( (BDRunAction *) run_action )->SetIO(fIO);
   SetUserAction(run_action);
}
//______________________________________________________________________________
void BDActionInitialization::Build() const
{
   // much simpler way, when there's no custom classes that need access 
   // to specific data 
   // SetUserAction(new BDPrimaryGeneratorAction);
   // SetUserAction(new BDRunAction);
   // SetUserAction(new BDEventAction);

   // For the messenger, need to get a pointer to the event generator that is 
   // a data member of PrimaryGeneratorAction 
   G4VUserPrimaryGeneratorAction *gen_action = new BDPrimaryGeneratorAction(); 
   fMessenger->SetEvGen( ((BDPrimaryGeneratorAction *)gen_action)->GetEvGen() ); 
   SetUserAction(gen_action);

   // RunAction: we need to set its IO pointer to the one located here 
   G4UserRunAction* run_action = new BDRunAction();
   ( (BDRunAction *) run_action )->SetIO(fIO);
   SetUserAction(run_action);

   // EventAction: we need to set its IO pointer to the one located here 
   G4UserEventAction *event_action = new BDEventAction(); 
   ( (BDEventAction *)event_action )->SetIO(fIO); 
   SetUserAction(event_action);
}  

