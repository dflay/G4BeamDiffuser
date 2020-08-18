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
/// \file BDEventAction.cc
/// \brief Implementation of the BDEventAction class

#include "BDEventAction.hh"
#include "BDCalorimeterSD.hh"
#include "BDCalorHit.hh"
#include "BDAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include <iostream> 

//______________________________________________________________________________
BDEventAction::BDEventAction()
 : G4UserEventAction(),
   fAbsHCID(-1),
   fGapHCID(-1),
   fDiffHCID(-1)
{

}
//______________________________________________________________________________
BDEventAction::~BDEventAction()
{

}
//______________________________________________________________________________
BDCalorHitsCollection* 
BDEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<BDCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("BDEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    
//______________________________________________________________________________
void BDEventAction::PrintEventStatistics(
                              G4double Edep, G4double TrackLength,
                              G4int layer) const
{

  // DFlay modification 
  // turn off printing gap info, since we removed it from the code 

  // print event statistics
  char msg[200]; 
  sprintf(msg,"[BDEventAction]: Beam Diffuser: E (tot) = %.3lf MeV, Track Len (tot) = %.3lf mm, Layer = %d",
          Edep/CLHEP::MeV,TrackLength/CLHEP::mm,layer); 
  G4cout << msg << G4endl;
  // G4cout << "   Diffuser: total energy: " << std::setw(7) << G4BestUnit(Edep, "Energy")
  //        << "       total track length: " << std::setw(7) << G4BestUnit(TrackLength, "Length")
  //        << "             layer number: " << std::setw(7) << layer 
  //    << G4endl; 
     // << "        Gap: total energy: " 
     // << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     // << "       total track length: " 
     // << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     // << G4endl;
}
//______________________________________________________________________________
void BDEventAction::BeginOfEventAction(const G4Event* /*event*/)
{

}
//______________________________________________________________________________
void BDEventAction::EndOfEventAction(const G4Event* event)
{

  // FIXME: Make sure ALL HITS are written!
  //        Look at G4SBSIO and GESBSGEMoutput for guidance  

  // DFlay modification 
  // turn off printing gap info, since we removed it from the code 
  
  // Get hits collections IDs (only once)
  // if(fAbsHCID  == -1) fAbsHCID  = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitsCollection");
  // if(fGapHCID  == -1) fGapHCID  = G4SDManager::GetSDMpointer()->GetCollectionID("GapHitsCollection");
  if(fDiffHCID == -1) fDiffHCID = G4SDManager::GetSDMpointer()->GetCollectionID("DiffuserHitsCollection");

  // Get hits collections
  // auto absoHC = GetHitsCollection(fAbsHCID, event);
  // auto gapHC  = GetHitsCollection(fGapHCID, event);
  auto diffHC = GetHitsCollection(fDiffHCID, event);

  int NHits = diffHC->entries();

  // Get hit with total values
  // auto absoHit = (*absoHC)[absoHC->entries()-1];
  // auto gapHit  = (*gapHC)[gapHC->entries()-1];
  auto diffHit = (*diffHC)[diffHC->entries()-1];   // THIS IS THE LAST HIT! 
 
  // Print per event (modulo n)
  auto eventID     = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     
    // PrintEventStatistics(
    //   absoHit->GetEdep(), absoHit->GetTrackLength(),
    //   gapHit->GetEdep(), gapHit->GetTrackLength());
    // PrintEventStatistics(diffHit->GetEdep(),diffHit->GetTrackLength(),diffHit->GetLayer());
    diffHC->PrintAllHits(); 
 
  }  
  
  // Fill histograms, ntuple

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill ntuple in specific units 
  G4double EDEP = diffHit->GetEdep()/CLHEP::GeV; 
  G4double ETOT = diffHit->GetTotalEnergy()/CLHEP::GeV; 
  G4double TLEN = diffHit->GetTrackLength()/CLHEP::mm;
  G4int layer   = diffHit->GetLayer();  

  G4ThreeVector pos = diffHit->GetPos();
  G4double x = pos.getX()/CLHEP::mm; 
  G4double y = pos.getY()/CLHEP::mm; 
  G4double z = pos.getZ()/CLHEP::mm;

  G4ThreeVector labPos = diffHit->GetLabPos();
  G4double lx = labPos.getX()/CLHEP::mm; 
  G4double ly = labPos.getY()/CLHEP::mm; 
  G4double lz = labPos.getZ()/CLHEP::mm;

  G4ThreeVector mom = diffHit->GetMomentum();
  G4double px = mom.getX()/(CLHEP::GeV/CLHEP::c_light); 
  G4double py = mom.getY()/(CLHEP::GeV/CLHEP::c_light); 
  G4double pz = mom.getZ()/(CLHEP::GeV/CLHEP::c_light);
 
  // fill histograms
  analysisManager->FillH1(0,EDEP);
  // analysisManager->FillH1(1, gapHit->GetEdep());
  analysisManager->FillH1(1,TLEN);
  // analysisManager->FillH1(3, gapHit->GetTrackLength());
  
  // fill ntuple
  analysisManager->FillNtupleDColumn(0 ,EDEP );
  // analysisManager->FillNtupleDColumn(1, gapHit->GetEdep());
  analysisManager->FillNtupleDColumn(1 ,TLEN );
  analysisManager->FillNtupleDColumn(2 ,ETOT );
  // analysisManager->FillNtupleDColumn(3, gapHit->GetTrackLength());
  analysisManager->FillNtupleDColumn(3 ,x    ); 
  analysisManager->FillNtupleDColumn(4 ,y    ); 
  analysisManager->FillNtupleDColumn(5 ,z    );
  analysisManager->FillNtupleDColumn(6 ,lx   ); 
  analysisManager->FillNtupleDColumn(7 ,ly   ); 
  analysisManager->FillNtupleDColumn(8 ,lz   ); 
  analysisManager->FillNtupleDColumn(9 ,px   ); 
  analysisManager->FillNtupleDColumn(10,py   ); 
  analysisManager->FillNtupleDColumn(11,pz   ); 
  analysisManager->FillNtupleDColumn(12,layer); 
  analysisManager->AddNtupleRow();  
}  

