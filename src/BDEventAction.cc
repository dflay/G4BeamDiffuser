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

#define  _L_UNIT CLHEP::m
#define  _E_UNIT CLHEP::GeV
#define  _T_UNIT CLHEP::ns 

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

  // int NHits = diffHC->entries();

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

  // fill out custom output structure 
  BDoutput bd;
  FillBDData(event,diffHC,bd); 
 
  fIO->SetBDData("BD",bd);  
 
}  
//______________________________________________________________________________
void BDEventAction::FillBDData(const G4Event *evt,BDCalorHitsCollection *hc,BDoutput &out){
   // fill the output class with hit data

   out.Clear();   // clear previous event data 

   int nstep=0,trackID=0,bdID=0;
   double w=0; 
   std::map< int,std::set<int> > tracks_layers;           // key = BeamDiffuser layer ID, value = set of unique tracks with edep in layer 
   std::map< int,std::map<int,int> > nsteps_track_layer;  // number of steps by track/layer  
   std::map< int,std::map<int,int> > pid;                 // particle type
   // std::map< int,std::map<int,int> > mid;                 // material/medium type (?)   
   std::map< int,std::map<int,double> > x,y,z,t,p,edep;   
   std::map< int,std::map<int,double> > xg,yg,zg;   
   // std::map< int,std::map<int,double> > beta;   

   // loop over all "hits" (i.e., individual tracking steps)
   int NHits = (int)hc->entries();  
   for(int i=0;i<NHits;i++){
      // get track ID and BeamDiffuser layer ID 
      trackID = (*hc)[i]->GetTrackID();
      bdID = (*hc)[i]->GetLayer();
      // now we examine the track
      std::pair<std::set<int>::iterator, bool> track = tracks_layers[bdID].insert(trackID);
      if( track.second ){
	 // new track in this layer, first step
	 nsteps_track_layer[bdID][trackID] = 1;
         // time of hit 
         t[bdID][trackID]    = (*hc)[i]->GetHitTime(); 
         // positional data (local coordinates of detector)  
	 x[bdID][trackID]    = (*hc)[i]->GetPos().x(); 
	 y[bdID][trackID]    = (*hc)[i]->GetPos().y(); 
	 z[bdID][trackID]    = (*hc)[i]->GetPos().z();
         // positional data (global or lab coordinates) 
	 xg[bdID][trackID]   = (*hc)[i]->GetLabPos().x(); 
	 yg[bdID][trackID]   = (*hc)[i]->GetLabPos().y(); 
	 zg[bdID][trackID]   = (*hc)[i]->GetLabPos().z();
         // energy and momentum 
         edep[bdID][trackID] = (*hc)[i]->GetEdep();
         p[bdID][trackID]    = (*hc)[i]->GetMom();  // momentum (magnitude) 
	 // beta[bdID][trackID] = (*hc)[i]->GetBeta();
         // PID 
	 pid[bdID][trackID]  = (*hc)[i]->GetPID();
	 // mid[bdID][trackID]  = (*hc)[i]->GetMID();
      }else{
	 // existing track in this layer, additional step; increment sums and averages
	 nstep = nsteps_track_layer[bdID][trackID];
	 w     = (double)nstep/( (double)(nstep+1) );
	 // the coordinates below represent local BD hit coordinates 
	 x[bdID][trackID] = x[bdID][trackID]*w +( (*hc)[i]->GetPos().x() )*(1.0-w);
	 y[bdID][trackID] = y[bdID][trackID]*w +( (*hc)[i]->GetPos().y() )*(1.0-w);
	 z[bdID][trackID] = z[bdID][trackID]*w +( (*hc)[i]->GetPos().z() )*(1.0-w);
         // the global coordinates 
	 xg[bdID][trackID] = xg[bdID][trackID]*w +( (*hc)[i]->GetLabPos().x() )*(1.0-w);
	 yg[bdID][trackID] = yg[bdID][trackID]*w +( (*hc)[i]->GetLabPos().y() )*(1.0-w);
	 zg[bdID][trackID] = zg[bdID][trackID]*w +( (*hc)[i]->GetLabPos().z() )*(1.0-w);
	 // for edep, we do the sum:
	 edep[bdID][trackID] += (*hc)[i]->GetEdep();
         // increment 
	 nsteps_track_layer[bdID][trackID]++;
      }
   } 

   bdID = 0; 
   trackID = 0;

   G4TrajectoryContainer *trajectorylist = evt->GetTrajectoryContainer(); // for particle history information

   // for particle history details.  don't utilize until we move this over to g4sbs
   // int MIDtemp=0,TIDtemp=0,PIDtemp=0,hitidx=0,nbouncetemp=0; 
   // std::set<int> TIDs_unique; //all unique track IDs involved in BD hits in this event (for filling particle history tree)

   // now accumulate data into output class 
   for(std::map<int,std::set<int> >::iterator hit=tracks_layers.begin(); hit!=tracks_layers.end(); hit++ ){
      std::set<int> tracklist = hit->second;
      bdID = hit->first;
      for(std::set<int>::iterator track=tracklist.begin(); track!=tracklist.end(); track++ ){
	 trackID = *track;
	 out.plane.push_back( bdID );
	 out.t.push_back( t[bdID][trackID]/_T_UNIT );
         // coordinates in detector system
	 out.x.push_back( (-y[bdID][trackID])/_L_UNIT );
	 out.y.push_back( (x[bdID][trackID])/_L_UNIT );
	 out.z.push_back( z[bdID][trackID]/_L_UNIT );
         // coordinates in the hall 
	 out.xg.push_back( xg[bdID][trackID]/_L_UNIT );
	 out.yg.push_back( yg[bdID][trackID]/_L_UNIT );
	 out.zg.push_back( zg[bdID][trackID]/_L_UNIT );
	 out.trid.push_back( trackID );
	 out.pid.push_back( pid[bdID][trackID] );
	 // out.mid.push_back( mid[bdID][trackID] );
	 out.p.push_back( p[bdID][trackID]/_E_UNIT );
	 // out.beta.push_back( beta[bdID][trackID]/_E_UNIT );
	 out.edep.push_back( edep[bdID][trackID]/_E_UNIT );
	 out.nhits++;
	 // if( trajectorylist ){ 
	 //    // fill Particle History, starting with the particle itself 
         //    // and working all the way back to primary particles:
	 //    MIDtemp = mid[gemID][trackID];
	 //    TIDtemp = trackID;
	 //    PIDtemp = pid[gemID][trackID];
	 //    hitidx = out.nhits;
	 //    nbouncetemp = 0;
	 //    do {
	 //       G4Trajectory *trajectory = (G4Trajectory*) (*trajectorylist)[TrajectoryIndex[TIDtemp]];
	 //       PIDtemp = trajectory->GetPDGEncoding();
	 //       MIDtemp = MotherTrackIDs[TIDtemp];
	 //       std::pair<set<int>::iterator, bool > newtrajectory = TIDs_unique.insert( TIDtemp );
	 //       if( newtrajectory.second ){ 
	 //          // this trajectory does not yet exist in the 
         //          // particle history of this detector for this event. Add it:
	 //          out.ParticleHistory.PID.push_back( PIDtemp );
	 //          out.ParticleHistory.MID.push_back( MIDtemp );
	 //          out.ParticleHistory.TID.push_back( TIDtemp );
	 //          // for hitindex: of course, this means that if a trajectory is 
         //          // involved in multiple hits in this detector, this variable 
         //          // will point to the first hit encountered only!
	 //          out.ParticleHistory.hitindex.push_back( hitidx ); 
	 //          out.ParticleHistory.nbounce.push_back( nbouncetemp );
	 //          out.ParticleHistory.vx.push_back( (trajectory->GetPoint(0)->GetPosition() ).x()/_L_UNIT );
	 //          out.ParticleHistory.vy.push_back( (trajectory->GetPoint(0)->GetPosition() ).y()/_L_UNIT );
	 //          out.ParticleHistory.vz.push_back( (trajectory->GetPoint(0)->GetPosition() ).z()/_L_UNIT );
	 //          out.ParticleHistory.px.push_back( (trajectory->GetInitialMomentum() ).x()/_E_UNIT );
	 //          out.ParticleHistory.py.push_back( (trajectory->GetInitialMomentum() ).y()/_E_UNIT );
	 //          out.ParticleHistory.pz.push_back( (trajectory->GetInitialMomentum() ).z()/_E_UNIT );
	 //          out.ParticleHistory.npart++;
	 //       }
	 //       TIDtemp = MIDtemp;
	 //       nbouncetemp++;
	 //    } while( MIDtemp!=0 );
	 // }
      }

   }

}

