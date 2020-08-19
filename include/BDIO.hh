#ifndef BEAM_DIFFUSER_IO_HH
#define BEAM_DIFFUSER_IO_HH

// custom input/output handler separate from the G4 default analyzerManager
// designed to mimic what is implemented in g4sbs to ease the integration process 

#include <cstdlib>
#include <map>
#include <iterator>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h" 
#include "TROOT.h"
#include "TObject.h"
#include "THashTable.h"
#include "TClonesArray.h"

#include "G4Run.hh"

#include "BDoutput.hh"

#define MAXHITDATA 2000

// typedef struct hit {
//   Double_t p[MAXHITDATA], edep[MAXHITDATA];
//   Double_t x[MAXHITDATA], y[MAXHITDATA], z[MAXHITDATA], t[MAXHITDATA];
//   Double_t xg[MAXHITDATA], yg[MAXHITDATA], zg[MAXHITDATA];
//   Int_t gid[MAXHITDATA],trkid[MAXHITDATA],trid[MAXHITDATA];
//   Int_t mid[MAXHITDATA], pid[MAXHITDATA];
//   Int_t ndata;
// } hit_t; 

class BDIO { 

   public:
      BDIO(const char *fn="BD.root");
      ~BDIO();

      void SetFileName(const char *fn)   { strcpy(fFilename,fn); }
      
      void SetBDData(G4String SDname,BDoutput data);
      void FillTree();
      void WriteTree();
      void InitializeTree();
      void BranchBD(G4String SDname="BD");  

   private:
      TFile *fFile;
      TTree *fTree;

      char fFilename[255]; 

      std::map<G4String,BDoutput> BDdata;  

};  

#endif 
