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

class BDIO { 

   public:
      BDIO(const char *fn="BD.root");
      ~BDIO();

      void SetFileName(const char *fn)   { strcpy(fFilename,fn); }
 
      void Initialize();
      void FillTree();
      void Write();
      void CloseFile();
     
      void SetBDData(G4String SDname,BDoutput data);
      void BranchBD(G4String SDname="BD");  

   private:
      TFile *fFile;
      TTree *fTree;

      char fFilename[255]; 

      std::map<G4String,BDoutput> BDdata;  

};  

#endif 
