#include "BDIO.hh"
//______________________________________________________________________________
BDIO::BDIO(const char *fn){
   fFile = nullptr; 
   fTree = nullptr; 
   strcpy(fFilename,fn); 
}
//______________________________________________________________________________
BDIO::~BDIO(){
   if(fTree) { delete fTree; }
   fTree = nullptr;

   if(fFile) { delete fFile; }
   fFile = nullptr;
}
//______________________________________________________________________________
void BDIO::Initialize(){
   // if a file exists, close it and delete the pointer 
   if(fFile){
      fFile->Close();
      delete fFile;
      fFile = nullptr; 
   }

   // create a new file 
   fFile = new TFile(fFilename,"RECREATE");

   // delete existing tree if necessary
   if(fTree){ 
      delete fTree; 
      fTree = nullptr; 
   }

   // create tree
   fTree = new TTree("T","BD Test");

   // define branches 
   BranchBD(); 
}
//______________________________________________________________________________
void BDIO::FillTree(){
   if(!fTree){ 
      fprintf(stderr, "Error %s: %s line %d - Trying to fill non-existant tree\n", __PRETTY_FUNCTION__, __FILE__, __LINE__ );
      return;
   }
   fTree->Fill();
}
//______________________________________________________________________________
void BDIO::Write(){
   // check pointers
   assert(fFile);
   std::cout << "[BDIO]: FILE IS OK!" << std::endl;
   assert(fTree);
   std::cout << "[BDIO]: TREE IS OK!" << std::endl;
   if( !fFile->IsOpen() ){
      G4cerr << "ERROR: " << __FILE__ << " line " << __LINE__ << ": TFile not open" << G4endl;
      exit(1);
   }
   // cd to top level of file and write the tree
   fFile->cd();
   fTree->Write("T", TObject::kOverwrite);
}
//______________________________________________________________________________
void BDIO::CloseFile(){
   // reset the branch addressing, delete pointer
   fTree->ResetBranchAddresses();
   delete fTree;
   fTree = nullptr;
   // close file and delete pointer
   // fFile->Map();
   fFile->Close();
   delete fFile;
   fFile = nullptr;
}
//______________________________________________________________________________
void BDIO::BranchBD(G4String SDname){
   // create the branches for the BeamDiffuser 
   TString branch_name;
   TString branch_prefix = SDname.data();
   branch_prefix.ReplaceAll("/",".");
   // define branches
   fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(BDdata[SDname].nhits) );
   fTree->Branch( branch_name.Format("%s.hit.plane", branch_prefix.Data() ), &(BDdata[SDname].plane) );
   fTree->Branch( branch_name.Format("%s.hit.x"    , branch_prefix.Data() ), &(BDdata[SDname].x)     );
   fTree->Branch( branch_name.Format("%s.hit.y"    , branch_prefix.Data() ), &(BDdata[SDname].y)     );
   fTree->Branch( branch_name.Format("%s.hit.z"    , branch_prefix.Data() ), &(BDdata[SDname].z)     );
   fTree->Branch( branch_name.Format("%s.hit.t"    , branch_prefix.Data() ), &(BDdata[SDname].t)     );
   fTree->Branch( branch_name.Format("%s.hit.xg"   , branch_prefix.Data() ), &(BDdata[SDname].xg)    );
   fTree->Branch( branch_name.Format("%s.hit.yg"   , branch_prefix.Data() ), &(BDdata[SDname].yg)    );
   fTree->Branch( branch_name.Format("%s.hit.zg"   , branch_prefix.Data() ), &(BDdata[SDname].zg)    );
   fTree->Branch( branch_name.Format("%s.hit.trid" , branch_prefix.Data() ), &(BDdata[SDname].trid)  );
   fTree->Branch( branch_name.Format("%s.hit.mid"  , branch_prefix.Data() ), &(BDdata[SDname].mid)   );
   fTree->Branch( branch_name.Format("%s.hit.pid"  , branch_prefix.Data() ), &(BDdata[SDname].pid)   );
   fTree->Branch( branch_name.Format("%s.hit.p"    , branch_prefix.Data() ), &(BDdata[SDname].p)     );
   fTree->Branch( branch_name.Format("%s.hit.edep" , branch_prefix.Data() ), &(BDdata[SDname].edep)  );
   fTree->Branch( branch_name.Format("%s.hit.beta" , branch_prefix.Data() ), &(BDdata[SDname].beta)  );
}
//______________________________________________________________________________
void BDIO::SetBDData(G4String SDname,BDoutput data){
   BDdata[SDname] = data;
}
