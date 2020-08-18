#include "BDMessenger.hh"
#include "BDEventGen.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
// #include "G4UIcmdWithAnInteger.hh"
// #include "G4UIcmdWithAString.hh"
// #include "G4UIcmdWithADouble.hh"
// #include "G4UIcmdWithABool.hh"
// #include "G4UIcmdWith3Vector.hh"
//______________________________________________________________________________
BDMessenger::BDMessenger(){
   // create the top level directory 
   fDirectory = new G4UIdirectory("/b4/"); 
   fDirectory->SetGuidance("Commands specific to the beam"); 
   // raster
   // - horizontal (x) 
   fRasterXcmd = new G4UIcmdWithADoubleAndUnit("/b4/rasterx",this);
   fRasterXcmd->SetGuidance("Horizontal (x) raster size");
   fRasterXcmd->SetParameterName("rasterx",false); // must provide input 
   // - vertical (y)
   fRasterYcmd = new G4UIcmdWithADoubleAndUnit("/b4/rastery",this);
   fRasterYcmd->SetGuidance("Vertical (y) raster size");
   fRasterYcmd->SetParameterName("rastery",false); // must provide input
   // pointing 
   // - horizontal (x)  
   fBeamPointXcmd = new G4UIcmdWithADoubleAndUnit("/b4/beampointx",this); 
   fBeamPointXcmd->SetGuidance("Set beam pointing along the horizontal (x) direction"); 
   fBeamPointXcmd->SetParameterName("beampointx",false);  // must provide input 
   // - vertical (y)
   fBeamPointYcmd = new G4UIcmdWithADoubleAndUnit("/b4/beampointy",this); 
   fBeamPointYcmd->SetGuidance("Set beam pointing along the vertical (y) direction"); 
   fBeamPointYcmd->SetParameterName("beampointy",false);  // must provide input 
   // beam energy 
   fBeamEcmd = new G4UIcmdWithADoubleAndUnit("/b4/beamE",this); 
   fBeamEcmd->SetGuidance("Set beam energy");
   fBeamEcmd->SetParameterName("beamE",false); 
   // beam energy stdev  
   fBeamESIGcmd = new G4UIcmdWithADoubleAndUnit("/b4/beamEsig",this); 
   fBeamESIGcmd->SetGuidance("Set beam energy standard deviation");
   fBeamESIGcmd->SetParameterName("beamEsig",false); 
}
//______________________________________________________________________________
BDMessenger::~BDMessenger(){
   delete fDirectory;
   delete fRasterXcmd;
   delete fRasterYcmd;
   delete fBeamPointXcmd;
   delete fBeamPointYcmd;
   delete fBeamEcmd; 
   delete fBeamESIGcmd; 
}
//______________________________________________________________________________
void BDMessenger::SetNewValue(G4UIcommand *cmd,G4String newValue){
   // default to zero 
   G4double be=0,bes=0; 
   G4double rx=0,ry=0;
   G4double px=0,py=0;
   // process commands
   if(cmd==fRasterXcmd){
      rx = fRasterXcmd->GetNewDoubleValue(newValue); 
      fEventGen->SetBeamRasterX(rx);  
   } 
   if(cmd==fRasterYcmd){
      ry = fRasterYcmd->GetNewDoubleValue(newValue);
      fEventGen->SetBeamRasterY(ry); 
   } 
   if(cmd==fBeamPointXcmd){
      px = fBeamPointXcmd->GetNewDoubleValue(newValue); 
      fEventGen->SetBeamPointingX(px); 
   } 
   if(cmd==fBeamPointYcmd){
      py = fBeamPointYcmd->GetNewDoubleValue(newValue);
      fEventGen->SetBeamPointingY(py); 
   } 
   if(cmd==fBeamEcmd){
      be = fBeamEcmd->GetNewDoubleValue(newValue); 
      fEventGen->SetBeamEnergy(be); 
   }
   if(cmd==fBeamESIGcmd){
      bes = fBeamESIGcmd->GetNewDoubleValue(newValue); 
      fEventGen->SetBeamEnergySigma(bes); 
   }
}
