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
}
//______________________________________________________________________________
BDMessenger::~BDMessenger(){
   delete fDirectory;
   delete fRasterXcmd;
   delete fRasterYcmd;
   delete fBeamPointXcmd;
   delete fBeamPointYcmd;
}
//______________________________________________________________________________
void BDMessenger::SetNewValue(G4UIcommand *cmd,G4String newValue){
   // default to zero 
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
}
