#ifndef BD_MESSENGER_HH
#define BD_MESSENGER_HH

// messenger class for custom macro options

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"

#ifdef __APPLE__
#include <unistd.h>
#endif 

// forward class declarations
class BDEventGen;
class G4UIdirectory; 
class G4UIcmdWithADoubleAndUnit; 

class BDMessenger : public G4UImessenger {

   public: 
      BDMessenger();
      ~BDMessenger();

      void SetEvGen(BDEventGen *f) { fEventGen = f;}
      void SetNewValue(G4UIcommand *cmd,G4String newValue);

    private:
      BDEventGen *fEventGen;

      // directory 
      G4UIdirectory *fDirectory; 

      // beam raster 
      G4UIcmdWithADoubleAndUnit *fRasterXcmd; 
      G4UIcmdWithADoubleAndUnit *fRasterYcmd;
      // beam pointing 
      G4UIcmdWithADoubleAndUnit *fBeamPointXcmd;  
      G4UIcmdWithADoubleAndUnit *fBeamPointYcmd; 
      // beam energy 
      G4UIcmdWithADoubleAndUnit *fBeamEcmd; 
      // beam energy stdev  
      G4UIcmdWithADoubleAndUnit *fBeamESIGcmd; 

      // beam angular alignment 
      G4UIcmdWithADoubleAndUnit *fBeamAngleXcmd;  
      G4UIcmdWithADoubleAndUnit *fBeamAngleYcmd; 
      G4UIcmdWithADoubleAndUnit *fBeamAngleZcmd;
      G4UIcmdWithADoubleAndUnit *fBeamAngleSpreadcmd;

}; 

#endif 
