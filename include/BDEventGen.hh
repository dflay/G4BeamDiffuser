#ifndef BDEVENTGEN_HH
#define BDEVENTGEN_HH

#include "globals.hh"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4SystemOfUnits.hh"

class BDEventGen { 

   public: 
      BDEventGen();
      ~BDEventGen();

      void SetBeamEnergy(double E)    { fBeamE     = E; fBeamP = G4ThreeVector(0.,0.,E); }
      void SetBeamCurrent(double v)   { fBeamCur   = v;    }  
      void SetBeamRasterX(double v)   { fRasterX   = v;    } 
      void SetBeamRasterY(double v)   { fRasterY   = v;    }
      void SetBeamPointingX(double v) { fPointingX = v;    } 
      void SetBeamPointingY(double v) { fPointingY = v;    } 

      double GetBeamEnergy()    const { return fBeamE;     }
      double GetBeamCurrent()   const { return fBeamCur;   }
      double GetBeamRasterX()   const { return fRasterX;   } 
      double GetBeamRasterY()   const { return fRasterY;   } 
      double GetBeamPointingX() const { return fPointingX; } 
      double GetBeamPointingY() const { return fPointingY; }

      int GenerateEvent();

      G4ThreeVector GetVertex() const { return fVertex;    } 
      G4ThreeVector GetBeamP()  const { return fBeamP;     }  

   private:
      double fBeamE;         // beam energy 
      double fBeamCur;       // beam current  
      double fRasterX;       // beam raster x
      double fRasterY;       // beam raster y
      double fPointingX;     // beam pointing x 
      double fPointingY;     // beam pointing y

      G4ThreeVector fVertex; // generated vertex based on raster and pointing
      G4ThreeVector fBeamP;  // beam momentum  

}; 

#endif 
