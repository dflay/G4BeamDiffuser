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

      void Print();
      void Initialize();  
      void SetBeamEnergy(double E)           { fBeamE_in  = E; fBeamP = G4ThreeVector(0.,0.,E); }
      void SetBeamEnergySigma(double sig)    { fBeamE_sig = sig;  }
      void SetBeamCurrent(double v)          { fBeamCur   = v;    }  
      void SetBeamRasterX(double v)          { fRasterX   = v;    } 
      void SetBeamRasterY(double v)          { fRasterY   = v;    }
      void SetBeamPointingX(double v)        { fPointingX = v;    } 
      void SetBeamPointingY(double v)        { fPointingY = v;    } 

      double GetBeamEnergy()           const { return fBeamE;     }
      double GetBeamCurrent()          const { return fBeamCur;   }
      double GetBeamRasterX()          const { return fRasterX;   } 
      double GetBeamRasterY()          const { return fRasterY;   } 
      double GetBeamPointingX()        const { return fPointingX; } 
      double GetBeamPointingY()        const { return fPointingY; }

      // beam angle
      void SetZOrigin(double v)              { fZ0 = v; }
      void SetBeamAngleX(double v)           { fBeamAngleX = v;    } 
      void SetBeamAngleY(double v)           { fBeamAngleY = v;    } 
      void SetBeamAngleZ(double v)           { fBeamAngleZ = v;    } 
      void SetBeamAngleSpread(double v)      { fBeamAngleSpread = v; } 
 
      double GetBeamAngleX()           const { return fBeamAngleX; } 
      double GetBeamAngleY()           const { return fBeamAngleY; }
      double GetBeamAngleZ()           const { return fBeamAngleZ; }

      int GenerateEvent();

      G4ThreeVector GetVertex() const { return fVertex;    } 
      G4ThreeVector GetBeamP()  const { return fBeamP;     }  

   private:
      double fBeamE_in;      // beam energy (input) 
      double fBeamE_sig;     // beam energy stdev 
      double fBeamE;         // beam energy (random generated)  
      double fBeamCur;       // beam current  
      double fRasterX;       // beam raster x
      double fRasterY;       // beam raster y
      double fPointingX;     // beam pointing x 
      double fPointingY;     // beam pointing y
      double fX;             // beam x (random generated) 
      double fY;             // beam y (random generated)
      double fZ0;            // beam origin (z)  
      bool fInit; 

      // beam angle 
      double fBeamAngleX; 
      double fBeamAngleY; 
      double fBeamAngleZ;
      double fBeamAngleSpread;
      // smear widths
      double fbax_sig; 
      double fbay_sig; 
      double fbaz_sig; 

      G4ThreeVector fVertex; // generated vertex based on raster and pointing
      G4ThreeVector fBeamP;  // beam momentum  

}; 

#endif 
