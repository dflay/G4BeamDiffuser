#ifndef BD_DIFFUSER_PARAMETERISATION_HH
#define BD_DIFFUSER_PARAMETERISATION_HH

// Beam Diffuser
// - parameterization of the Hall A or C diffuser device 
// - repeated volumes of varying size, separated by some distance 

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

namespace diffuser {
   enum colorIndex {
      kWhite   = 0, 
      kGray    = 1, 
      kGrey    = 2, 
      kBlack   = 3, 
      kBrown   = 4, 
      kRed     = 5, 
      kGreen   = 6, 
      kBlue    = 7, 
      kCyan    = 8, 
      kMagenta = 9, 
      kYellow  = 10  
   };
} 

class BDParameterisation: public G4VPVParameterisation {

   public:
      BDParameterisation(char Hall='A',G4ThreeVector r0=G4ThreeVector(0,0,0) );
      virtual ~BDParameterisation();

      // position and rotation 
      void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;

      // size of the object 
      virtual void ComputeDimensions(G4Tubs &plate,const G4int copyNo,const G4VPhysicalVolume *physVol) const;

      // shape of the object 
      // virtual G4VSolid *ComputeSolid(const G4int copyNo,G4VPhysicalVolume *physVol);

      // NOTE: Don't need this since the shape does not change 
      // material, sensitivity, visual attributes 
      // G4VTouchable should not be used for ordinary parameterisation
      // virtual G4Material *ComputeMaterial(const G4int copyNo,G4VPhysicalVolume *physVol,
      //                                     const G4VTouchable *parentTouch=0);
 
      // initialize private member arrays based on scalar members
      void InitParameters();      

      // setters
      void SetHall(char h)        { fHall       = h; }  
      void SetGap(double v)       { fGap        = v; }  
      void SetWidth(double v)     { fWidth      = v; }  
      void SetRadiusMin(double v) { fRadius_min = v; }
      void SetNLayers(int v)      { fNLayers    = v; }  

      // getters 
      char GetHall()            const { return fHall;       } 
      G4int GetNLayers()        const { return fNLayers;    }  
      G4double GetGap()         const { return fGap;        }  
      G4double GetWidth()       const { return fWidth;      }  
      G4double GetRadiusMin()   const { return fRadius_min; }  
      G4double GetRadiusMax()   const { return fRadius_max; }  

   private:
      G4ThreeVector fR0;      // origin of device relative to mother volume
      double fTotalThickness; // default z-offset to center the assembly in its logical volume 
      double fGap;            // separation of the plates
      double fWidth;          // width of a plate 
      double fRadius_min;     // inner radius of a plate 
      double fRadius_max;     // outer radius of a plate (derived from width and inner radius)  
      double *fThickness;     // plate thicknesses
      double *fStartPhi;      // start angles 
      double *fDeltaPhi;      // step angles
      int    *fColor;         // plate colors  
      int fNLayers;           // number of layers 
      char fHall;             // A or C 

};

#endif 
