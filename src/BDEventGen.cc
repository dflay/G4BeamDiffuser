#include "BDEventGen.hh"
//______________________________________________________________________________
BDEventGen::BDEventGen(){
   fInit       = false;
   fBeamE      = 0.;  // generated beam energy
   fBeamE_in   = 0.;  // input beam energy
   fBeamE_sig  = 0.;
   fBeamCur    = 0.;
   fRasterX    = 0.;
   fRasterY    = 0.;
   fPointingX  = 0.;
   fPointingY  = 0.;
   fBeamAngleX = 0.;
   fBeamAngleY = 0.;
   fBeamAngleZ = 0.;
   fBeamAngleSpread = 0.1*CLHEP::perCent; 
   fbax_sig    = 1.*CLHEP::mrad;
   fbay_sig    = 1.*CLHEP::mrad;
   fbaz_sig    = 1.*CLHEP::mrad;
   fX          = 0.;
   fY          = 0.;
   fVertex     = G4ThreeVector(0,0,0);
   fBeamP      = G4ThreeVector(0,0,0);  
}
//______________________________________________________________________________
BDEventGen::~BDEventGen(){

}
//______________________________________________________________________________
void BDEventGen::Initialize(){
   // set to some reasonable values
   // for some reason the input values from the macro aren't used all the time... 
   fBeamE_in   = 11.*CLHEP::GeV; 
   fBeamE_sig  = 20.*CLHEP::MeV; 
   fRasterX    = 1.*CLHEP::mm; 
   fRasterY    = 1.*CLHEP::mm; 
   fPointingX  = 0.*CLHEP::mm; 
   fPointingY  = 0.*CLHEP::mm;
   fBeamAngleX = 0.*CLHEP::rad;
   fBeamAngleY = 0.*CLHEP::rad;
   fBeamAngleZ = 0.*CLHEP::rad;
   std::cout << "[BDEventGen]: Now initialized." << std::endl;
   fInit       = true;  
}
//______________________________________________________________________________
int BDEventGen::GenerateEvent(){
   // init values if neeeded 
   if(!fInit) Initialize(); 
   // compute the vertex based on the beam raster and pointing 
   fX = fPointingX + CLHEP::RandFlat::shoot(-fRasterX/2.,fRasterX/2.);
   fY = fPointingY + CLHEP::RandFlat::shoot(-fRasterY/2.,fRasterY/2.);
   // set in vertex vector
   fVertex.setX(fX); 
   fVertex.setY(fY); 
   fVertex.setZ(0);   // z location is locked to zero (set in PrimaryGenerator class)  
   // set beam energy 
   fBeamE = CLHEP::RandGauss::shoot(fBeamE_in,fBeamE_sig);
 
   // beam angle; we project back from a random spot on the dump plane
   // G4double bd_size = 100*CLHEP::mm; // defines a radius of 5 mm  
   // G4double bd_x = CLHEP::RandFlat::shoot(-bd_size/2.,bd_size/2.);
   // G4double bd_y = CLHEP::RandFlat::shoot(-bd_size/2.,bd_size/2.);

   // G4double bd_L = 32.*CLHEP::m + fZ0;  // account for distance to z origin.  BD is at 32 m downstream of target 
   // fBeamAngleX = atan(bd_x/bd_L); 
   // fBeamAngleY = atan(bd_y/bd_L);
   // fBeamAngleZ = 0;               // do we really need this?   

   // randomize by some percent around chosen value 
   G4double pct  = fBeamAngleSpread;  
   G4double bd_L = 32.*CLHEP::m + fZ0;  // account for distance to z origin.  BD is at 32 m downstream of target 
   G4double bd_x = bd_L*tan(fBeamAngleX); 
   G4double bd_y = bd_L*tan(fBeamAngleY);
   bd_x = CLHEP::RandGauss::shoot(bd_x,bd_x*pct); 
   bd_y = CLHEP::RandGauss::shoot(bd_y,bd_y*pct);

   // compute new angles 
   fBeamAngleX = atan(bd_x/bd_L); 
   fBeamAngleY = atan(bd_y/bd_L);
   fBeamAngleZ = 0;               // do we really need this?   

   // beam angle (smear by some amount, only if angle is nonzero) 
   // if(fBeamAngleX!=0) fBeamAngleX = CLHEP::RandGauss::shoot(fBeamAngleX,fbax_sig);  
   // if(fBeamAngleY!=0) fBeamAngleY = CLHEP::RandGauss::shoot(fBeamAngleY,fbay_sig);  
   // if(fBeamAngleZ!=0) fBeamAngleZ = CLHEP::RandGauss::shoot(fBeamAngleZ,fbaz_sig);  

   if(fBeamE<0.5*CLHEP::MeV){
      return 1;
   }else{
      Print();
      return 0; 
   }
}
//______________________________________________________________________________
void BDEventGen::Print(){
   char msg[200];
   std::cout << "[BDEventGen]: Generating an event!" << std::endl;
   sprintf(msg,"[BDEventGen]: raster: x = %.1lf mm, y = %.1lf mm",fRasterX,fRasterY); 
   std::cout << msg << std::endl;
   sprintf(msg,"[BDEventGen]: pointing: x = %.1lf mm, y = %.1lf mm",fPointingX,fPointingY); 
   std::cout << msg << std::endl;
   sprintf(msg,"[BDEventGen]: beam angle: x = %.3lf mrad, y = %.3lf mrad, z = %.3lf mrad",
           fBeamAngleX/CLHEP::mrad,fBeamAngleY/CLHEP::mrad,fBeamAngleZ/CLHEP::mrad); 
   std::cout << msg << std::endl;

   sprintf(msg,"[BDEventGen]: random x = %.3lf mm, y = %.3lf mm",fX/CLHEP::mm,fY/CLHEP::mm); 
   std::cout << msg << std::endl;
   sprintf(msg,"[BDEventGen]: E_in  = %.3lf GeV, E_sig = %.3lf MeV, random beam energy = %.3lf GeV",
	 fBeamE_in/CLHEP::GeV,fBeamE_sig/CLHEP::MeV,fBeamE/CLHEP::GeV); 
   std::cout << msg << std::endl;
   std::cout << "-----------------------------------------------" << std::endl;
}
