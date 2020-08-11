#include "BDEventGen.hh"
//______________________________________________________________________________
BDEventGen::BDEventGen(){
   fBeamE     = 0.;
   fBeamCur   = 0.;
   fRasterX   = 0.;
   fRasterY   = 0.;
   fPointingX = 0.;
   fPointingY = 0.;
   fVertex    = G4ThreeVector(0,0,0);
   fBeamP     = G4ThreeVector(0,0,0);  
}
//______________________________________________________________________________
BDEventGen::~BDEventGen(){

}
//______________________________________________________________________________
int BDEventGen::GenerateEvent(){
   std::cout << "[BDEventGen]: Generating an event!" << std::endl;
   char msg[200];
   sprintf(msg,"[BDEventGen]: raster: x = %.1lf mm, y = %.1lf mm",fRasterX,fRasterY); 
   std::cout << msg << std::endl;
   sprintf(msg,"[BDEventGen]: pointing: x = %.1lf mm, y = %.1lf mm",fPointingX,fPointingY); 
   std::cout << msg << std::endl;
   // compute the vertex based on the beam raster and pointing 
   G4double vx = fPointingX + CLHEP::RandFlat::shoot(-fRasterX/2.,fRasterX/2.);
   G4double vy = fPointingY + CLHEP::RandFlat::shoot(-fRasterY/2.,fRasterY/2.);
   // set in vertex vector
   fVertex.setX(vx); 
   fVertex.setY(vy); 
   fVertex.setZ(0);   // z location is locked to zero (set in PrimaryGenerator class)  
   return 0; 
}
