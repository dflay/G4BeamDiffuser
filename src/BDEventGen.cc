#include "BDEventGen.hh"
//______________________________________________________________________________
BDEventGen::BDEventGen(){
   fInit      = false;
   fBeamE     = 0.;  // generated beam energy
   fBeamE_in  = 0.;  // input beam energy
   fBeamE_sig = 0.;
   fBeamCur   = 0.;
   fRasterX   = 0.;
   fRasterY   = 0.;
   fPointingX = 0.;
   fPointingY = 0.;
   fX         = 0.;
   fY         = 0.;
   fVertex    = G4ThreeVector(0,0,0);
   fBeamP     = G4ThreeVector(0,0,0);  
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

   if(fBeamE<0.5*CLHEP::MeV){
      return 1;
   }else{
      // Print();
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
   sprintf(msg,"[BDEventGen]: random x = %.3lf mm, y = %.3lf mm",fX/CLHEP::mm,fY/CLHEP::mm); 
   std::cout << msg << std::endl;
   sprintf(msg,"[BDEventGen]: E_in  = %.3lf GeV, E_sig = %.3lf MeV, random beam energy = %.3lf GeV",
	 fBeamE_in/CLHEP::GeV,fBeamE_sig/CLHEP::MeV,fBeamE/CLHEP::GeV); 
   std::cout << msg << std::endl;
   std::cout << "-----------------------------------------------" << std::endl;
}
