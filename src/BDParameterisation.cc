#include "BDParameterisation.hh"
//______________________________________________________________________________
BDParameterisation::BDParameterisation(char Hall,G4ThreeVector r0){
   fHall    = Hall; // hall label 
   fR0      = r0;   // origin of part (relative to mother) 
   InitParameters();
}
//______________________________________________________________________________
BDParameterisation::~BDParameterisation(){
   delete[] fThickness; 
   delete[] fStartPhi; 
   delete[] fDeltaPhi; 
   delete[] fColor; 
}
//______________________________________________________________________________
void BDParameterisation::InitParameters(){

   if(fHall!='A' && fHall!='C'){
      std::cout << "[BDParameterisation::InitParameters]: Invalid Hall = " 
                << fHall << std::endl;
      exit(1);
   }

   // initialize materials  
   auto nistManager = G4NistManager::Instance();
   nistManager->FindOrBuildMaterial("G4_Al");

   // define some constants 
   G4double inch        = 25.4*mm; 
   G4double thk_031     = 0.03125*inch;
   G4double thk_062     = 0.06250*inch;
   G4double thk_094     = 0.09375*inch;
   G4double thk_100     = 0.10000*inch; 
   G4double thk_125     = 0.12500*inch; 
   G4double startPhi_30 = 255.*deg; 
   G4double startPhi_60 = 255.*deg;
   G4double startPhi_90 = 225.*deg;
   G4double deltaPhi_30 = 30.*deg; 
   G4double deltaPhi_60 = 60.*deg; 
   G4double deltaPhi_90 = 90.*deg; 

   fGap        = 0.195*inch; 
   fWidth      = 2.*inch; 
   fRadius_min = 5.*inch; 
   fRadius_max = fRadius_min + fWidth; 

   // make arrays large enough for all the layers we may need 
   const int NP_MAX = 16; 
   fThickness = new double[NP_MAX];
   fStartPhi  = new double[NP_MAX];
   fDeltaPhi  = new double[NP_MAX];
   fColor     = new int[NP_MAX];
   for(int i=0;i<NP_MAX;i++){
      fThickness[i] = 0;
      fStartPhi[i]  = 0;
      fDeltaPhi[i]  = 0;
      fColor[i]     = 0;
   } 

   // set thickness array 
   if(fHall=='A'){
      fNLayers = 15;
      // thicknesses  
      fThickness[0]  = thk_125; fThickness[1]  = thk_125; fThickness[2]  = thk_100; fThickness[3]  = thk_125;
      fThickness[4]  = thk_125; fThickness[5]  = thk_100; fThickness[6]  = thk_125; fThickness[7]  = thk_100; 
      fThickness[8]  = thk_125; fThickness[9]  = thk_100; fThickness[10] = thk_125; fThickness[11] = thk_125; 
      fThickness[12] = thk_100; fThickness[13] = thk_125; fThickness[14] = thk_125;
      // start angles 
      fStartPhi[0]   = startPhi_30; fStartPhi[1]  = startPhi_30; fStartPhi[2]  = startPhi_30; fStartPhi[3]  = startPhi_60; 
      fStartPhi[4]   = startPhi_60; fStartPhi[5]  = startPhi_60; fStartPhi[6]  = startPhi_90; fStartPhi[7]  = startPhi_90;
      fStartPhi[8]   = startPhi_90; fStartPhi[9]  = startPhi_60; fStartPhi[10] = startPhi_60; fStartPhi[11] = startPhi_60;
      fStartPhi[12]  = startPhi_30; fStartPhi[13] = startPhi_30; fStartPhi[14] = startPhi_30;   
      // step angles 
      fDeltaPhi[0]   = deltaPhi_30; fDeltaPhi[1]  = deltaPhi_30; fDeltaPhi[2]  = deltaPhi_30; fDeltaPhi[3]  = deltaPhi_60; 
      fDeltaPhi[4]   = deltaPhi_60; fDeltaPhi[5]  = deltaPhi_60; fDeltaPhi[6]  = deltaPhi_90; fDeltaPhi[7]  = deltaPhi_90;
      fDeltaPhi[8]   = deltaPhi_90; fDeltaPhi[9]  = deltaPhi_60; fDeltaPhi[10] = deltaPhi_60; fDeltaPhi[11] = deltaPhi_60;
      fDeltaPhi[12]  = deltaPhi_30; fDeltaPhi[13] = deltaPhi_30; fDeltaPhi[14] = deltaPhi_30;   
      // colors 
      fColor[0]      = diffuser::BLUE;    fColor[1]  = diffuser::BLUE;    fColor[2]  = diffuser::MAGENTA; fColor[3]  = diffuser::BLUE;
      fColor[4]      = diffuser::BLUE;    fColor[5]  = diffuser::MAGENTA; fColor[6]  = diffuser::BLUE;    fColor[7]  = diffuser::MAGENTA;
      fColor[8]      = diffuser::BLUE;    fColor[9]  = diffuser::MAGENTA; fColor[10] = diffuser::BLUE;    fColor[11] = diffuser::BLUE;
      fColor[12]     = diffuser::MAGENTA; fColor[13] = diffuser::BLUE;    fColor[14] = diffuser::BLUE;
   }else if(fHall=='C'){
      fNLayers = 16;  
      // thicknesses  
      fThickness[0]  = thk_125; fThickness[1]  = thk_125; fThickness[2]  = thk_125; fThickness[3]  = thk_125;
      fThickness[4]  = thk_125; fThickness[5]  = thk_125; fThickness[6]  = thk_094; fThickness[7]  = thk_094; 
      fThickness[8]  = thk_094; fThickness[9]  = thk_094; fThickness[10] = thk_094; fThickness[11] = thk_094; 
      fThickness[12] = thk_062; fThickness[13] = thk_062; fThickness[14] = thk_031; fThickness[15] = thk_031; 
      // start angles 
      fStartPhi[0]   = startPhi_30; fStartPhi[1]  = startPhi_30; fStartPhi[2]  = startPhi_30; fStartPhi[3]  = startPhi_30; 
      fStartPhi[4]   = startPhi_60; fStartPhi[5]  = startPhi_60; fStartPhi[6]  = startPhi_90; fStartPhi[7]  = startPhi_90;
      fStartPhi[8]   = startPhi_90; fStartPhi[9]  = startPhi_60; fStartPhi[10] = startPhi_60; fStartPhi[11] = startPhi_60;
      fStartPhi[12]  = startPhi_30; fStartPhi[13] = startPhi_30; fStartPhi[14] = startPhi_30; fStartPhi[15] = startPhi_30;
      // step angles 
      fDeltaPhi[0]   = deltaPhi_30; fDeltaPhi[1]  = deltaPhi_30; fDeltaPhi[2]  = deltaPhi_30; fDeltaPhi[3]  = deltaPhi_30; 
      fDeltaPhi[4]   = deltaPhi_60; fDeltaPhi[5]  = deltaPhi_60; fDeltaPhi[6]  = deltaPhi_90; fDeltaPhi[7]  = deltaPhi_90;
      fDeltaPhi[8]   = deltaPhi_90; fDeltaPhi[9]  = deltaPhi_60; fDeltaPhi[10] = deltaPhi_60; fDeltaPhi[11] = deltaPhi_60;
      fDeltaPhi[12]  = deltaPhi_30; fDeltaPhi[13] = deltaPhi_30; fDeltaPhi[14] = deltaPhi_30; fDeltaPhi[15] = deltaPhi_30;
      // colors 
      fColor[0]      = diffuser::BLUE;   fColor[1]  = diffuser::BLUE;   fColor[2]  = diffuser::BLUE;  fColor[3]  = diffuser::BLUE;
      fColor[4]      = diffuser::BLUE;   fColor[5]  = diffuser::BLUE;   fColor[6]  = diffuser::GREEN; fColor[7]  = diffuser::GREEN;
      fColor[8]      = diffuser::GREEN;  fColor[9]  = diffuser::GREEN;  fColor[10] = diffuser::GREEN; fColor[11] = diffuser::GREEN;
      fColor[12]     = diffuser::YELLOW; fColor[13] = diffuser::YELLOW; fColor[14] = diffuser::RED;   fColor[15] = diffuser::RED;
   }
}
//______________________________________________________________________________
void BDParameterisation::ComputeTransformation(const G4int copyNo, 
      G4VPhysicalVolume *physVol) const{
   // z coordinate: distance along z to *center* of the plate  
   // sum over previous layers 
   G4double Ls=0;
   for(int i=0;i<copyNo;i++) Ls += fThickness[i]; 
   // put it all together 
   G4double x  = fR0.x(); 
   G4double y  = fR0.y(); 
   G4double z0 = fR0.z();
   G4double z  = z0 + Ls + (double)(copyNo-1)*fGap + 0.5*fThickness[copyNo]; 
   // set the 3-vector
   G4ThreeVector P = G4ThreeVector(x,y,z); 
   physVol->SetTranslation(P); // set position 
   physVol->SetRotation(0);    // no rotation  
}
//______________________________________________________________________________
void BDParameterisation::ComputeDimensions(G4Tubs &plate,
      const G4int copyNo,const G4VPhysicalVolume *physVol) const{
   plate.SetInnerRadius(fRadius_min);
   plate.SetOuterRadius(fRadius_max);
   plate.SetZHalfLength(fThickness[copyNo]/2.);
   plate.SetStartPhiAngle(fStartPhi[copyNo]); 
   plate.SetDeltaPhiAngle(fDeltaPhi[copyNo]);

   // determine color by copy number 
   G4VisAttributes *vis = new G4VisAttributes(); 
   if(fColor[copyNo]==diffuser::RED     ) vis->SetColour( G4Colour::Red()     );  
   if(fColor[copyNo]==diffuser::YELLOW  ) vis->SetColour( G4Colour::Yellow()  );  
   if(fColor[copyNo]==diffuser::GREEN   ) vis->SetColour( G4Colour::Green()   );  
   if(fColor[copyNo]==diffuser::BLUE    ) vis->SetColour( G4Colour::Blue()    ); 
   if(fColor[copyNo]==diffuser::MAGENTA ) vis->SetColour( G4Colour::Magenta() ); 

   // set properties 
   physVol->GetLogicalVolume()->SetVisAttributes(vis); 
   // physVol->GetLogicalVolume()->SetMaterial(theMaterial);  
 
}
// //______________________________________________________________________________
// G4VSolid *BDParameterisation::ComputeSolid(const G4int copyNo, 
//       G4VPhysicalVolume *physVol){
//    char name[200];
//    sprintf(name,"plate_%02d",copyNo); 
//    G4VSolid *solid = new G4Tubs(name,
//                                 fRadius_min,fRadius_max,fThickness[copyNo]/2.,
//                                 fStartPhi[copyNo],fDeltaPhi[copyNo]); 
//    return solid; 
// }
// //______________________________________________________________________________
// G4Material *BDParameterisation::ComputeMaterial(const G4int copyNo, 
//       G4VPhysicalVolume *physVol,const G4VTouchable *parentTouch){
//    // each plate is made of aluminum  
//    G4Material *theMaterial = G4Material::GetMaterial("G4_Al");  
// 
//    // determine color by copy number 
//    G4VisAttributes *vis = new G4VisAttributes(); 
//    if(fColor[copyNo]==diffuser::RED     ) vis->SetColour( G4Colour::Red()     );  
//    if(fColor[copyNo]==diffuser::YELLOW  ) vis->SetColour( G4Colour::Yellow()  );  
//    if(fColor[copyNo]==diffuser::GREEN   ) vis->SetColour( G4Colour::Green()   );  
//    if(fColor[copyNo]==diffuser::BLUE    ) vis->SetColour( G4Colour::Blue()    ); 
//    if(fColor[copyNo]==diffuser::MAGENTA ) vis->SetColour( G4Colour::Magenta() ); 
// 
//    // set properties 
//    physVol->GetLogicalVolume()->SetVisAttributes(vis); 
//    physVol->GetLogicalVolume()->SetMaterial(theMaterial);  
// 
//    return theMaterial;
// }
