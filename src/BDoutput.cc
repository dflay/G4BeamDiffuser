#include "BDoutput.hh"
//______________________________________________________________________________
BDoutput::BDoutput(){
   Clear();
}
//______________________________________________________________________________
BDoutput::~BDoutput(){

}
//______________________________________________________________________________
void BDoutput::Clear(){
   nhits = 0;
   plane.clear();
   trid.clear();
   t.clear();
   x.clear();
   y.clear();
   z.clear();
   xg.clear();
   yg.clear();
   zg.clear();
   p.clear();
   pid.clear();
   mid.clear();
   beta.clear();
   edep.clear();
}
