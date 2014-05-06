#ifndef KineUtils_hh
#define KineUtils_hh

#include <TObject.h>
#include <TVector3.h>

class KineUtils 
{
public:

  virtual ~KineUtils(){}
  
  //  calculate phi from x, y (between 0 and 2*pi)
  static float phi( float x, float y );

  //  adjust phi1 to the value within +/- phi of phi2
  static void adjust( float& phi1, float phi2 );

  //  calculate phi1-phi2 keeping value between 0 and pi
  static float dPhi( float ph11, float phi2 );

  //  calculate phi1-phi2 keeping value between 0 and pi
  static float dR( float eta1, float eta2, float ph11, float phi2 );

  //  calculate theta from eta
  static float theta( float eta_ );

  //  calculate eta from theta
  static float eta( float theta_ );

  //  calculate eta from x, y, z
  static float eta( float xx_, float yy_, float zz_ );

  // calculate et from e and eta
  static float et( float e_, float eta_ );

  // calculate e from et and eta
  static float e( float et_, float eta_ );

  // calculate rapidity from E, pz
  static float y( float E, float pz );

  static float distanceToEnvelop( float x0, float y0, float z0,  
				 float vx, float vy, float vz,
				 float rmax, float zmax );

  static float straightToEnvelop( float& eta_, float& phi_, 
				  float eta0, float phi0,
				  float x0, float y0, float z0, 
				  float rmax, float zmax );

  static float straightToEnvelop( float& eta_, float& phi_, 
				  float x0, float y0, float z0, 
				  float vx, float vy, float vz, 
				  float rmax, float zmax );
  
  static void helixToEnvelop( float& eta, float& phi, 
			      float pt0, float eta0, float phi0, 
			      float x0, float y0, float z0, 
			      float rmax, float zmax, 
			      float B );			      

  //d0
  static float d0(const TVector3& pv,const TVector3& vtx, TVector3& p);


  ClassDef( KineUtils, 0 )

};
#endif
