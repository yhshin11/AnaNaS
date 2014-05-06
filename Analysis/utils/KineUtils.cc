#include "Analysis/utils/KineUtils.hh"

#include <cmath>
#include <iostream>
#include <cassert>
using namespace std;

#include "Analysis/utils/Constants.hh"

ClassImp( KineUtils )

float 
KineUtils::phi( float x, float y )
{
  float phi_ =atan2(y, x);
  return (phi_>=0) ?  phi_ : phi_ + 2*Constants::pi;
}

float 
KineUtils::dPhi( float phi1, float phi2 )
{
  float phi1_= phi( cos(phi1), sin(phi1) );
  float phi2_= phi( cos(phi2), sin(phi2) );
  float dphi_= phi1_-phi2_;
  if( dphi_> Constants::pi ) dphi_-=2*Constants::pi;
  if( dphi_<-Constants::pi ) dphi_+=2*Constants::pi;
  return dphi_;
}

void 
KineUtils::adjust( float& phi1, float phi2 )
{
  phi1 = phi2 + dPhi( phi1, phi2 );
}

float 
KineUtils::theta( float eta_ )
{
  float th_ =  2.0*atan(exp(-eta_));
  if( th_<0 ) th_ += 2*Constants::pi;
  return th_;

}

float 
KineUtils::eta( float theta_ )
{
  if( sin(theta_/2.)==0 ) return 10000.*cos(theta_/2.);
  return -log(tan(theta_/2.0));
}

float 
KineUtils::eta( float xx_, float yy_, float zz_ )
{
  float rr_ = sqrt( xx_*xx_+yy_*yy_ );
  if( rr_==0 ) return (zz_<0)?-10000:10000;
  float theta_=Constants::pi/2.;
  if( fabs(zz_)>0 ) theta_ -= atan( zz_/rr_ );
  return eta( theta_);
}

float 
KineUtils::y( float E, float pz )
{
  return 0.5 * log ((E+pz+1.e-10)/(E-pz+1.e-10));
}

float
KineUtils::dR( float eta1, float eta2, float phi1, float phi2 )
{
  float dEta_ = eta2-eta1;
  float dPhi_ = dPhi( phi1, phi2 );
  float dR_ = sqrt( dEta_*dEta_ + dPhi_*dPhi_ );

  return dR_;
}

float
KineUtils::et( float e_, float eta_ )
{
  float theta_ = theta( eta_ );
  return e_*sin(theta_);
}

float
KineUtils::e( float et_, float eta_ )
{
  float theta_ = theta( eta_ );
  float stheta_ = sin(theta_);
  assert( stheta_>0 );
  return et_/stheta_;
}

float
KineUtils::distanceToEnvelop( float x0, float y0, float z0, 
			      float vx, float vy, float vz,
			      float rmax, float zmax )
{
  float rmax2 = rmax*rmax;
  float r2 = x0*x0 + y0*y0;

  // distance to the envelop
  float dz_=10000;
  if( vz>0 )
    {
      dz_=(zmax-z0)/vz;
    }
  else if( vz<0 )
    {
      dz_=(-zmax-z0)/vz;
    }
  if( dz_<0 ) return -1;

  float aa_ = vx*vx + vy*vy;
  assert( aa_>0 );
  float bb_ = 2*(x0*vx+y0*vy);
  float cc_ = r2-rmax2;
  if( aa_*cc_>0 ) return -1;
  float D_ = bb_*bb_-4*aa_*cc_;
  assert( D_>0 );
  float dr_ = (-bb_+sqrt(D_))/(2*aa_);
  if( dr_<0 ) return -1;

  float dd_ = dr_;
  if( dr_>dz_ ) dd_ = dz_;

  return dd_;
}

float
KineUtils::straightToEnvelop( float& eta_, float& phi_, 
		      float eta0, float phi0,
		      float x0, float y0, float z0, 
		      float rmax, float zmax )
{
  float theta_ = theta(eta0);
  float ctheta  = cos(theta_);
  float stheta  = sin(theta_);
  
  float cphi  = cos( phi0 );
  float sphi  = sin( phi0 );
  
  float vx_ = cphi*stheta;
  float vy_ = sphi*stheta;
  float vz_ = ctheta;
  return
    straightToEnvelop( eta_, phi_, x0, y0, z0, vx_, vy_, vz_, rmax, zmax );
}

float
KineUtils::straightToEnvelop( float& eta_, float& phi_, 
			      float x0, float y0, float z0, 
			      float vx, float vy, float vz,
			      float rmax, float zmax )
{
  float rmax2 = rmax*rmax;

  float r2 = x0*x0 + y0*y0;
  if( r2>=rmax2 ) return 0;
  if( fabs(z0)>=zmax ) return 0;
  
  float dd_ = distanceToEnvelop(  x0, y0, z0, 
				  vx, vy, vz,
				  rmax, zmax   );
  if( dd_<0 ) return 0;

  // intersection tangent/envelop
  float xx_ = x0 + dd_*vx;
  float yy_ = y0 + dd_*vy;
  float zz_ = z0 + dd_*vz;
  
  phi_ = phi( xx_, yy_ );
  eta_ = eta( xx_, yy_, zz_ );

  return dd_;
}

void
KineUtils::helixToEnvelop( float& etaAtSurf, float& phiAtSurf,
			   float pt0, float eta0, float phi0, 
			   float x0, float y0, float z0, 
			   float rmax, float zmax, 
			   float B )
{
  //
  // Warning: pt0 is signed
  //
  float rmax2 = rmax*rmax;

  // helix...
  float theta = 2*atan( exp( -eta0 ) );
  if( theta<0 ) theta += 2*Constants::pi;
  float lambda = Constants::pi/2 -theta;
  float tanl   = tan(lambda);
  float ctheta  = cos(theta);
  float stheta  = sin(theta);
  float ttheta  = tan( theta );
  assert( ttheta!=0 );
  assert( stheta>0 );
  float phi_   = phi0;
  float cphi_  = cos( phi_ );
  float sphi_  = sin( phi_ );
  float pt    = fabs(pt0);
  assert( pt!=0 );

  float chg   = pt0/pt;
  float vx =  cphi_;
  float vy =  sphi_;
  float ux = -chg*sphi_;
  float uy =  chg*cphi_;

  // leading cosines at vertex
  float vx0 =  vx*stheta;
  float vy0 =  vy*stheta;
  float vz0 =  ctheta;

  // distance to envelop
  float dd_ = straightToEnvelop( etaAtSurf, phiAtSurf, 
				 x0, y0, z0, vx0, vy0, vz0, rmax, zmax );

  // radius of curvature
  float R  =  pt/(0.003*B);
  
  // center of the trajectory
  float cx =  x0 - R*ux;
  float cy =  y0 - R*uy;
  float cz =  z0;
      
  // angle at initial point along the trajectory
  float a = dd_*stheta/R;
  int npt = 0;

  //  cout << "eta0/phi0/tanl=" << eta0 << "/" << phi0 << "/" << tanl << endl;
  //  cout << "x0/y0/z0=" << x0 << "/" << y0 << "/" << z0 << endl;

  float x_ = x0;
  float y_ = y0;
  float z_ = y0;
  float vxa_ = vx0;
  float vya_ = vy0;
  float vza_ = vz0;

  //  unsigned maxpt=5000;
  while( dd_>0.01 )
    {
      //      if( npt>4 ) break;

      float cosa = cos(a);
      float sina = sin(a);
      float uxa_ = cosa*ux + sina*vx;
      float uya_ = cosa*uy + sina*vy;

      // the position
      x_ = cx + R*uxa_;
      y_ = cy + R*uya_;
      z_ = cz + a*R*tanl;

      // the direction
      vxa_ =  chg*uya_*stheta;
      vya_ = -chg*uxa_*stheta;
      vza_ =  ctheta;

      float r2 = x_*x_ + y_*y_;
      //      cout << "\nnpt=" << npt << endl;
      //      cout << "x/y=" << x_ << "/" << y_ << endl;
      //      cout << "r/rmax=" << sqrt(r2) << "/" << rmax << endl;
      //      cout << "z/zmax=" << z_ << "/" << zmax << endl;
      //      cout << "vx/vy/vz=" << vxa_ << "/" << vya_ << "/" << vza_ << endl;
      //      cout << "eta/phi=" << etaAtSurf << "/" << phiAtSurf << endl;
      
      //      cout << "--> x_/y_/z_ " << x_ << "/" << y_ << "/" << z_ << endl;

      //      cout << "straight --> eta/phi " << etaAtSurf << "/" << phiAtSurf << endl;
      etaAtSurf = eta( x_, y_, z_ );
      phiAtSurf = phi( x_, y_ );
      //      cout << "helix    --> eta/phi " << etaAtSurf << "/" << phiAtSurf << endl;
 
      if( r2>rmax2 )      break;
      if( fabs(z_)>zmax ) break;

      if( a>3*Constants::pi/2. ) break;  // against loopers

      
      dd_ = straightToEnvelop(  etaAtSurf, phiAtSurf, 
				x_, y_, z_, vxa_, vya_, vza_, rmax, zmax );
      if( dd_<0 ) break;

      a += dd_*stheta/R;

      //      cout << "a=" << a << endl;
      //      cout << "dd=" << dd_ << endl;

      npt++;
    }
}  



float
KineUtils::d0(const TVector3& pv, const TVector3& vtx, TVector3& p) {

  return ( -(vtx.X()-pv.X() )*p.Py() + (vtx.Py()-pv.Y() )*p.Px() ) / p.Pt();

}
