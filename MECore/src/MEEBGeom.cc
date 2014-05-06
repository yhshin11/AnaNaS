#include <iostream>
#include <stdlib.h>
#include <string>
#include <assert.h>
using namespace std;

#include "ME.hh"
#include "MEEBGeom.hh"

#include <TGraph.h>

ClassImp(MEEBGeom)

int 
MEEBGeom::barrel( EBGlobalCoord ieta, EBGlobalCoord iphi )
{
  int iz=1;
  if( ieta<0 ) iz=-1;
  ieta *= iz;
  assert( ieta>0 && ieta<=85 );
  if( iphi<0 ) iphi+=360;
  assert( iphi>0 && iphi<=360 );
  return iz;
}

int 
MEEBGeom::sm( EBGlobalCoord ieta, EBGlobalCoord iphi )
{
  int iz=1;
  if( ieta<0 ) iz=-1;
  ieta *= iz;
  assert( ieta>0 && ieta<=85 );

  //  int iphi_ = iphi+10;
  int iphi_ = iphi;
  if( iphi_>360 ) iphi_-=360;
  assert( iphi_>0 && iphi_<=360 );

  int ism = (iphi_-1)/20 + 1;
  assert( ism>=1 && ism<=18 );
  //  if( iz==1 ) ism += 18;
  if( iz==-1 ) ism += 18;

  return ism;
}

int 
MEEBGeom::dcc( EBGlobalCoord ieta, EBGlobalCoord iphi )
{
  int ism = sm( ieta, iphi );
  return dccFromSm( ism );
}

int
MEEBGeom::dccFromSm( int ism )
{
  assert( ism>=1 && ism<=36 );
  int iz=1;
  if( ism>18 ) iz=-1;
  if( iz==-1  ) ism-=18;
  assert( ism>=1 && ism<=18 );
  int idcc = 9+ism;
  if( iz==+1  ) idcc+=18;
  return idcc;
}

int
MEEBGeom::smFromDcc( int idcc )
{
  if( idcc>600 ) idcc-=600;  // also works with FEDids
  assert( idcc>=10 && idcc<=45 );
  int ism=idcc-9;
  if( ism>18 ) ism-=18;
  else         ism+=18;
  return ism;
}

TString
MEEBGeom::smName( int ism )
{
  assert( ism>=1 && ism<=36 );
  TString out = "EB+";
  if( ism>18 )
    {
      out = "EB-";
      ism -= 18;
    }
  out += ism;
  return out;
}

int 
MEEBGeom::lmmod( EBGlobalCoord ieta, EBGlobalCoord iphi )
{
  pair<EBLocalCoord,EBLocalCoord> ixy = localCoord( ieta, iphi );
  return lm_channel( ixy.first/5, ixy.second/5 );
}

int 
MEEBGeom::tt( EBGlobalCoord ieta, EBGlobalCoord iphi )
{
  pair<EBLocalCoord,EBLocalCoord> ixy = localCoord( ieta, iphi );
  return tt_channel( ixy.first/5, ixy.second/5 );
}

int 
MEEBGeom::crystal( EBGlobalCoord ieta, EBGlobalCoord iphi )
{
  pair<EBLocalCoord,EBLocalCoord> ixy = localCoord( ieta, iphi );
  return crystal_channel( ixy.first, ixy.second );
}

int 
MEEBGeom::side( EBGlobalCoord ieta, EBGlobalCoord iphi )
{
  int ilmmod = lmmod( ieta, iphi );
  return (ilmmod%2==0)?1:0;
}

int 
MEEBGeom::lmr( EBGlobalCoord ieta, EBGlobalCoord iphi )
{
  int idcc  = dcc( ieta, iphi );
  int ism   = idcc-9;
  int iside = side( ieta, iphi );
  int ilmr = 1 + 2*(ism-1) + iside;
  return ilmr;
}

pair< MEEBGeom::EBLocalCoord, MEEBGeom::EBLocalCoord >    
MEEBGeom::localCoord( MEEBGeom::EBGlobalCoord ieta, MEEBGeom::EBGlobalCoord iphi )
{
  int iz=1;
  if( ieta<0 ) iz=-1;
  ieta *= iz;
  assert( ieta>0 && ieta<=85 );

  //  int iphi_ = iphi+10;
  int iphi_ = iphi;
  if( iphi_>360 ) iphi_-=360;
  assert( iphi_>0 && iphi_<=360 );

  int ix = ieta-1;
  
  int iy = (iphi_-1)%20;
  if( iz==-1 ) iy = 19-iy;
  //***** if( iz==1 ) iy = 19-iy;

  return pair< EBLocalCoord, EBLocalCoord >(ix,iy);  
}

pair< MEEBGeom::EBLocalCoord, MEEBGeom::EBLocalCoord >    
MEEBGeom::localCoord( int icr )
{
  assert( icr>=1 && icr<=1700 );
  int ix = (icr-1)/20;
  int iy = 19 - (icr-1)%20;
  return pair< EBLocalCoord, EBLocalCoord >(ix,iy);  
}

pair< MEEBGeom::EBGlobalCoord, MEEBGeom::EBGlobalCoord >    
MEEBGeom::globalCoord( int ism, MEEBGeom::EBLocalCoord ix, MEEBGeom::EBLocalCoord iy )
{
  assert( ism>=1 && ism<=36 );
  assert( ix>=0 && ix<85 );
  assert( iy>=0 && iy<20 );
  //  int iz=-1;
  int iz=1;
  if( ism>18 ) 
    {
      iz=-1;
      ism -= 18;
    }
  if( iz==-1 ) iy = 19-iy;
  //**** if( iz==1 ) iy = 19-iy;

  int ieta = ix+1;
  ieta *= iz;
  //  int iphi = -9 + iy + 20*(ism-1);
  int iphi = 1 + iy + 20*(ism-1);

  return pair< EBGlobalCoord, EBGlobalCoord >(ieta,iphi);  
}

pair< float, float >    
MEEBGeom::globalCoord( int ism, float x, float y )
{
  assert( ism>=1 && ism<=36 );
  int iz=1;
  if( ism>18 ) 
    {
      iz=-1;
      ism -= 18;
    }
  if( iz==-1 ) y = 19-y;

  float eta = x+1;
  eta *= iz;
  //  float phi = -9 + y + 20*(ism-1);
  float phi = 1 + y + 20*(ism-1);

  return pair< float, float >(eta,phi);  
}

pair< MEEBGeom::EBGlobalCoord, MEEBGeom::EBGlobalCoord >    
MEEBGeom::globalCoord( int ism, int icr )
{
  assert( ism>=1 && ism<=36 );
  assert( icr>=1 && icr<=1700 );

  int ix = (icr-1)/20;
  int iy = 19 - (icr-1)%20;

  return globalCoord( ism, ix, iy );  
}

int 
MEEBGeom::lm_channel( EBTTLocalCoord iX, EBTTLocalCoord iY )
{
  string 
    str_
    (
     // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
     "-01-02-02-02-02-04-04-04-04-06-06-06-06-08-08-08-08" // 3
     "-01-02-02-02-02-04-04-04-04-06-06-06-06-08-08-08-08" // 2
     "-01-03-03-03-03-05-05-05-05-07-07-07-07-09-09-09-09" // 1
     "-01-03-03-03-03-05-05-05-05-07-07-07-07-09-09-09-09" // 0
     // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
     );

  int iym, ixm, nc, il, ic, ii;
  iym=4;
  ixm=17;
  int iX_ = iX+1;
  int iY_ = iY+1;
  nc=3;		  
  il=iym-iY_;
  ic=nc*(iX_-1);
  ii=il*nc*ixm+ic;
  if( str_.substr(ii,1).find("-")==string::npos ) return -1;
  return atoi( str_.substr(ii+1,2).c_str() );
}

int 
MEEBGeom::tt_type( EBTTLocalCoord iX, EBTTLocalCoord iY )
{
  string 
    str_
    (
     // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
     "-01-01-01-02-02-01-01-02-02-01-01-02-02-01-01-02-02" // 3
     "-01-01-01-02-02-01-01-02-02-01-01-02-02-01-01-02-02" // 2
     "-01-01-01-02-02-01-01-02-02-01-01-02-02-01-01-02-02" // 1
     "-01-01-01-02-02-01-01-02-02-01-01-02-02-01-01-02-02" // 0
     // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
     );

  int iym, ixm, nc, il, ic, ii;
  iym=4;
  ixm=17;
  int iX_ = iX+1;
  int iY_ = iY+1;
  nc=3;		  
  il=iym-iY_;
  ic=nc*(iX_-1);
  ii=il*nc*ixm+ic;
  if( str_.substr(ii,1).find("-")==string::npos ) return -1;
  return atoi( str_.substr(ii+1,2).c_str() );
}

int
MEEBGeom::hv_channel( EBTTLocalCoord iX, EBTTLocalCoord iY )
{
  string 
    str_
    (
     // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
     "-01-03-05-07-09-11-13-15-17-19-21-23-25-27-29-31-33" // 3
     "-01-03-05-07-09-11-13-15-17-19-21-23-25-27-29-31-33" // 2
     "-02-04-06-08-10-12-14-16-18-20-22-24-26-28-30-32-34" // 1
     "-02-04-06-08-10-12-14-16-18-20-22-24-26-28-30-32-34" // 0
     // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
     );

  int iym, ixm, nc, il, ic, ii;
  iym=4;
  ixm=17;
  int iX_ = iX+1;
  int iY_ = iY+1;
  nc=3;		  
  il=iym-iY_;
  ic=nc*(iX_-1);
  ii=il*nc*ixm+ic;
  if( str_.substr(ii,1).find("-")==string::npos ) return -1;
  return atoi( str_.substr(ii+1,2).c_str() );
}

int
MEEBGeom::lv_channel( EBTTLocalCoord iX, EBTTLocalCoord iY )
{
  string 
    str_
    (
     // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
     "-01-02-02-04-04-06-06-08-08-10-10-12-12-14-14-16-16" // 3
     "-01-02-02-04-04-06-06-08-08-10-10-12-12-14-14-16-16" // 2
     "-01-03-03-05-05-07-07-09-09-11-11-13-13-15-15-17-17" // 1
     "-01-03-03-05-05-07-07-09-09-11-11-13-13-15-15-17-17" // 0
     // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
     );

  int iym, ixm, nc, il, ic, ii;
  iym=4;
  ixm=17;
  int iX_ = iX+1;
  int iY_ = iY+1;
  nc=3;		  
  il=iym-iY_;
  ic=nc*(iX_-1);
  ii=il*nc*ixm+ic;
  if( str_.substr(ii,1).find("-")==string::npos ) return -1;
  return atoi( str_.substr(ii+1,2).c_str() );
}

int 
MEEBGeom::tt_channel( EBTTLocalCoord iX, EBTTLocalCoord iY )
{

  int itt =  4*iX+4-iY;  // :)

  return itt;

}

int 
MEEBGeom::crystal_channel( EBLocalCoord ix, EBLocalCoord iy )
{
  // "Test beam numbering"
  int icr =  20*ix + 19-iy + 1;  // :)

  return icr;
}

int 
MEEBGeom::electronic_channel( EBLocalCoord ix, EBLocalCoord iy )
{
  int iX = ix/5;
  int iY = iy/5;
  int itt  = tt_channel( iX, iY );
  int type = tt_type( iX, iY );

  int iVFE = ix%5+1;
  int islot = iy%5+1;
  if( iVFE%2==1 ) islot = 6-islot;
  int icr = 5*(iVFE-1)+(islot-1);
  // rotate for type-1 towers
  if( type==1 )
    {
      icr = 24-icr;
    }
  icr += 25*(itt-1);

  return icr;
}

TGraph*
MEEBGeom::getGraphBoundary(  int type, int num, bool global )
{
  int ism_=0;
  if( type==iSuperModule )
    {
      ism_ = num; 
    }
  else if( type==iLMRegion )
    {
      ism_ = (num-1)/2+1;
      if( ism_>18 ) 
	ism_-=18;
      else
	ism_+=18;
    }
  else
      abort();

  int ism=1;
  if( global ) ism = ism_;

  //  list< pair< float, float > > l;
  // getBoundary( l, type, num, global, ism );
  //  int n = l.size();
  //  if( n==0 ) return 0;
  
  float ix[100];
  float iy[100];
  int ixmin =  0;
  int ixmax = 84;
  int iymin =  0;
  int iymax = 19;

  int n(0);
  if( type==iSuperModule )
    {
      n=5;
      ix[0] = ixmin-0.5;
      iy[0] = iymin-0.5;
      ix[1] = ixmin-0.5;
      iy[1] = iymax+0.5;
      ix[2] = ixmax+0.5;
      iy[2] = iymax+0.5;
      ix[3] = ixmax+0.5;
      iy[3] = iymin-0.5;
      ix[4] = ix[0];
      iy[4] = iy[0];
    }
  else if( type==iLMRegion )
    {
      int iside = num;
      if( global )
	{
	  iside = (num-1)%2;
	}
      
      if( iside==1 )
	{
	  n=5;
	  ix[0] = ixmin+5-0.5;
	  iy[0] = iymin+10-0.5;
	  ix[1] = ixmin+5-0.5;
	  iy[1] = iymax+0.5;
	  ix[2] = ixmax+0.5;
	  iy[2] = iymax+0.5;
	  ix[3] = ixmax+0.5;
	  iy[3] = iymin+10-0.5;
	  ix[4] = ix[0];
	  iy[4] = iy[0];
	}
      else
	{
	  n=7;
	  ix[0] = ixmin-0.5;
	  iy[0] = iymin-0.5;
	  ix[1] = ixmin-0.5;
	  iy[1] = iymax+0.5;
	  ix[2] = ixmin+5-0.5;
	  iy[2] = iymax+0.5;
	  ix[3] = ixmin+5-0.5;
	  iy[3] = iymax-10+0.5;
	  ix[4] = ixmax+0.5;
	  iy[4] = iymax-10+0.5;
	  ix[5] = ixmax+0.5;
	  iy[5] = iymin-0.5;
	  ix[6] = ix[0];
	  iy[6] = iy[0];
	}
    }

  if( global )
    {
      for( int ii=0; ii<n; ii ++ )
	{
	  pair<float,float> xy = globalCoord( ism, ix[ii], iy[ii] ); 
	  ix[ii] = xy.first;
	  iy[ii] = xy.second;
	}
    }

//   int ii=0;
//   list< pair< float, float > >::const_iterator l_it;      
//   for( l_it=l.begin(); l_it!=l.end(); l_it++ )
//     {
//       //      cout << "[" << l_it->first << "," << l_it->second << "]" << endl;
//       ix[ii] = l_it->first;
//       iy[ii] = l_it->second;
//       ii++;
//     }


//  assert( ii==n );
  return new TGraph( n, ix, iy );
}

pair< int, int >
MEEBGeom::pn( int ilmmod )
{
  switch( ilmmod )
    {
    case   1: return pair<int,int>(  0,  5 );
    case   2: return pair<int,int>(  1,  6 );
    case   3: return pair<int,int>(  1,  6 );
    case   4: return pair<int,int>(  2,  7 );
    case   5: return pair<int,int>(  2,  7 );
    case   6: return pair<int,int>(  3,  8 );
    case   7: return pair<int,int>(  3,  8 );
    case   8: return pair<int,int>(  4,  9 );
    case   9: return pair<int,int>(  4,  9 );
    default:
      abort();
    }
  return pair<int,int>(-1,-1);
}

pair<int,int> 
MEEBGeom::memFromLmr( int ilmr )
{
  pair< int, int > dccAndSide_ = ME::dccAndSide( ilmr );
  int idcc  = dccAndSide_.first;
  return pair<int,int>( idcc, idcc );
}

vector<int> 
MEEBGeom::lmmodFromLmr( int ilmr )
{
  pair< int, int > dccAndSide_ = ME::dccAndSide( ilmr );
  int iside = dccAndSide_.second;
  vector< int > vec;
  for( int ilmmod=1; ilmmod<=9; ilmmod++ )
    {
      if( (ilmmod+iside)%2==1 )  vec.push_back(ilmmod);
    }
  return vec;
}


