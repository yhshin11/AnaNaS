#ifndef MEEEDisplay_hh
#define MEEEDisplay_hh

//
// Authors  : Gautier Hamel de Monchenault and Julie Malcles, Saclay 
//

#include <vector>
#include <map>
#include <list>

#include "MEEEGeom.hh"

class MEEEDisplay
{
  // static functions
public:

  static MEEEGeom::EtaPhiPoint getNode( MEEEGeom::SuperCrysCoord iX, 
					MEEEGeom::SuperCrysCoord iY, 
					int iz, int jx, int jy );
  static MEEEGeom::EtaPhiPoint center( int ix, int iy, int iz );
  
  static void drawEEGlobal();
  static void drawEELocal( int isect );

  static void drawEE( float etamin=-1000, float etamax=+1000, float phimin=-1000, float phimax=1000 );
  static void drawSC( MEEEGeom::SuperCrysCoord iX, 
		      MEEEGeom::SuperCrysCoord iY, 
		      int iz, float shift=0.,
		      float etamin=-1000, float etamax=+1000, float phimin=-1000, float phimax=1000 
		      );

  static void drawXtal( MEEEGeom::CrysCoord ix, 
			MEEEGeom::CrysCoord iy, 
			int iz, int color=kBlue, float shift=0,
			float etamin=-1000, float etamax=+1000, float phimin=-1000, float phimax=1000 );
  static TPolyLine* getXtalPolyLine( MEEEGeom::CrysCoord ix, 
				     MEEEGeom::CrysCoord iy, 
				     int iz, float shift=0,
				     float etamin=-1000, float etamax=+1000, float phimin=-1000, float phimax=1000 ); 
  static TPolyLine* getSCPolyLine( MEEEGeom::SuperCrysCoord iX, 
				   MEEEGeom::SuperCrysCoord iY, 
				   int iz, float shift=0,
				   float etamin=-1000, float etamax=+1000, float phimin=-1000, float phimax=1000 ); 
  static void drawRz( int lc=kGreen+3, int fc=kGreen-9 );

  static int bkgColor;
  static int lineColor;
  static int lineWidth;

  static void refresh();

  virtual ~MEEEDisplay() {}

private:

  static void sc_nodes( int itype, std::vector< int >& jx, std::vector< int >& jy );
  static void set( int iX, int iY, int jx, int jy, float eta, float phioverpi );
  static void setFirstQuadrant();

  static std::map< int, MEEEGeom::EtaPhiPoint > _pointMap;

  static void setRzXtals();
  static std::map< int, TPolyLine* > _rzXtals;

  static std::list<TObject*> _list;
  static void registerTObject( TObject* );

  ClassDef(MEEEDisplay,0) // MEEEDisplay -- Monitoring utility for survey of Ecal

    };

#endif

